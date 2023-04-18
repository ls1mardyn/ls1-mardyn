/*
 * ExtendedProfileSampling.cpp
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#include "ExtendedProfileSampling.h"

#include <math.h>
#include <limits>
#include <algorithm>
#include <map>

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "utils/Random.h"
#include "utils/FileUtils.h"
#include "Simulation.h"
#include "longRange/LongRangeCorrection.h"


ExtendedProfileSampling::ExtendedProfileSampling() {}

void ExtendedProfileSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    _numBinsGlobal = static_cast<unsigned int>(_globalBoxLength[1]/_binwidth);
    if (_globalBoxLength[1]/_binwidth != static_cast<float>(_numBinsGlobal)) {
        global_log->error() << "[ExtendedProfileSampling] Can not divide domain without remainder! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }
    _slabVolume = _globalBoxLength[0]*_globalBoxLength[2]*_binwidth;

    if (_slabVolume < 1e-12) {
        global_log->error() << "[ExtendedProfileSampling] Slab volume too small (<1e-12)!" << std::endl;
        Simulation::exit(-1);
    }

    _numComps = domain->getNumberOfComponents();

     // Entry per component and bin; 0 represents all components combined
    _lenVector = (_singleComp) ? _numBinsGlobal : _numBinsGlobal * (_numComps+1);

    resizeVectors();
    resetVectors();

    _cellProcessor = _simulation.getCellProcessor();
    _particlePairsHandler = std::make_shared<ParticlePairs2PotForceAdapter>(*domain);
    // MolID is maximum possible number minus rank to prevent duplicate IDs
    // Always insert molecule of first component
    const unsigned long molID = std::numeric_limits<unsigned long>::max() - static_cast<unsigned long>(domainDecomp->getRank());
    _mTest = Molecule(molID, &(_simulation.getEnsemble()->getComponents()->at(0)));
}

void ExtendedProfileSampling::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("binwidth", _binwidth);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("stop", _stopSampling);
    xmlconfig.getNodeValue("singlecomponent", _singleComp);
    xmlconfig.getNodeValue("highermoments", _sampleHigherMoms);

    xmlconfig.getNodeValue("chemicalpotential/@enable", _sampleChemPot);
    xmlconfig.getNodeValue("chemicalpotential/lattice", _lattice);
    xmlconfig.getNodeValue("chemicalpotential/factorNumTest", _factorNumTest);
    xmlconfig.getNodeValue("chemicalpotential/samplefrequency", _samplefrequency);

    std::string strCids;
    if (xmlconfig.getNodeValue("chemicalpotential/cids", strCids)) {
        // Parse string of comma-separated values into vector
        std::stringstream ssCids(strCids);
        for (unsigned int i; ssCids >> i;) {
            _cidsTest.push_back(i);
            if (ssCids.peek() == ',' || ssCids.peek() == ' ') {
                ssCids.ignore();
            } else if ((i <= 0) || (i > _numComps)) {
                global_log->warning() << "[ExtendedProfileSampling] cid " << i << " is not valid! cids must be between 1 and " << _numComps << std::endl;
            }
        }
    } else {
        _cidsTest.push_back(1);  // Default
    }
    for(auto const& cid : _cidsTest) { strCids += std::to_string(cid) + ", "; }
    strCids.erase(strCids.length()-2);  // Remove last ", "

    std::string insMethod;
    if (_lattice) {
        insMethod = "in a lattice";
    } else {
        insMethod = "randomly";
    }

    global_log->info() << "[ExtendedProfileSampling] Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    global_log->info() << "[ExtendedProfileSampling] Binwidth: " << _binwidth << std::endl;
    if (_singleComp) {
        global_log->info() << "[ExtendedProfileSampling] All components treated as single one" << std::endl;
    } else {
        global_log->info() << "[ExtendedProfileSampling] All components sampled individually" << std::endl;
    }
    if (_sampleHigherMoms) {
        global_log->info() << "[ExtendedProfileSampling] Sampling of higher moments enabled" << std::endl;
    } else {
        global_log->info() << "[ExtendedProfileSampling] Sampling of higher moments disabled" << std::endl;
    }
    if (_sampleChemPot) {
        global_log->info() << "[ExtendedProfileSampling] Sampling of chemical potential enabled with a sampling frequency of " << _samplefrequency << std::endl;
        global_log->info() << "[ExtendedProfileSampling] " << _factorNumTest << " * numParticles will be inserted " << insMethod << std::endl;
        global_log->info() << "[ExtendedProfileSampling] Inserting particles with cids = " << strCids << std::endl;
        if (_samplefrequency > _writeFrequency) {
            global_log->warning() << "[ExtendedProfileSampling] Sample frequency (" << _samplefrequency << ") is greater than write frequency ("
                                  << _writeFrequency << ")! " << std::endl;
        }
        if ((_cidsTest.size() > 1) and (_singleComp)) {
            global_log->warning() << "[ExtendedProfileSampling] <singlecomponent> and <cids> set! Output may not include all specified components" << std::endl;
        }
    } else {
        global_log->info() << "[ExtendedProfileSampling] Sampling of chemical potential disabled" << std::endl;
    }
}

void ExtendedProfileSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    // Sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep <= _startSampling) or (simstep > _stopSampling)) {
        return;
    }

    // Do not write or sample data directly after (re)start in the first step
    if ( simstep == _simulation.getNumInitTimesteps() ) {
        return;
    }
    
    // Variables per step
    std::array<double, 3> regionLowCorner {};
    std::array<double, 3> regionHighCorner {};
    std::array<double, 3> regionSize {};
    for (unsigned short d = 0; d < 3; d++) {
        regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
        regionSize[d] = regionHighCorner[d] - regionLowCorner[d];
    }


    CommVar<std::vector<unsigned long>> numMolecules_step;
    CommVar<std::vector<double>> mass_step;
    CommVar<std::vector<double>> ekin_step;                        // Including drift energy
    CommVar<std::vector<double>> epot_step;
    CommVar<std::vector<double>> orientation_step;
    CommVar<std::vector<double>> chemPot_step;
    CommVar<std::vector<unsigned long>> countNTest_step;
    std::array<CommVar<std::vector<double>>, 3> ekinVect_step;
    std::array<CommVar<std::vector<double>>, 3> velocityVect_step;
    std::array<CommVar<std::vector<double>>, 3> virialVect_step;
    std::array<CommVar<std::vector<double>>, 3> forceVect_step;
    std::array<CommVar<std::vector<double>>, 3> energyfluxVect_step;

    std::array<std::vector<double>, 3> veloDrift_step_global;       // Drift velocity
    std::vector<double> temperature_step_global (_lenVector, 0.0);  // Required for sampling of chem. pot.

    CommVar<std::vector<double>> hmDelta_step;
    std::array<CommVar<std::vector<double>>, 3> hmHeatflux_step;
    std::array<CommVar<std::vector<double>>, 9> hmPressure_step;
    std::array<CommVar<std::vector<double>>, 9> hmR_step;
    std::array<CommVar<std::vector<double>>, 27> hmM_step;

    if (_sampleHigherMoms) {

        hmDelta_step.local.resize(_lenVector);

        hmDelta_step.global.resize(_lenVector);

        std::fill(hmDelta_step.local.begin(), hmDelta_step.local.end(), 0.0f);

        std::fill(hmDelta_step.global.begin(), hmDelta_step.global.end(), 0.0f);

        for (unsigned short d = 0; d < 3; d++) {
            hmHeatflux_step.at(d).local.resize(_lenVector);
            hmHeatflux_step.at(d).global.resize(_lenVector);

            std::fill(hmHeatflux_step.at(d).local.begin(),  hmHeatflux_step.at(d).local.end(), 0.0f);
            std::fill(hmHeatflux_step.at(d).global.begin(), hmHeatflux_step.at(d).global.end(), 0.0f);
        }

        for (unsigned short d = 0; d < 9; d++) {
            hmPressure_step.at(d).local.resize(_lenVector);
            hmR_step.at(d).local.resize(_lenVector);
            hmM_step.at(d).local.resize(_lenVector);
            hmM_step.at(d+9u).local.resize(_lenVector);
            hmM_step.at(d+18u).local.resize(_lenVector);

            hmPressure_step.at(d).global.resize(_lenVector);
            hmR_step.at(d).global.resize(_lenVector);
            hmM_step.at(d).global.resize(_lenVector);
            hmM_step.at(d+9u).global.resize(_lenVector);
            hmM_step.at(d+18u).global.resize(_lenVector);

            std::fill(hmPressure_step.at(d).local.begin(), hmPressure_step.at(d).local.end(), 0.0f);
            std::fill(hmR_step.at(d).local.begin(),        hmR_step.at(d).local.end(), 0.0f);
            std::fill(hmM_step.at(d).local.begin(),        hmM_step.at(d).local.end(), 0.0f);
            std::fill(hmM_step.at(d+9u).local.begin(),     hmM_step.at(d+9u).local.end(), 0.0f);
            std::fill(hmM_step.at(d+18u).local.begin(),    hmM_step.at(d+18u).local.end(), 0.0f);

            std::fill(hmPressure_step.at(d).global.begin(), hmPressure_step.at(d).global.end(), 0.0f);
            std::fill(hmR_step.at(d).global.begin(),        hmR_step.at(d).global.end(), 0.0f);
            std::fill(hmM_step.at(d).global.begin(),        hmM_step.at(d).global.end(), 0.0f);
            std::fill(hmM_step.at(d+9u).global.begin(),     hmM_step.at(d+9u).global.end(), 0.0f);
            std::fill(hmM_step.at(d+18u).global.begin(),    hmM_step.at(d+18u).global.end(), 0.0f);
        }

    } else {
        hmDelta_step.local.clear();
        hmDelta_step.global.clear();
        for (unsigned short d = 0; d < 3; d++) {
            hmHeatflux_step.at(d).local.clear();
            hmHeatflux_step.at(d).global.clear();
        }
        for (unsigned short d = 0; d < 9; d++) {
            hmPressure_step.at(d).local.clear();
            hmR_step.at(d).local.clear();
            hmM_step.at(d).local.clear();
            hmM_step.at(d+9u).local.clear();
            hmM_step.at(d+18u).local.clear();
            hmPressure_step.at(d).global.clear();
            hmR_step.at(d).global.clear();
            hmM_step.at(d).global.clear();
            hmM_step.at(d+9u).global.clear();
            hmM_step.at(d+18u).global.clear();
        }
    }

    numMolecules_step.local.resize(_lenVector);
    mass_step.local.resize(_lenVector);
    ekin_step.local.resize(_lenVector);
    epot_step.local.resize(_lenVector);
    orientation_step.local.resize(_lenVector);
    chemPot_step.local.resize(_lenVector);
    countNTest_step.local.resize(_lenVector);

    numMolecules_step.global.resize(_lenVector);
    mass_step.global.resize(_lenVector);
    ekin_step.global.resize(_lenVector);
    epot_step.global.resize(_lenVector);
    orientation_step.global.resize(_lenVector);
    chemPot_step.global.resize(_lenVector);
    countNTest_step.global.resize(_lenVector);
    
    for (unsigned short d = 0; d < 3; d++) {
        ekinVect_step.at(d).local.resize(_lenVector);
        velocityVect_step.at(d).local.resize(_lenVector);
        virialVect_step.at(d).local.resize(_lenVector);
        forceVect_step.at(d).local.resize(_lenVector);
        energyfluxVect_step.at(d).local.resize(_lenVector);

        ekinVect_step.at(d).global.resize(_lenVector);
        velocityVect_step.at(d).global.resize(_lenVector);
        virialVect_step.at(d).global.resize(_lenVector);
        forceVect_step.at(d).global.resize(_lenVector);
        energyfluxVect_step.at(d).global.resize(_lenVector);

        veloDrift_step_global.at(d).resize(_lenVector);
    }

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);
    std::fill(mass_step.local.begin(),        mass_step.local.end(), 0.0f);
    std::fill(ekin_step.local.begin(),        ekin_step.local.end(), 0.0f);
    std::fill(epot_step.local.begin(),         epot_step.local.end(), 0.0f);
    std::fill(orientation_step.local.begin(),  orientation_step.local.end(), 0.0f);
    std::fill(chemPot_step.local.begin(),      chemPot_step.local.end(), 0.0f);
    std::fill(countNTest_step.local.begin(),   countNTest_step.local.end(), 0ul);

    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);
    std::fill(mass_step.global.begin(),        mass_step.global.end(), 0.0f);
    std::fill(ekin_step.global.begin(),        ekin_step.global.end(), 0.0f);
    std::fill(epot_step.global.begin(),         epot_step.global.end(), 0.0f);
    std::fill(orientation_step.global.begin(),  orientation_step.global.end(), 0.0f);
    std::fill(chemPot_step.global.begin(),      chemPot_step.global.end(), 0.0f);
    std::fill(countNTest_step.global.begin(),   countNTest_step.global.end(), 0ul);
    
    for (unsigned short d = 0; d < 3; d++) {
        std::fill(ekinVect_step.at(d).local.begin(),      ekinVect_step.at(d).local.end(), 0.0f);
        std::fill(velocityVect_step.at(d).local.begin(),   velocityVect_step.at(d).local.end(), 0.0f);
        std::fill(virialVect_step.at(d).local.begin(),     virialVect_step.at(d).local.end(), 0.0f);
        std::fill(forceVect_step.at(d).local.begin(),      forceVect_step.at(d).local.end(), 0.0f);
        std::fill(energyfluxVect_step.at(d).local.begin(), energyfluxVect_step.at(d).local.end(), 0.0f);

        std::fill(ekinVect_step.at(d).global.begin(),      ekinVect_step.at(d).global.end(), 0.0f);
        std::fill(velocityVect_step.at(d).global.begin(),   velocityVect_step.at(d).global.end(), 0.0f);
        std::fill(virialVect_step.at(d).global.begin(),     virialVect_step.at(d).global.end(), 0.0f);
        std::fill(forceVect_step.at(d).global.begin(),      forceVect_step.at(d).global.end(), 0.0f);
        std::fill(energyfluxVect_step.at(d).global.begin(), energyfluxVect_step.at(d).global.end(), 0.0f);

        std::fill(veloDrift_step_global.at(d).begin(),      veloDrift_step_global.at(d).end(), 0.0f);
    }

    // Calculate drift as it is needed first
    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double ry = pit->r(1);
        const unsigned int index = std::min(_numBinsGlobal, static_cast<unsigned int>(ry/_binwidth));  // Index of bin

        std::vector<unsigned int> cids = {0};  // add velocities to "all components" (0) and respective component
        if (!_singleComp) { cids.push_back(pit->componentid() + 1); }

        for (unsigned int cid : cids) {
            const uint32_t indexCID = cid*_numBinsGlobal + index;
            numMolecules_step.local.at(indexCID) ++;
            velocityVect_step.at(0).local.at(indexCID) += pit->v(0);
            velocityVect_step.at(1).local.at(indexCID) += pit->v(1);
            velocityVect_step.at(2).local.at(indexCID) += pit->v(2);
        }
    }

#ifdef ENABLE_MPI
    MPI_Allreduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[0].local.data(), velocityVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[1].local.data(), velocityVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[2].local.data(), velocityVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global.at(i) = numMolecules_step.local.at(i);
        velocityVect_step[0].global.at(i) = velocityVect_step[0].local.at(i);
        velocityVect_step[1].global.at(i) = velocityVect_step[1].local.at(i);
        velocityVect_step[2].global.at(i) = velocityVect_step[2].local.at(i);
    }
#endif

    for (unsigned long i = 0; i < _lenVector; i++) {
        if (numMolecules_step.global.at(i) > 0ul) {
            veloDrift_step_global[0].at(i) = velocityVect_step[0].global.at(i) / numMolecules_step.global.at(i);
            veloDrift_step_global[1].at(i) = velocityVect_step[1].global.at(i) / numMolecules_step.global.at(i);
            veloDrift_step_global[2].at(i) = velocityVect_step[2].global.at(i) / numMolecules_step.global.at(i);
        }
    }

    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double ry = pit->r(1);
        const unsigned int index = std::min(_numBinsGlobal, static_cast<unsigned int>(ry/_binwidth));  // Index of bin

        const double veloX = pit->v(0);
        const double veloY = pit->v(1);
        const double veloZ = pit->v(2);
        const double mass = pit->mass();
        const double epot = pit->U_pot();

        std::vector<unsigned int> cids = {0};  // add quantities to "all components" (0) and respective component
        if (!_singleComp) { cids.push_back(pit->componentid() + 1); }

        for (unsigned int cid : cids) {
            const uint32_t indexCID = cid*_numBinsGlobal + index;
            mass_step.local.at(indexCID) += pit->mass();
            ekin_step.local.at(indexCID) += pit->U_kin();
            epot_step.local.at(indexCID) += epot;

            Quaternion q_pit = pit->q();
            // Molecule should be elongated (ie. main axis) in z direction
            std::array<double, 3> pos_rot = q_pit.rotate({0., 0., 1.});
            // Orientation order parameter, cf. Eq. 29 in Mecke2001
            // cos^2(theta) = dy^2 = pos_rot[1]^2 since length of rotated vector (hypotenuse) is 1
            orientation_step.local.at(indexCID) += 3.*pos_rot[1]*pos_rot[1] - 1.;
            // Angle between normal vector of y plane and z axis of molecule
            // const double angle = 90 - (std::atan(abs(pos_rot[1])/(std::sqrt(pos_rot[0]*pos_rot[0]+pos_rot[2]*pos_rot[2])))*(180/3.14159265));

            ekinVect_step[0].local.at(indexCID) += 0.5*mass*veloX*veloX;
            ekinVect_step[1].local.at(indexCID) += 0.5*mass*veloY*veloY;
            ekinVect_step[2].local.at(indexCID) += 0.5*mass*veloZ*veloZ;
            virialVect_step[0].local.at(indexCID) += pit->Vi(0);
            virialVect_step[1].local.at(indexCID) += pit->Vi(1);
            virialVect_step[2].local.at(indexCID) += pit->Vi(2);
            forceVect_step[0].local.at(indexCID) += pit->F(0);
            forceVect_step[1].local.at(indexCID) += pit->F(1);
            forceVect_step[2].local.at(indexCID) += pit->F(2);
            energyfluxVect_step[0].local.at(indexCID) += (pit->U_kin() + epot)*veloX + (pit->Vi(0)*veloX + pit->Vi(3)*veloY + pit->Vi(4)*veloZ);
            energyfluxVect_step[1].local.at(indexCID) += (pit->U_kin() + epot)*veloY + (pit->Vi(6)*veloX + pit->Vi(1)*veloY + pit->Vi(5)*veloZ);
            energyfluxVect_step[2].local.at(indexCID) += (pit->U_kin() + epot)*veloZ + (pit->Vi(7)*veloX + pit->Vi(8)*veloY + pit->Vi(2)*veloZ);

            if (_sampleHigherMoms) {
                const double veloCorrX = veloX - veloDrift_step_global[0].at(index);
                const double veloCorrY = veloY - veloDrift_step_global[1].at(index);
                const double veloCorrZ = veloZ - veloDrift_step_global[2].at(index);
                const double veloCorrSqrt = veloCorrX*veloCorrX + veloCorrY*veloCorrY + veloCorrZ*veloCorrZ;  // Squared velocity of particle without drift
                const std::array<double, 3> velos = {veloCorrX, veloCorrY, veloCorrZ};

                hmDelta_step.local.at(indexCID) += veloCorrSqrt*veloCorrSqrt;

                for (unsigned short i = 0; i < 3; i++) {

                    hmHeatflux_step[i].local.at(indexCID)     += 0.5*veloCorrSqrt*velos[i];

                    for (unsigned short j = 0; j < 3; j++) {
                        if (i == j) {  // Trace elements
                            hmPressure_step[3*i+j].local.at(indexCID) += velos[i]*velos[j] - (1./3.)*veloCorrSqrt;  // Pressure; cxcx, cxcy, cxcz, cycx, cycy, cycz, czcx, czcy, czcz
                            hmR_step[3u*i+j].local.at(indexCID)        += (velos[i]*velos[j] - (1./3.)*veloCorrSqrt)*veloCorrSqrt;   // R; cxcx, cxcy, cxcz, cycx, cycy, cycz, czcx, czcy, czcz
                        } else {
                            hmPressure_step[3u*i+j].local.at(indexCID) += velos[i]*velos[j];
                            hmR_step[3u*i+j].local.at(indexCID)        += velos[i]*velos[j]*veloCorrSqrt;
                        }
                    }
                }

                double m[3][3][3] = {0.0};
                for (unsigned short i = 0; i < 3; i++) {
                    for (unsigned short j = 0; j < 3; j++) {
                        for (unsigned short k = 0; k < 3; k++) {
                            m[i][j][k] = velos[i]*velos[j]*velos[k];
                        }
                    }
                }

                const double mSum[3] = {
                    m[0][0][0] + m[0][1][1] + m[0][2][2],
                    m[1][0][0] + m[1][1][1] + m[1][2][2],
                    m[2][0][0] + m[2][1][1] + m[2][2][2]
                };

                hmM_step[0].local.at(indexCID) += m[0][0][0] - 0.6*(mSum[0]);  // m: cxcxcx
                hmM_step[1].local.at(indexCID) += m[0][0][1] - 0.2*(mSum[1]);  // m: cxcxcy
                hmM_step[2].local.at(indexCID) += m[0][0][2] - 0.2*(mSum[2]);  // m: cxcxcz
                hmM_step[3].local.at(indexCID) += m[0][1][0] - 0.2*(mSum[1]);  // m: cxcycx
                hmM_step[4].local.at(indexCID) += m[0][1][1] - 0.2*(mSum[0]);  // m: cxcycy
                hmM_step[5].local.at(indexCID) += m[0][1][2];                  // m: cxcycz
                hmM_step[6].local.at(indexCID) += m[0][2][0] - 0.2*(mSum[2]);  // m: cxczcx
                hmM_step[7].local.at(indexCID) += m[0][2][1];                  // m: cxczcy
                hmM_step[8].local.at(indexCID) += m[0][2][2] - 0.2*(mSum[0]);  // m: cxczcz

                hmM_step[9].local.at(indexCID) += m[1][0][0] - 0.2*(mSum[1]);  // m: cycxcx
                hmM_step[10].local.at(indexCID) += m[1][0][1] - 0.2*(mSum[0]); // m: cycxcy
                hmM_step[11].local.at(indexCID) += m[1][0][2];                 // m: cycxcz
                hmM_step[12].local.at(indexCID) += m[1][1][0] - 0.2*(mSum[0]); // m: cycycx
                hmM_step[13].local.at(indexCID) += m[1][1][1] - 0.6*(mSum[1]); // m: cycycy
                hmM_step[14].local.at(indexCID) += m[1][1][2] - 0.2*(mSum[2]); // m: cycycz
                hmM_step[15].local.at(indexCID) += m[1][2][0];                 // m: cyczcx
                hmM_step[16].local.at(indexCID) += m[1][2][1] - 0.2*(mSum[2]); // m: cyczcy
                hmM_step[17].local.at(indexCID) += m[1][2][2] - 0.2*(mSum[1]); // m: cyczcz

                hmM_step[18].local.at(indexCID) += m[2][0][0] - 0.2*(mSum[2]); // m: czcxcx
                hmM_step[19].local.at(indexCID) += m[2][0][1];                 // m: czcxcy
                hmM_step[20].local.at(indexCID) += m[2][0][2] - 0.2*(mSum[0]); // m: czcxcz
                hmM_step[21].local.at(indexCID) += m[2][1][0];                 // m: czcycx
                hmM_step[22].local.at(indexCID) += m[2][1][1] - 0.2*(mSum[2]); // m: czcycy
                hmM_step[23].local.at(indexCID) += m[2][1][2] - 0.2*(mSum[1]); // m: czcycz
                hmM_step[24].local.at(indexCID) += m[2][2][0] - 0.2*(mSum[0]); // m: czczcx
                hmM_step[25].local.at(indexCID) += m[2][2][1] - 0.2*(mSum[1]); // m: czczcy
                hmM_step[26].local.at(indexCID) += m[2][2][2] - 0.6*(mSum[2]); // m: czczcz
            }
        }
    }

    // Calculate chemical potential
    // if sampling of chem. pot. is activated
    if (_sampleChemPot and
        // only conducted every _samplefrequency step
        ((simstep - _startSampling) % _samplefrequency == 0 ) and
        // only in second half of sampling interval in order to be able to calculate a correct temperature
        ((simstep - _startSampling) % _writeFrequency > 0.5*_writeFrequency)) {

        // Root calculates temperature (without drift) per bin over all processes in second half of sampling interval for chem. pot. sampling
        // and propagates it to other processes
        if (domainDecomp->getRank() == 0) {
            for (unsigned long index = 0; index < _lenVector; index++) {
                const double v_x = _velocityVect_accum[0].at(index) / _countSamples.at(index);
                const double v_y = _velocityVect_accum[1].at(index) / _countSamples.at(index);
                const double v_z = _velocityVect_accum[2].at(index) / _countSamples.at(index);
                double v_drift_sqr = v_x*v_x + v_y*v_y + v_z*v_z;
                temperature_step_global.at(index) = (2*_ekin_accum.at(index) - v_drift_sqr*_mass_accum.at(index)) / _doftotal_accum.at(index);
            }
        }
#ifdef ENABLE_MPI
        MPI_Bcast(temperature_step_global.data(), _lenVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        // Calculate number of test particles per bin
        std::vector<double> dX(_numBinsGlobal, 0.0);
        std::vector<double> dY(_numBinsGlobal, 0.0);
        std::vector<double> dZ(_numBinsGlobal, 0.0);
        std::vector<unsigned long> nX(_numBinsGlobal, 0ul);
        std::vector<unsigned long> nY(_numBinsGlobal, 0ul);
        std::vector<unsigned long> nZ(_numBinsGlobal, 0ul);
        unsigned long nTestGlobal {0ul};

        for (unsigned int i = 0; i < _numBinsGlobal; i++) {
            // Make sure, number of test particles is never zero and number of test particles per direction is at least 1
            const unsigned long nTest = std::max(1ul,static_cast<unsigned long>(_factorNumTest*numMolecules_step.global.at(i)));

            nY.at(i) = std::max(1.0,std::pow((nTest*_binwidth*_binwidth)/(_globalBoxLength[0]*_globalBoxLength[2]),(1./3.)));
            dY.at(i) = std::min(static_cast<float>(_binwidth/nY.at(i)), static_cast<float>(0.5*regionSize[1]));

            nX.at(i) = std::max(1.0,std::pow((nTest*_globalBoxLength[0]*_globalBoxLength[0])/(_binwidth*_globalBoxLength[2]),(1./3.)));
            dX.at(i) = std::min(static_cast<float>(_globalBoxLength[0]/nX.at(i)), static_cast<float>(0.5*regionSize[0]));

            nZ.at(i) = std::max(1ul,nTest/(nX.at(i)*nY.at(i)));
            dZ.at(i) = std::min(static_cast<float>(_globalBoxLength[2]/nZ.at(i)), static_cast<float>(0.5*regionSize[2]));
            
            nTestGlobal += nTest;
        }

        if (_lattice) {
            // Insert particles in lattice structure and sample chem. pot.

            // Index of bin in which the left region boundary (y-dir) is in; std::min if particle position is precisely at right boundary
            const unsigned int idxStart = std::min(_numBinsGlobal, static_cast<unsigned int>(regionLowCorner[1]/_binwidth));
            double rY = regionLowCorner[1]+0.5*dY.at(idxStart);
            while (rY < regionHighCorner[1]) {
                const unsigned int index = std::min(_numBinsGlobal, static_cast<unsigned int>(rY/_binwidth));  // Index of bin

                double rX = regionLowCorner[0]+0.5*dX.at(idxStart);
                while (rX < regionHighCorner[0]) {

                    double rZ = regionLowCorner[2]+0.5*dZ.at(idxStart);
                    while (rZ < regionHighCorner[2]) {
                        if (temperature_step_global.at(index) > 1e-9) {
                            _mTest.setr(0,rX);
                            _mTest.setr(1,rY);
                            _mTest.setr(2,rZ);
                            for (auto cid : _cidsTest) {
                                const uint32_t indexCID = cid*_numBinsGlobal + index;
                                _mTest.setComponent(&(_simulation.getEnsemble()->getComponents()->at(cid-1)));
                                const double deltaUpot = particleContainer->getEnergy(_particlePairsHandler.get(), &_mTest, *_cellProcessor)
                                                        + 2.0*_simulation.getLongRangeCorrection()->getUpotCorr(&_mTest);
                                double chemPot = exp(-deltaUpot/temperature_step_global.at(index));  // Global temperature of all components
                                if (std::isfinite(chemPot)) {
                                    // For component 0 (single component)
                                    chemPot_step.local.at(index) += chemPot;
                                    countNTest_step.local.at(index)++;
                                    // Also add to specific component record if no single-component-sampling
                                    if (!_singleComp) {
                                        chemPot_step.local.at(indexCID) += chemPot;
                                        countNTest_step.local.at(indexCID)++;
                                    }
#ifndef NDEBUG
                                    std::cout << "[ExtendedProfileSampling] Rank " << domainDecomp->getRank() << " : Inserting molecule at x,y,z = "
                                            << _mTest.r(0) << " , " << _mTest.r(1) << " , " << _mTest.r(2)
                                            << " ; cid = " << _mTest.componentid()
                                            << " ; chemPot = " << chemPot << " ; dU = " << deltaUpot << " ; T = " << temperature_step_global.at(index) << " ; index = " << index << std::endl;
#endif
                                }
                            }
                        }
                        rZ += dZ.at(index);
                    }
                    rX += dX.at(index);
                }
                rY += dY.at(index);
            }
        } else {
            // Random insertion
            // NOTE: This differs from the lattice method as it does not take the local density into account
            std::unique_ptr<Random> rnd(new Random());

            // Share of volume of present rank from whole domain
            const float domainShare = (regionSize[0]*regionSize[1]*regionSize[2])/(_globalBoxLength[0]*_globalBoxLength[1]*_globalBoxLength[2]); 
            const unsigned long nTest = static_cast<unsigned long>(domainShare*nTestGlobal);

            for (unsigned long i = 0; i < nTest; i++) {
                const double rX = regionLowCorner[0] + rnd->rnd()*regionSize[0];
                const double rY = regionLowCorner[1] + rnd->rnd()*regionSize[1];
                const double rZ = regionLowCorner[2] + rnd->rnd()*regionSize[2];
                const unsigned int index = std::min(_numBinsGlobal, static_cast<unsigned int>(rY/_binwidth));  // Index of bin
                if (temperature_step_global.at(index) > 1e-9) {
                    _mTest.setr(0,rX);
                    _mTest.setr(1,rY);
                    _mTest.setr(2,rZ);
                    for (auto cid : _cidsTest) {
                        const uint32_t indexCID = cid*_numBinsGlobal + index;
                        _mTest.setComponent(&(_simulation.getEnsemble()->getComponents()->at(cid-1)));
                        const double deltaUpot = particleContainer->getEnergy(_particlePairsHandler.get(), &_mTest, *_cellProcessor)
                                                + _simulation.getLongRangeCorrection()->getUpotCorr(&_mTest);
                        double chemPot = exp(-deltaUpot/temperature_step_global.at(index));
                        if (std::isfinite(chemPot)) {
                            // For component 0 (single component)
                            chemPot_step.local.at(index) += chemPot;
                            countNTest_step.local.at(index)++;
                            // Also add to specific component record if no single-component-sampling
                            if (!_singleComp) {
                                chemPot_step.local.at(indexCID) += chemPot;
                                countNTest_step.local.at(indexCID)++;
                            }
#ifndef NDEBUG
                            std::cout << "[ExtendedProfileSampling] Rank " << domainDecomp->getRank() << " : Inserting molecule at x,y,z = "
                                    << _mTest.r(0) << " , " << _mTest.r(1) << " , " << _mTest.r(2)
                                    << " ; cid = " << _mTest.componentid()
                                    << " ; chemPot = " << chemPot << " ; dU = " << deltaUpot << " ; T = " << temperature_step_global.at(index) << " ; index = " << index << std::endl;
#endif
                        }
                    }
                }
            }
        }
    }

    // Gather quantities of all processes. Note: MPI_Reduce instead of MPI_Allreduce!
#ifdef ENABLE_MPI
    MPI_Reduce(mass_step.local.data(), mass_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin_step.local.data(), ekin_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(epot_step.local.data(), epot_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(orientation_step.local.data(), orientation_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[0].local.data(), virialVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[1].local.data(), virialVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[2].local.data(), virialVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(forceVect_step[0].local.data(), forceVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(forceVect_step[1].local.data(), forceVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(forceVect_step[2].local.data(), forceVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[0].local.data(), ekinVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[1].local.data(), ekinVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[2].local.data(), ekinVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(energyfluxVect_step[0].local.data(), energyfluxVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(energyfluxVect_step[1].local.data(), energyfluxVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(energyfluxVect_step[2].local.data(), energyfluxVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(chemPot_step.local.data(), chemPot_step.global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(countNTest_step.local.data(), countNTest_step.global.data(), _numBinsGlobal, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (_sampleHigherMoms) {
        MPI_Reduce(hmDelta_step.local.data(), hmDelta_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        for (unsigned short d = 0; d < 3; d++) {
            MPI_Reduce(hmHeatflux_step[d].local.data(), hmHeatflux_step[d].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        for (unsigned short d = 0; d < 9; d++) {
            MPI_Reduce(hmPressure_step[d].local.data(), hmPressure_step[d].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(hmR_step[d].local.data(), hmR_step[d].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(hmM_step[d].local.data(), hmM_step[d].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(hmM_step[d+9u].local.data(), hmM_step[d+9u].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(hmM_step[d+18u].local.data(), hmM_step[d+18u].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        mass_step.global.at(i) = mass_step.local.at(i);
        ekin_step.global.at(i) = ekin_step.local.at(i);
        epot_step.global.at(i) = epot_step.local.at(i);
        orientation_step.global.at(i) = orientation_step.local.at(i);
        virialVect_step[0].global.at(i) = virialVect_step[0].local.at(i);
        virialVect_step[1].global.at(i) = virialVect_step[1].local.at(i);
        virialVect_step[2].global.at(i) = virialVect_step[2].local.at(i);
        forceVect_step[0].global.at(i) = forceVect_step[0].local.at(i);
        forceVect_step[1].global.at(i) = forceVect_step[1].local.at(i);
        forceVect_step[2].global.at(i) = forceVect_step[2].local.at(i);
        ekinVect_step[0].global.at(i) = ekinVect_step[0].local.at(i);
        ekinVect_step[1].global.at(i) = ekinVect_step[1].local.at(i);
        ekinVect_step[2].global.at(i) = ekinVect_step[2].local.at(i);
        energyfluxVect_step[0].global.at(i) = energyfluxVect_step[0].local.at(i);
        energyfluxVect_step[1].global.at(i) = energyfluxVect_step[1].local.at(i);
        energyfluxVect_step[2].global.at(i) = energyfluxVect_step[2].local.at(i);
        chemPot_step.global.at(i) = chemPot_step.local.at(i);
        countNTest_step.global.at(i) = countNTest_step.local.at(i);
        if (_sampleHigherMoms) {
            hmDelta_step.global.at(i) = hmDelta_step.local.at(i);
            for (unsigned short d = 0; d < 3; d++) {
                hmHeatflux_step[d].global.at(i) = hmHeatflux_step[d].local.at(i);
            }
            for (unsigned short d = 0; d < 9; d++) {
                hmPressure_step[d].global.at(i) = hmPressure_step[d].local.at(i);
                hmR_step[d].global.at(i)    = hmR_step[d].local.at(i);
                hmM_step[d].global.at(i)    = hmM_step[d].local.at(i);
                hmM_step[d+9u].global.at(i)  = hmM_step[d+9u].local.at(i);
                hmM_step[d+18u].global.at(i) = hmM_step[d+18u].local.at(i);
            }
        }
    }
#endif

    // Only root knows real quantities (MPI_Reduce instead of MPI_Allreduce)
    // Accumulate data
    if (domainDecomp->getRank() == 0) {
        for (unsigned long i = 0; i < _lenVector; i++) {
            const unsigned long numMols = numMolecules_step.global.at(i);
            const unsigned int cid = i/_numBinsGlobal;
            unsigned int dof_rot {0};
            unsigned int dof_total {0};
            if (!_singleComp) {
                if (cid == 0) {
                    for (unsigned long cj = 0; cj < _numComps; cj++) {
                        dof_rot = _simulation.getEnsemble()->getComponent(cj)->getRotationalDegreesOfFreedom();
                        dof_total += (3 + dof_rot)*numMolecules_step.global.at((cj+1)*_numBinsGlobal + i);
                    }
                } else {
                    dof_rot = _simulation.getEnsemble()->getComponent(cid-1)->getRotationalDegreesOfFreedom();
                    dof_total = (3 + dof_rot)*numMols;
                }
            } else {
                // For single component sampling, the rot. DOF of component 0 is taken
                dof_rot = _simulation.getEnsemble()->getComponent(0)->getRotationalDegreesOfFreedom();
                dof_total = (3 + dof_rot)*numMols;
            }
            
            const double ViX = virialVect_step[0].global.at(i);
            const double ViY = virialVect_step[1].global.at(i);
            const double ViZ = virialVect_step[2].global.at(i);

            _doftotal_accum.at(i)                += dof_total;
            _numMolecules_accum.at(i)            += numMols;
            _mass_accum.at(i)                    += mass_step.global.at(i);
            _ekin_accum.at(i)                    += ekin_step.global.at(i);
            _epot_accum.at(i)                    += epot_step.global.at(i);
            _orientation_accum.at(i)             += orientation_step.global.at(i);
            _virial_accum.at(i)                  += ViX + ViY + ViZ;
            _chemPot_accum.at(i)                 += chemPot_step.global.at(i);
            _countNTest_accum.at(i)              += countNTest_step.global.at(i);

            _ekinVect_accum[0].at(i)             += ekinVect_step[0].global.at(i);
            _ekinVect_accum[1].at(i)             += ekinVect_step[1].global.at(i);
            _ekinVect_accum[2].at(i)             += ekinVect_step[2].global.at(i);
            _velocityVect_accum[0].at(i)         += veloDrift_step_global[0].at(i);
            _velocityVect_accum[1].at(i)         += veloDrift_step_global[1].at(i);
            _velocityVect_accum[2].at(i)         += veloDrift_step_global[2].at(i);
            _virialVect_accum[0].at(i)           += ViX;
            _virialVect_accum[1].at(i)           += ViY;
            _virialVect_accum[2].at(i)           += ViZ;
            _forceVect_accum[0].at(i)            += forceVect_step[0].global.at(i);
            _forceVect_accum[1].at(i)            += forceVect_step[1].global.at(i);
            _forceVect_accum[2].at(i)            += forceVect_step[2].global.at(i);
            _energyfluxVect_accum[0].at(i)       += energyfluxVect_step[0].global.at(i) / _slabVolume;
            _energyfluxVect_accum[1].at(i)       += energyfluxVect_step[1].global.at(i) / _slabVolume;
            _energyfluxVect_accum[2].at(i)       += energyfluxVect_step[2].global.at(i) / _slabVolume;

            _countSamples.at(i)++;

            if (_sampleHigherMoms) {
                _hmDelta_accum.at(i)             += hmDelta_step.global.at(i) / _slabVolume;
                for (unsigned short d = 0; d < 3; d++) {
                    _hmHeatflux_accum[d].at(i)   += hmHeatflux_step[d].global.at(i) / _slabVolume;
                }
                for (unsigned short d = 0; d < 9; d++) {
                    _hmPressure_accum[d].at(i)   += hmPressure_step[d].global.at(i) / _slabVolume;
                    _hmR_accum[d].at(i)          += hmR_step[d].global.at(i) / _slabVolume;
                    _hmM_accum[d].at(i)          += hmM_step[d].global.at(i) / _slabVolume;
                    _hmM_accum[d+9u].at(i)        += hmM_step[d+9u].global.at(i) / _slabVolume;
                    _hmM_accum[d+18u].at(i)       += hmM_step[d+18u].global.at(i) / _slabVolume;
                }
            }
        }
    }

    // Write out data every _writeFrequency step
    if ( (simstep - _startSampling) % _writeFrequency == 0 ) {

        if (domainDecomp->getRank() == 0) {
            unsigned long numOutputs = (_singleComp) ? 1ul : (_numComps+1);
            {
            // Write output file
            std::stringstream ss;
            ss << std::setw(9) << std::setfill('0') << simstep;
            const std::string fname = "ExtendedProfileSampling_TS"+ss.str()+".dat";
            std::ofstream ofs;
            ofs.open(fname, std::ios::out);
            ofs << setw(24) << "pos";                           // Bin position
            for (unsigned long cid = 0; cid < numOutputs; cid++) {
                ofs << setw(22) << "numParts[" << cid << "]"        // Average number of molecules in bin per step
                    << setw(22) << "rho["      << cid << "]"        // Density
                    << setw(22) << "T["        << cid << "]"        // Temperature without drift (i.e. "real" temperature)
                    << setw(22) << "ekin["     << cid << "]"        // Kinetic energy including drift
                    << setw(22) << "epot["     << cid << "]"        // Potential energy
                    << setw(22) << "orientation[" << cid << "]"     // Orientation order parameter, cf. Eq. 29 in Mecke2001
                    << setw(22) << "p["        << cid << "]"        // Pressure
                    << setw(22) << "chemPot_res[" << cid << "]"     // Chemical potential as known as mu_tilde (equals the ms2 value)
                    << setw(22) << "numTest["     << cid << "]"     // Number of inserted test particles per sample step
                    << setw(22) << "T_x["      << cid << "]"        // Temperature in x-direction
                    << setw(22) << "T_y["      << cid << "]"        // Temperature in y-direction
                    << setw(22) << "T_z["      << cid << "]"        // Temperature in z-direction
                    << setw(22) << "v_x["      << cid << "]"        // Drift velocity in x-direction
                    << setw(22) << "v_y["      << cid << "]"        // Drift velocity in y-direction
                    << setw(22) << "v_z["      << cid << "]"        // Drift velocity in z-direction
                    << setw(22) << "p_x["      << cid << "]"        // Pressure in x-direction
                    << setw(22) << "p_y["      << cid << "]"        // Pressure in y-direction
                    << setw(22) << "p_z["      << cid << "]"        // Pressure in z-direction
                    << setw(22) << "F_x["      << cid << "]"        // Force in x-direction
                    << setw(22) << "F_y["      << cid << "]"        // Force in y-direction
                    << setw(22) << "F_z["      << cid << "]"        // Force in z-direction
                    << setw(22) << "jEF_x["    << cid << "]"        // Energy flux in x-direction
                    << setw(22) << "jEF_y["    << cid << "]"        // Energy flux in y-direction
                    << setw(22) << "jEF_z["    << cid << "]"        // Energy flux in z-direction
                    << setw(22) << "numSamples["  << cid << "]";    // Number of samples (<= _writeFrequency)
            }
            ofs << std::endl;

            for (unsigned long idx = 0; idx < _numBinsGlobal; idx++) {
                ofs << FORMAT_SCI_MAX_DIGITS << (idx+0.5)*_binwidth;
                for (unsigned long cid = 0; cid < numOutputs; cid++) {
                    unsigned long i = idx + cid*_numBinsGlobal;
                    double numMolsPerStep {std::nan("0")}; // Not an int as particles change bin during simulation
                    double rho {std::nan("0")};
                    double T {std::nan("0")};
                    double ekin {std::nan("0")};
                    double epot {std::nan("0")};
                    double orderparam {std::nan("0")};
                    double p {std::nan("0")};
                    double chemPot_res {std::nan("0")};
                    double numTest {std::nan("0")};
                    double T_x {std::nan("0")};
                    double T_y {std::nan("0")};
                    double T_z {std::nan("0")};
                    double v_x {std::nan("0")};
                    double v_y {std::nan("0")};
                    double v_z {std::nan("0")};
                    double p_x {std::nan("0")};
                    double p_y {std::nan("0")};
                    double p_z {std::nan("0")};
                    double F_x {std::nan("0")};
                    double F_y {std::nan("0")};
                    double F_z {std::nan("0")};
                    double jEF_x {std::nan("0")};
                    double jEF_y {std::nan("0")};
                    double jEF_z {std::nan("0")};
                    double numSamples {std::nan("0")};
                    if ((_countSamples.at(i) > 0ul) and (_doftotal_accum.at(i)) > 0ul) {
                        const unsigned long countSamples = _countSamples.at(i);
                        const double numMols_accum = static_cast<double>(_numMolecules_accum.at(i));

                        numMolsPerStep = numMols_accum/countSamples;
                        rho         = numMolsPerStep           / _slabVolume;
                        v_x         = _velocityVect_accum[0].at(i)   / countSamples;
                        v_y         = _velocityVect_accum[1].at(i)   / countSamples;
                        v_z         = _velocityVect_accum[2].at(i)   / countSamples;

                        double v_drift_sqr = v_x*v_x + v_y*v_y + v_z*v_z;

                        T           = (2*_ekin_accum.at(i) - v_drift_sqr*_mass_accum.at(i)) / _doftotal_accum.at(i);
                        ekin        = _ekin_accum.at(i) / numMols_accum;
                        epot        = _epot_accum.at(i) / numMols_accum;
                        orderparam  = 0.5*_orientation_accum.at(i) / numMols_accum;;
                        p           = rho * ( (_virial_accum.at(i))/(3.0*numMols_accum) + T);

                        T_x         = (2*_ekinVect_accum[0].at(i) - (v_x*v_x)*_mass_accum.at(i)) / numMols_accum;
                        T_y         = (2*_ekinVect_accum[1].at(i) - (v_y*v_y)*_mass_accum.at(i)) / numMols_accum;
                        T_z         = (2*_ekinVect_accum[2].at(i) - (v_z*v_z)*_mass_accum.at(i)) / numMols_accum;
                        p_x         = rho * ( _virialVect_accum[0].at(i)/numMols_accum + T);
                        p_y         = rho * ( _virialVect_accum[1].at(i)/numMols_accum + T);
                        p_z         = rho * ( _virialVect_accum[2].at(i)/numMols_accum + T);
                        F_x         = _forceVect_accum[0].at(i)      / numMols_accum;
                        F_y         = _forceVect_accum[1].at(i)      / numMols_accum;
                        F_z         = _forceVect_accum[2].at(i)      / numMols_accum;
                        jEF_x       = _energyfluxVect_accum[0].at(i) / countSamples;
                        jEF_y       = _energyfluxVect_accum[1].at(i) / countSamples;
                        jEF_z       = _energyfluxVect_accum[2].at(i) / countSamples;

                        numSamples  = countSamples;
                    }
                    if ((_chemPot_accum.at(i) > 0.0) and (_countNTest_accum.at(i) > 0ul)) {
                        numTest     = static_cast<double>(_countNTest_accum.at(i)*_samplefrequency) / (2*_writeFrequency);
                        chemPot_res = -log(_chemPot_accum.at(i)/_countNTest_accum.at(i)) + log(rho);  // Implemented in accordance to ms2
                    }
                    ofs << FORMAT_SCI_MAX_DIGITS << numMolsPerStep
                        << FORMAT_SCI_MAX_DIGITS << rho
                        << FORMAT_SCI_MAX_DIGITS << T
                        << FORMAT_SCI_MAX_DIGITS << ekin
                        << FORMAT_SCI_MAX_DIGITS << epot
                        << FORMAT_SCI_MAX_DIGITS << orderparam
                        << FORMAT_SCI_MAX_DIGITS << p
                        << FORMAT_SCI_MAX_DIGITS << chemPot_res
                        << FORMAT_SCI_MAX_DIGITS << numTest
                        << FORMAT_SCI_MAX_DIGITS << T_x
                        << FORMAT_SCI_MAX_DIGITS << T_y
                        << FORMAT_SCI_MAX_DIGITS << T_z
                        << FORMAT_SCI_MAX_DIGITS << v_x
                        << FORMAT_SCI_MAX_DIGITS << v_y
                        << FORMAT_SCI_MAX_DIGITS << v_z
                        << FORMAT_SCI_MAX_DIGITS << p_x
                        << FORMAT_SCI_MAX_DIGITS << p_y
                        << FORMAT_SCI_MAX_DIGITS << p_z
                        << FORMAT_SCI_MAX_DIGITS << F_x
                        << FORMAT_SCI_MAX_DIGITS << F_y
                        << FORMAT_SCI_MAX_DIGITS << F_z
                        << FORMAT_SCI_MAX_DIGITS << jEF_x
                        << FORMAT_SCI_MAX_DIGITS << jEF_y
                        << FORMAT_SCI_MAX_DIGITS << jEF_z
                        << FORMAT_SCI_MAX_DIGITS << numSamples;
                }
                ofs << std::endl;
            }
            ofs.close();
            }

            // Write output for higher moments
            if (_sampleHigherMoms) {
                std::map<std::string, unsigned short> dirs { {"x", 0}, {"y", 1}, {"z", 2}, };
                std::stringstream ss;
                ss << std::setw(9) << std::setfill('0') << simstep;
                const std::string fname = "ExtendedProfileSampling_HigherMoments_TS"+ss.str()+".dat";
                std::ofstream ofs;
                ofs.open(fname, std::ios::out);
                ofs << setw(24) << "pos";                           // Bin position
                for (unsigned long cid = 0; cid < numOutputs; cid++) {
                    ofs << setw(22) << "delta[" << cid << "]";
                    // key k_i and value v_i
                    for (auto& [k_i, v_i] : dirs) {
                        ofs << setw(20) << "q_" << k_i << "[" << cid << "]";
                    }
                    for (auto& [k_i, v_i] : dirs) {
                        for (auto& [k_j, v_j] : dirs) {
                            ofs << setw(19) << "p_" << k_i << k_j << "[" << cid << "]";
                        }
                    }
                    for (auto& [k_i, v_i] : dirs) {
                        for (auto& [k_j, v_j] : dirs) {
                            ofs << setw(19) << "R_" << k_i << k_j << "[" << cid << "]";
                        }
                    }
                    for (auto& [k_i, v_i] : dirs) {
                        for (auto& [k_j, v_j] : dirs) {
                            for (auto& [k_k, v_k] : dirs) {
                                ofs << setw(18) << "m_" << k_i << k_j << k_k << "[" << cid << "]";
                            }
                        }
                    }
                }
                ofs << std::endl;

                for (unsigned long idx = 0; idx < _numBinsGlobal; idx++) {
                    ofs << FORMAT_SCI_MAX_DIGITS << (idx+0.5)*_binwidth;
                    for (unsigned long cid = 0; cid < numOutputs; cid++) {
                        unsigned long i = idx + cid*_numBinsGlobal;
                        double delta {0.0};
                        std::array<double, 3> q = {0.0};
                        std::array<double, 9> p = {0.0};
                        std::array<double, 9> R = {0.0};
                        std::array<double, 27> m = {0.0};
                        if (_countSamples.at(i) > 0ul) {
                            delta = _hmDelta_accum.at(i)/_countSamples.at(i);
                            for (unsigned short d = 0; d < 3; d++) {
                                q[d] = _hmHeatflux_accum[d].at(i)/_countSamples.at(i);
                            }
                            for (unsigned short d = 0; d < 9; d++) {
                                p[d] = _hmPressure_accum[d].at(i)/_countSamples.at(i);
                                R[d] = _hmR_accum[d].at(i)/_countSamples.at(i);
                                m[d]    = _hmM_accum[d].at(i)/_countSamples.at(i);
                                m[d+9u]  = _hmM_accum[d+9u].at(i)/_countSamples.at(i);
                                m[d+18u] = _hmM_accum[d+18u].at(i)/_countSamples.at(i);
                            }
                        }
                        ofs << FORMAT_SCI_MAX_DIGITS << delta;
                        // key k_i and value v_i
                        for (auto& [k_i, v_i] : dirs) {
                            ofs << FORMAT_SCI_MAX_DIGITS << q[v_i];
                        }
                        for (auto& [k_i, v_i] : dirs) {
                            for (auto& [k_j, v_j] : dirs) {
                                ofs << FORMAT_SCI_MAX_DIGITS << p[3u*v_i+v_j];
                            }
                        }
                        for (auto& [k_i, v_i] : dirs) {
                            for (auto& [k_j, v_j] : dirs) {
                                ofs << FORMAT_SCI_MAX_DIGITS << R[3u*v_i+v_j];
                            }
                        }
                        for (auto& [k_i, v_i] : dirs) {
                            for (auto& [k_j, v_j] : dirs) {
                                for (auto& [k_k, v_k] : dirs) {
                                    ofs << FORMAT_SCI_MAX_DIGITS << m[9u*v_i+3u*v_j+v_k];
                                }
                            }
                        }
                    }
                    ofs << std::endl;
                }
                ofs.close();
            }
        }

        // Reset vectors to zero
        resetVectors();
    }
}

// Resize vectors
void ExtendedProfileSampling::resizeVectors() {

    _numMolecules_accum.resize(_lenVector);
    _doftotal_accum.resize(_lenVector);
    _mass_accum.resize(_lenVector);
    _ekin_accum.resize(_lenVector);
    _epot_accum.resize(_lenVector);
    _orientation_accum.resize(_lenVector);
    _virial_accum.resize(_lenVector);
    _chemPot_accum.resize(_lenVector);
    _countNTest_accum.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        _ekinVect_accum.at(d).resize(_lenVector);
        _velocityVect_accum.at(d).resize(_lenVector);
        _virialVect_accum.at(d).resize(_lenVector);
        _forceVect_accum.at(d).resize(_lenVector);
        _energyfluxVect_accum.at(d).resize(_lenVector);
    }

    _countSamples.resize(_lenVector);

    if (_sampleHigherMoms) {
        _hmDelta_accum.resize(_lenVector);
        for (unsigned short d = 0; d < 3; d++) {
            _hmHeatflux_accum.at(d).resize(_lenVector);
        }
        for (unsigned short d = 0; d < 9; d++) {
            _hmPressure_accum.at(d).resize(_lenVector);
            _hmR_accum.at(d).resize(_lenVector);
            _hmM_accum.at(d).resize(_lenVector);
            _hmM_accum.at(d+9u).resize(_lenVector);
            _hmM_accum.at(d+18u).resize(_lenVector);
        }
    } else {
        _hmDelta_accum.clear();
        for (unsigned short d = 0; d < 3; d++) {
            _hmHeatflux_accum.at(d).clear();
        }
        for (unsigned short d = 0; d < 9; d++) {
            _hmPressure_accum.at(d).clear();
            _hmR_accum.at(d).clear();
            _hmM_accum.at(d).clear();
            _hmM_accum.at(d+9u).clear();
            _hmM_accum.at(d+18u).clear();
        }
    }

}

// Fill vectors with zeros
void ExtendedProfileSampling::resetVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    std::fill(_doftotal_accum.begin(), _doftotal_accum.end(), 0ul);
    std::fill(_mass_accum.begin(), _mass_accum.end(), 0.0f);
    std::fill(_ekin_accum.begin(), _ekin_accum.end(), 0.0f);
    std::fill(_epot_accum.begin(), _epot_accum.end(), 0.0f);
    std::fill(_orientation_accum.begin(), _orientation_accum.end(), 0.0f);
    std::fill(_virial_accum.begin(), _virial_accum.end(), 0.0f);
    std::fill(_chemPot_accum.begin(), _chemPot_accum.end(), 0.0f);
    std::fill(_countNTest_accum.begin(), _countNTest_accum.end(), 0ul);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(_ekinVect_accum.at(d).begin(), _ekinVect_accum.at(d).end(), 0.0f);
        std::fill(_velocityVect_accum.at(d).begin(), _velocityVect_accum.at(d).end(), 0.0f);
        std::fill(_virialVect_accum.at(d).begin(), _virialVect_accum.at(d).end(), 0.0f);
        std::fill(_forceVect_accum.at(d).begin(), _forceVect_accum.at(d).end(), 0.0f);
        std::fill(_energyfluxVect_accum.at(d).begin(), _energyfluxVect_accum.at(d).end(), 0.0f);
    }

    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);

    if (_sampleHigherMoms) {
        std::fill(_hmDelta_accum.begin(), _hmDelta_accum.end(), 0.0f);
        for (unsigned short d = 0; d < 3; d++) {
            std::fill(_hmHeatflux_accum.at(d).begin(), _hmHeatflux_accum.at(d).end(), 0.0f);
        }
        for (unsigned short d = 0; d < 9; d++) {
            std::fill(_hmPressure_accum.at(d).begin(), _hmPressure_accum.at(d).end(), 0.0f);
            std::fill(_hmR_accum.at(d).begin(), _hmR_accum.at(d).end(), 0.0f);
            std::fill(_hmM_accum.at(d).begin(), _hmM_accum.at(d).end(), 0.0f);
            std::fill(_hmM_accum.at(d+9u).begin(), _hmM_accum.at(d+9u).end(), 0.0f);
            std::fill(_hmM_accum.at(d+18u).begin(), _hmM_accum.at(d+18u).end(), 0.0f);
        }
    }
}

// Get value of quantity at certain index; Mainly used for unit test
double ExtendedProfileSampling::getQuantity(DomainDecompBase* domainDecomp, std::string quantityName, unsigned long index) {
    if (index > _lenVector) {
        global_log->error() << "[ExtendedProfileSampling] Trying to get value but index (" << index << ") is greater than possible (" << _lenVector << ")" << std::endl;
        return 0.0;
    }
    // Only root knows real quantities (MPI_Reduce instead of MPI_Allreduce at accumulation)
    if (domainDecomp->getRank() != 0) {
        global_log->error() << "[ExtendedProfileSampling] Non-root process tried to get value" << std::endl;
        return 0.0;
    }

    if (_countSamples.at(index) == 0) {
        global_log->error() << "[ExtendedProfileSampling] Nothing sampled yet at given index" << std::endl;
        return 0.0;
    }

    if (quantityName == "T") {
        const double v_x = _velocityVect_accum[0].at(index) / _countSamples.at(index);
        const double v_y = _velocityVect_accum[1].at(index) / _countSamples.at(index);
        const double v_z = _velocityVect_accum[2].at(index) / _countSamples.at(index);
        double v_drift_sqr = v_x*v_x + v_y*v_y + v_z*v_z;
        return (2*_ekin_accum.at(index) - v_drift_sqr*_mass_accum.at(index)) / _doftotal_accum.at(index);
    } else if (quantityName == "rho") {
        return _numMolecules_accum.at(index) / (_slabVolume * _countSamples.at(index));
    } else if (quantityName == "ekin") {
        return _ekin_accum.at(index) / _countSamples.at(index);
    }
    
    global_log->error() << "[ExtendedProfileSampling] Quantity (" << quantityName << ") unknown!" << std::endl;
    return 0.0;
}
