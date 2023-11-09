/*
 * CylindricSampling.cpp
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#include "CylindricSampling.h"

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


CylindricSampling::CylindricSampling() {}

void CylindricSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    _distMax = 0.5*min(_globalBoxLength[0],_globalBoxLength[2]);  // Only sample particles within maximum radius fitting into box

    _numBinsGlobalHeight = static_cast<unsigned int>(_globalBoxLength[1]/_binwidth);
    _numBinsGlobalRadius = static_cast<unsigned int>(_distMax/_binwidth);

    if (_globalBoxLength[1]/_binwidth != static_cast<float>(_numBinsGlobalHeight)) {
        global_log->error() << "[CylindricSampling] Can not divide domain without remainder! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }
    if (_distMax/_binwidth != static_cast<float>(_numBinsGlobalRadius)) {
        global_log->error() << "[CylindricSampling] Can not divide domain without remainder! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }

     // Entry per bin; all components sampled as one
    _lenVector = _numBinsGlobalRadius * _numBinsGlobalHeight;

    resizeVectors();
    resetVectors();
}

void CylindricSampling::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("binwidth", _binwidth);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("stop", _stopSampling);


    global_log->info() << "[CylindricSampling] Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    global_log->info() << "[CylindricSampling] Binwidth: " << _binwidth << std::endl;
    global_log->info() << "[CylindricSampling] All components treated as single one" << std::endl;
}

void CylindricSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

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
    std::array<CommVar<std::vector<double>>, 3> ekinVect_step;
    std::array<CommVar<std::vector<double>>, 3> velocityVect_step;
    std::array<CommVar<std::vector<double>>, 3> virialVect_step;
    std::array<CommVar<std::vector<double>>, 3> forceVect_step;

    std::array<std::vector<double>, 3> veloDrift_step_global;       // Drift velocity

    numMolecules_step.local.resize(_lenVector);
    mass_step.local.resize(_lenVector);
    ekin_step.local.resize(_lenVector);

    numMolecules_step.global.resize(_lenVector);
    mass_step.global.resize(_lenVector);
    ekin_step.global.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        ekinVect_step[d].local.resize(_lenVector);
        velocityVect_step[d].local.resize(_lenVector);
        virialVect_step[d].local.resize(_lenVector);
        forceVect_step[d].local.resize(_lenVector);

        ekinVect_step[d].global.resize(_lenVector);
        velocityVect_step[d].global.resize(_lenVector);
        virialVect_step[d].global.resize(_lenVector);
        forceVect_step[d].global.resize(_lenVector);

        veloDrift_step_global[d].resize(_lenVector);
    }

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);
    std::fill(mass_step.local.begin(),        mass_step.local.end(), 0.0f);
    std::fill(ekin_step.local.begin(),        ekin_step.local.end(), 0.0f);

    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);
    std::fill(mass_step.global.begin(),        mass_step.global.end(), 0.0f);
    std::fill(ekin_step.global.begin(),        ekin_step.global.end(), 0.0f);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(ekinVect_step[d].local.begin(),      ekinVect_step[d].local.end(), 0.0f);
        std::fill(velocityVect_step[d].local.begin(),   velocityVect_step[d].local.end(), 0.0f);
        std::fill(virialVect_step[d].local.begin(),     virialVect_step[d].local.end(), 0.0f);
        std::fill(forceVect_step[d].local.begin(),      forceVect_step[d].local.end(), 0.0f);

        std::fill(ekinVect_step[d].global.begin(),      ekinVect_step[d].global.end(), 0.0f);
        std::fill(velocityVect_step[d].global.begin(),   velocityVect_step[d].global.end(), 0.0f);
        std::fill(virialVect_step[d].global.begin(),     virialVect_step[d].global.end(), 0.0f);
        std::fill(forceVect_step[d].global.begin(),      forceVect_step[d].global.end(), 0.0f);

        std::fill(veloDrift_step_global[d].begin(),      veloDrift_step_global[d].end(), 0.0f);
    }

    // Calculate drift as it is needed first
    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double ry = pit->r(1);
        const double distCenter = std::sqrt(std::pow(pit->r(0)-0.5*_globalBoxLength[0],2) + std::pow(pit->r(2)-0.5*_globalBoxLength[0],2));
        // Do not consider particles outside of most outer radius
        if (distCenter >= _distMax) { continue; }
        const unsigned int indexH = std::min(_numBinsGlobalHeight, static_cast<unsigned int>(ry/_binwidth));  // Index of bin of height
        const unsigned int indexR = std::min(_numBinsGlobalRadius, static_cast<unsigned int>(distCenter/_binwidth));  // Index of bin of radius
        const unsigned int index = _numBinsGlobalHeight*indexR + indexH;

        numMolecules_step.local[index]) ++;
        velocityVect_step[0].local[index]) += pit->v(0);
        velocityVect_step[1].local[index]) += pit->v(1);
        velocityVect_step[2].local[index]) += pit->v(2);
    }

#ifdef ENABLE_MPI
    MPI_Allreduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[0].local.data(), velocityVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[1].local.data(), velocityVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[2].local.data(), velocityVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global[i] = numMolecules_step.local[i];
        velocityVect_step[0].global[i] = velocityVect_step[0].local[i];
        velocityVect_step[1].global[i] = velocityVect_step[1].local[i];
        velocityVect_step[2].global[i] = velocityVect_step[2].local[i];
    }
#endif

    for (unsigned long i = 0; i < _lenVector; i++) {
        if (numMolecules_step.global[i] > 0ul) {
            veloDrift_step_global[0][i] = velocityVect_step[0].global[i] / numMolecules_step.global[i];
            veloDrift_step_global[1][i] = velocityVect_step[1].global[i] / numMolecules_step.global[i];
            veloDrift_step_global[2][i] = velocityVect_step[2].global[i] / numMolecules_step.global[i];
        }
    }

    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double ry = pit->r(1);
        const double distCenter = std::sqrt(std::pow(pit->r(0)-0.5*_globalBoxLength[0],2) + std::pow(pit->r(2)-0.5*_globalBoxLength[0],2));
        // Do not consider particles outside of most outer radius
        if (distCenter >= _distMax) { continue; }
        const unsigned int indexH = std::min(_numBinsGlobalHeight, static_cast<unsigned int>(ry/_binwidth));  // Index of bin of height
        const unsigned int indexR = std::min(_numBinsGlobalRadius, static_cast<unsigned int>(distCenter/_binwidth));  // Index of bin of radius
        const unsigned int index = _numBinsGlobalHeight*indexR + indexH;

        const double veloX = pit->v(0);
        const double veloY = pit->v(1);
        const double veloZ = pit->v(2);
        const double mass = pit->mass();

        mass_step.local[index]) += pit->mass();
        ekin_step.local[index]) += pit->U_kin();

        ekinVect_step[0].local[index]) += 0.5*mass*veloX*veloX;
        ekinVect_step[1].local[index]) += 0.5*mass*veloY*veloY;
        ekinVect_step[2].local[index]) += 0.5*mass*veloZ*veloZ;
        virialVect_step[0].local[index]) += pit->Vi(0);
        virialVect_step[1].local[index]) += pit->Vi(1);
        virialVect_step[2].local[index]) += pit->Vi(2);
        forceVect_step[0].local[index]) += pit->F(0);
        forceVect_step[1].local[index]) += pit->F(1);
        forceVect_step[2].local[index]) += pit->F(2);

    }

    // Gather quantities of all processes. Note: MPI_Reduce instead of MPI_Allreduce!
#ifdef ENABLE_MPI
    MPI_Reduce(mass_step.local.data(), mass_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin_step.local.data(), ekin_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[0].local.data(), virialVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[1].local.data(), virialVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[2].local.data(), virialVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(forceVect_step[0].local.data(), forceVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(forceVect_step[1].local.data(), forceVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(forceVect_step[2].local.data(), forceVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[0].local.data(), ekinVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[1].local.data(), ekinVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[2].local.data(), ekinVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        mass_step.global[i] = mass_step.local[i];
        ekin_step.global[i] = ekin_step.local[i];
        virialVect_step[0].global[i] = virialVect_step[0].local[i];
        virialVect_step[1].global[i] = virialVect_step[1].local[i];
        virialVect_step[2].global[i] = virialVect_step[2].local[i];
        forceVect_step[0].global[i] = forceVect_step[0].local[i];
        forceVect_step[1].global[i] = forceVect_step[1].local[i];
        forceVect_step[2].global[i] = forceVect_step[2].local[i];
        ekinVect_step[0].global[i] = ekinVect_step[0].local[i];
        ekinVect_step[1].global[i] = ekinVect_step[1].local[i];
        ekinVect_step[2].global[i] = ekinVect_step[2].local[i];
    }
#endif

    // Only root knows real quantities (MPI_Reduce instead of MPI_Allreduce)
    // Accumulate data
    if (domainDecomp->getRank() == 0) {
        for (unsigned long i = 0; i < _lenVector; i++) {
            const unsigned long numMols = numMolecules_step.global[i];
            unsigned int dof_rot {0};
            unsigned int dof_total {0};

            // For single component sampling, the rot. DOF of component 0 is taken
            dof_rot = _simulation.getEnsemble()->getComponent(0)->getRotationalDegreesOfFreedom();
            dof_total = (3 + dof_rot)*numMols;

            const double ViX = virialVect_step[0].global[i];
            const double ViY = virialVect_step[1].global[i];
            const double ViZ = virialVect_step[2].global[i];

            _doftotal_accum[i]                += dof_total;
            _numMolecules_accum[i]            += numMols;
            _mass_accum[i]                    += mass_step.global[i];
            _ekin_accum[i]                    += ekin_step.global[i];
            _virial_accum[i]                  += ViX + ViY + ViZ;

            _ekinVect_accum[0][i]             += ekinVect_step[0].global[i];
            _ekinVect_accum[1][i]             += ekinVect_step[1].global[i];
            _ekinVect_accum[2][i]             += ekinVect_step[2].global[i];
            _velocityVect_accum[0][i]         += veloDrift_step_global[0][i];
            _velocityVect_accum[1][i]         += veloDrift_step_global[1][i];
            _velocityVect_accum[2][i]         += veloDrift_step_global[2][i];
            _virialVect_accum[0][i]           += ViX;
            _virialVect_accum[1][i]           += ViY;
            _virialVect_accum[2][i]           += ViZ;
            _forceVect_accum[0][i]            += forceVect_step[0].global[i];
            _forceVect_accum[1][i]            += forceVect_step[1].global[i];
            _forceVect_accum[2][i]            += forceVect_step[2].global[i];

            if (numMols > 0ul) { _countSamples[i]++; }
        }
    }

    // Write out data every _writeFrequency step
    if ( (simstep - _startSampling) % _writeFrequency == 0 ) {

        if (domainDecomp->getRank() == 0) {
            // Write output file
            std::stringstream ss;
            ss << std::setw(9) << std::setfill('0') << simstep;
            const std::string fname = "CylindricSampling_TS"+ss.str()+".dat";
            std::ofstream ofs;
            ofs.open(fname, std::ios::out);
            ofs << setw(24) << "height"             // Bin position (height)
                << setw(24) << "radius"        // Bin position (radius)
                << setw(24) << "numParts"        // Average number of molecules in bin per step
                << setw(24) << "rho"        // Density
                << setw(24) << "T"        // Temperature without drift (i.e. "real" temperature)
                << setw(24) << "ekin"        // Kinetic energy including drift
                << setw(24) << "p"        // Pressure
                << setw(24) << "T_x"        // Temperature in x-direction
                << setw(24) << "T_y"        // Temperature in y-direction
                << setw(24) << "T_z"        // Temperature in z-direction
                << setw(24) << "v_x"        // Drift velocity in x-direction
                << setw(24) << "v_y"        // Drift velocity in y-direction
                << setw(24) << "v_z"        // Drift velocity in z-direction
                << setw(24) << "p_x"        // Pressure in x-direction
                << setw(24) << "p_y"        // Pressure in y-direction
                << setw(24) << "p_z"        // Pressure in z-direction
                << setw(24) << "F_x"        // Force in x-direction
                << setw(24) << "F_y"        // Force in y-direction
                << setw(24) << "F_z"        // Force in z-direction
                << setw(24) << "numSamples";    // Number of samples (<= _writeFrequency)
            ofs << std::endl;

            for (unsigned long i = 0; i < _lenVector; i++) {
                const unsigned int idxH = i % _numBinsGlobalHeight;  // Remainder
                const unsigned int idxR = i / _numBinsGlobalHeight; // Quotient
                ofs << FORMAT_SCI_MAX_DIGITS << (idxH+0.5)*_binwidth;  // Height bin
                ofs << FORMAT_SCI_MAX_DIGITS << (idxR+0.5)*_binwidth;  // Radius bin
                double numMolsPerStep {std::nan("0")}; // Not an int as particles change bin during simulation
                double rho {std::nan("0")};
                double T {std::nan("0")};
                double ekin {std::nan("0")};
                double p {std::nan("0")};
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
                double numSamples {std::nan("0")};
                if ((_countSamples[i] > 0ul) and (_doftotal_accum[i]) > 0ul) {
                    const unsigned long countSamples = _countSamples[i];
                    const double numMols_accum = static_cast<double>(_numMolecules_accum[i]);
                    const double slabVolume = 3.1415926536*_binwidth*(std::pow((idxR+1)*_binwidth,2)-std::pow(idxR*_binwidth,2));

                    numMolsPerStep = numMols_accum/_writeFrequency;
                    rho         = numMolsPerStep           / slabVolume;
                    v_x         = _velocityVect_accum[0][i]   / countSamples;
                    v_y         = _velocityVect_accum[1][i]   / countSamples;
                    v_z         = _velocityVect_accum[2][i]   / countSamples;

                    double v_drift_sqr = v_x*v_x + v_y*v_y + v_z*v_z;

                    T           = (2*_ekin_accum[i] - v_drift_sqr*_mass_accum[i]) / _doftotal_accum[i];
                    ekin        = _ekin_accum[i] / numMols_accum;
                    p           = rho * ( (_virial_accum[i])/(3.0*numMols_accum) + T);

                    T_x         = (2*_ekinVect_accum[0][i] - (v_x*v_x)*_mass_accum[i]) / numMols_accum;
                    T_y         = (2*_ekinVect_accum[1][i] - (v_y*v_y)*_mass_accum[i]) / numMols_accum;
                    T_z         = (2*_ekinVect_accum[2][i] - (v_z*v_z)*_mass_accum[i]) / numMols_accum;
                    p_x         = rho * ( _virialVect_accum[0][i]/numMols_accum + T);
                    p_y         = rho * ( _virialVect_accum[1][i]/numMols_accum + T);
                    p_z         = rho * ( _virialVect_accum[2][i]/numMols_accum + T);
                    F_x         = _forceVect_accum[0][i]      / numMols_accum;
                    F_y         = _forceVect_accum[1][i]      / numMols_accum;
                    F_z         = _forceVect_accum[2][i]      / numMols_accum;

                    numSamples  = countSamples;
                }
                ofs << FORMAT_SCI_MAX_DIGITS << numMolsPerStep
                    << FORMAT_SCI_MAX_DIGITS << rho
                    << FORMAT_SCI_MAX_DIGITS << T
                    << FORMAT_SCI_MAX_DIGITS << ekin
                    << FORMAT_SCI_MAX_DIGITS << p
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
                    << FORMAT_SCI_MAX_DIGITS << numSamples
                    << std::endl;
            }
            ofs.close();
        }

        // Reset vectors to zero
        resetVectors();
    }
}

// Resize vectors
void CylindricSampling::resizeVectors() {

    _numMolecules_accum.resize(_lenVector);
    _doftotal_accum.resize(_lenVector);
    _mass_accum.resize(_lenVector);
    _ekin_accum.resize(_lenVector);
    _virial_accum.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        _ekinVect_accum[d].resize(_lenVector);
        _velocityVect_accum[d].resize(_lenVector);
        _virialVect_accum[d].resize(_lenVector);
        _forceVect_accum[d].resize(_lenVector);
    }

    _countSamples.resize(_lenVector);
}

// Fill vectors with zeros
void CylindricSampling::resetVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    std::fill(_doftotal_accum.begin(), _doftotal_accum.end(), 0ul);
    std::fill(_mass_accum.begin(), _mass_accum.end(), 0.0f);
    std::fill(_ekin_accum.begin(), _ekin_accum.end(), 0.0f);
    std::fill(_virial_accum.begin(), _virial_accum.end(), 0.0f);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(_ekinVect_accum[d].begin(), _ekinVect_accum[d].end(), 0.0f);
        std::fill(_velocityVect_accum[d].begin(), _velocityVect_accum[d].end(), 0.0f);
        std::fill(_virialVect_accum[d].begin(), _virialVect_accum[d].end(), 0.0f);
        std::fill(_forceVect_accum[d].begin(), _forceVect_accum[d].end(), 0.0f);
    }

    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);
}

