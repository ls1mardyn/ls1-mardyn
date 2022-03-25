/*
 * ExtendedProfileSampling.cpp
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#include "ExtendedProfileSampling.h"

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "utils/Random.h"
#include "Simulation.h"

#include <math.h>


ExtendedProfileSampling::ExtendedProfileSampling()
    :
    // Control: general
    _binwidth(1.0f),
    _startSampling(0ul),
    _writeFrequency(10000ul),
    _stopSampling(1000000000ul),
    _singleComp(false),

    // Control: chemical potential
    _sampleChemPot(false),
    _lattice(true),
    _factorNumTest(4.0f),
    _samplefrequency(50ul)
{}

ExtendedProfileSampling::~ExtendedProfileSampling() {
    delete _particlePairsHandler; // ???
}

void ExtendedProfileSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    _numBinsGlobal = static_cast<uint16_t>(_globalBoxLength[1]/_binwidth);
    if (_globalBoxLength[1]/_binwidth != static_cast<float>(_numBinsGlobal)) {
        global_log->error() << "[ExtendedProfileSampling] Can not divide domain without remainder! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }
    _slabVolume = _globalBoxLength[0]*_globalBoxLength[2]*_binwidth;

    if (_slabVolume < 1e-12) {
        global_log->error() << "[ExtendedProfileSampling] Slab volume too small!" << std::endl;
        Simulation::exit(-1);
    }

     // Entry per component and bin; 0 represents all components combined
    _lenVector = (_singleComp) ? _numBinsGlobal : _numBinsGlobal * (domain->getNumberOfComponents()+1);

    resizeVectors();
    resetVectors();

    // ??? TODO: CAN CELL PROCESSOR CHANGE DURING SIMULATION???
    _cellProcessor = _simulation.getCellProcessor();
    // ??? TODO: CAN PP HANDLER CHANGE DURING SIMULATION???
    _particlePairsHandler = new ParticlePairs2PotForceAdapter(*domain);
    // MolID is maximum possible number minus rank to prevent duplicate IDs
    // Always insert molecule of first component
    const unsigned long molID = std::numeric_limits<unsigned long>::max() - static_cast<unsigned long>(domainDecomp->getRank());
    _mTest = Molecule(molID, &(_simulation.getEnsemble()->getComponents()->at(0)));
}

void ExtendedProfileSampling::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("binwidth", _binwidth);  // Default: 1.0
    xmlconfig.getNodeValue("start", _startSampling);  // Default: 0
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);  // Default: 10000
    xmlconfig.getNodeValue("stop", _stopSampling);  // Default: 1000000000
    xmlconfig.getNodeValue("singlecomponent", _singleComp);  // Default: false

    xmlconfig.getNodeValue("@enable", _sampleChemPot);  // Default: false
    bool bCP =   xmlconfig.getNodeValue("chemicalpotential/lattice", _lattice);  // Default: true
    bCP = bCP || xmlconfig.getNodeValue("chemicalpotential/factorNumTest", _factorNumTest);  // Default: 4.0
    bCP = bCP || xmlconfig.getNodeValue("chemicalpotential/samplefrequency", _samplefrequency);  // Default: 50

    if (bCP and !_sampleChemPot) {
        global_log->warning() << "[ExtendedProfileSampling] Sampling of chem. pot. disabled, but control values in config specified! " << std::endl;
    }

    string insMethod;
    if (_lattice) {
        insMethod = "in a lattice";
    } else {
        insMethod = "randomly";
    }

    global_log->info() << "[ExtendedProfileSampling] Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    global_log->info() << "[ExtendedProfileSampling] Binwidth: " << _binwidth << std::endl;
    if (_sampleChemPot) {
        global_log->info() << "[ExtendedProfileSampling] Sampling of chemical potential enabled with a sampling frequency of " << _samplefrequency << std::endl;
        global_log->info() << "[ExtendedProfileSampling] " << _factorNumTest << " * numParticles will be inserted " << insMethod << std::endl;
    }

    if (_samplefrequency > _writeFrequency) {
        global_log->warning() << "[ExtendedProfileSampling] Sample frequency is greater than write frequency! " << std::endl;
    }
}

// Needs to be called when halo cells are still existing (for sampling of chemical potential)
void ExtendedProfileSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    // Sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep <= _startSampling) or (simstep > _stopSampling)) {
        return;
    }

    // Sample only every _samplefrequency
	if ( (simstep - _startSampling) % _samplefrequency != 0 ) {
        return;
    }

    // Do not write or sample data directly after (re)start in the first step
	if ( simstep == global_simulation->getNumInitTimesteps() ) {
        return;
    }
    
    // Variables per step
    std::array<double, 3> regionLowCorner;
    std::array<double, 3> regionHighCorner;
    std::array<double, 3> regionSize;
    for (unsigned short d = 0; d < 3; d++) {
        regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
        regionSize[d] = regionHighCorner[d] - regionLowCorner[d];
    }

    CommVar<std::vector<unsigned long>> numMolecules_step;
    CommVar<std::vector<double>> ekin2_step;                        // Without drift energy
    CommVar<std::vector<double>> ekin2Trans_step;                   // Including drift energy
    CommVar<std::vector<double>> epot_step;
    CommVar<std::vector<double>> chemPot_step;
    CommVar<std::vector<unsigned long>> countNTest_step;
    std::array<CommVar<std::vector<double>>, 3> ekin2Vect_step;
    std::array<CommVar<std::vector<double>>, 3> velocityVect_step;
    std::array<CommVar<std::vector<double>>, 3> virialVect_step;
    std::array<CommVar<std::vector<double>>, 3> energyfluxVect_step;

    std::array<std::vector<double>, 3> veloDrift_step_global;       // Drift velocity per particle; global value as calculated with global values
    std::vector<double> temperature_step_global (_lenVector, 0.0);  // Required for sampling of chem. pot.

    numMolecules_step.local.resize(_lenVector);
    ekin2_step.local.resize(_lenVector);
    ekin2Trans_step.local.resize(_lenVector);
    epot_step.local.resize(_lenVector);
    chemPot_step.local.resize(_lenVector);
    countNTest_step.local.resize(_lenVector);
    
    for (unsigned short d = 0; d < 3; d++) {
        ekin2Vect_step.at(d).local.resize(_lenVector);
        velocityVect_step.at(d).local.resize(_lenVector);
        virialVect_step.at(d).local.resize(_lenVector);
        energyfluxVect_step.at(d).local.resize(_lenVector);

        veloDrift_step_global.at(d).resize(_lenVector);
    }

    numMolecules_step.global.resize(_lenVector);
    ekin2_step.global.resize(_lenVector);
    ekin2Trans_step.global.resize(_lenVector);
    epot_step.global.resize(_lenVector);
    chemPot_step.global.resize(_lenVector);
    countNTest_step.global.resize(_lenVector);
    
    for (unsigned short d = 0; d < 3; d++) {
        ekin2Vect_step.at(d).global.resize(_lenVector);
        velocityVect_step.at(d).global.resize(_lenVector);
        virialVect_step.at(d).global.resize(_lenVector);
        energyfluxVect_step.at(d).global.resize(_lenVector);
    }

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);
    std::fill(ekin2_step.local.begin(), epot_step.local.end(), 0.0f);
    std::fill(ekin2Trans_step.local.begin(), epot_step.local.end(), 0.0f);
    std::fill(epot_step.local.begin(), epot_step.local.end(), 0.0f);
    std::fill(chemPot_step.local.begin(), chemPot_step.local.end(), 0.0f);
    std::fill(countNTest_step.local.begin(), countNTest_step.local.end(), 0ul);
    
    for (unsigned short d = 0; d < 3; d++) {
        std::fill(ekin2Vect_step.at(d).local.begin(), ekin2Vect_step.at(d).local.end(), 0.0f);
        std::fill(velocityVect_step.at(d).local.begin(), velocityVect_step.at(d).local.end(), 0.0f);
        std::fill(virialVect_step.at(d).local.begin(), virialVect_step.at(d).local.end(), 0.0f);
        std::fill(energyfluxVect_step.at(d).local.begin(), energyfluxVect_step.at(d).local.end(), 0.0f);

        std::fill(veloDrift_step_global.at(d).begin(), veloDrift_step_global.at(d).end(), 0.0f);
    }

    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);
    std::fill(ekin2_step.global.begin(), epot_step.global.end(), 0.0f);
    std::fill(ekin2Trans_step.global.begin(), epot_step.global.end(), 0.0f);
    std::fill(epot_step.global.begin(), epot_step.global.end(), 0.0f);
    std::fill(chemPot_step.global.begin(), chemPot_step.global.end(), 0.0f);
    std::fill(countNTest_step.global.begin(), countNTest_step.global.end(), 0ul);
    
    for (unsigned short d = 0; d < 3; d++) {
        std::fill(ekin2Vect_step.at(d).global.begin(), ekin2Vect_step.at(d).global.end(), 0.0f);
        std::fill(velocityVect_step.at(d).global.begin(), velocityVect_step.at(d).global.end(), 0.0f);
        std::fill(virialVect_step.at(d).global.begin(), virialVect_step.at(d).global.end(), 0.0f);
        std::fill(energyfluxVect_step.at(d).global.begin(), energyfluxVect_step.at(d).global.end(), 0.0f);
    }

    // Calculate drift as it is needed first
    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double ry = pit->r(1);
        const uint16_t index = std::min(_numBinsGlobal, static_cast<uint16_t>(ry/_binwidth));  // Index of bin

        numMolecules_step.local.at(index) ++;
        velocityVect_step.at(0).local.at(index) += pit->v(0);
        velocityVect_step.at(1).local.at(index) += pit->v(1);
        velocityVect_step.at(2).local.at(index) += pit->v(2);

        if (!_singleComp) {
            const unsigned int cid = pit->componentid() + 1; // 0 represents all components
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
        const uint16_t index = std::min(_numBinsGlobal, static_cast<uint16_t>(ry/_binwidth));  // Index of bin

        const double veloCorrX = pit->v(0) - veloDrift_step_global[0].at(index);
        const double veloCorrY = pit->v(1) - veloDrift_step_global[1].at(index);
        const double veloCorrZ = pit->v(2) - veloDrift_step_global[2].at(index);
        const double veloX = pit->v(0);
        const double veloY = pit->v(1);
        const double veloZ = pit->v(2);
        const double mass = pit->mass();
        const double epot = /* pit->U_pot() */ 0.0;
        const double ekinX = mass * veloCorrX * veloCorrX + (pit->U_rot_2()/3.0); //??? Wie Rotation?
        const double ekinY = mass * veloCorrY * veloCorrY + (pit->U_rot_2()/3.0); //??? Wie Rotation?
        const double ekinZ = mass * veloCorrZ * veloCorrZ + (pit->U_rot_2()/3.0); //??? Wie Rotation?
        ekin2_step.local.at(index) += ekinX + ekinY + ekinZ;
        ekin2Trans_step.local.at(index) += pit->U_kin();
        epot_step.local.at(index) += epot;
        ekin2Vect_step[0].local.at(index) += ekinX;
        ekin2Vect_step[1].local.at(index) += ekinY;
        ekin2Vect_step[2].local.at(index) += ekinZ;
        virialVect_step[0].local.at(index) += pit->Vi(0);
        virialVect_step[1].local.at(index) += pit->Vi(1);
        virialVect_step[2].local.at(index) += pit->Vi(2);
        energyfluxVect_step[0].local.at(index) += (pit->U_kin() + epot)*veloX + (pit->Vi(0)*veloX + pit->Vi(3)*veloY + pit->Vi(4)*veloZ);
        energyfluxVect_step[1].local.at(index) += (pit->U_kin() + epot)*veloY + (pit->Vi(6)*veloX + pit->Vi(1)*veloY + pit->Vi(5)*veloZ);
        energyfluxVect_step[2].local.at(index) += (pit->U_kin() + epot)*veloZ + (pit->Vi(7)*veloX + pit->Vi(8)*veloY + pit->Vi(2)*veloZ);

        if (!_singleComp) {
            const unsigned int cid = pit->componentid() + 1; // 0 represents all components
            const uint32_t indexCID = cid*_numBinsGlobal + index;
            ekin2_step.local.at(indexCID) += ekinX + ekinY + ekinZ;
            ekin2Trans_step.local.at(indexCID) += pit->U_kin();
            epot_step.local.at(indexCID) += epot;
            ekin2Vect_step[0].local.at(indexCID) += ekinX;
            ekin2Vect_step[1].local.at(indexCID) += ekinY;
            ekin2Vect_step[2].local.at(indexCID) += ekinZ;
            virialVect_step[0].local.at(indexCID) += pit->Vi(0);
            virialVect_step[1].local.at(indexCID) += pit->Vi(1);
            virialVect_step[2].local.at(indexCID) += pit->Vi(2);
            energyfluxVect_step[0].local.at(indexCID) += (pit->U_kin() + epot)*veloX + (pit->Vi(0)*veloX + pit->Vi(3)*veloY + pit->Vi(4)*veloZ);
            energyfluxVect_step[1].local.at(indexCID) += (pit->U_kin() + epot)*veloY + (pit->Vi(6)*veloX + pit->Vi(1)*veloY + pit->Vi(5)*veloZ);
            energyfluxVect_step[2].local.at(indexCID) += (pit->U_kin() + epot)*veloZ + (pit->Vi(7)*veloX + pit->Vi(8)*veloY + pit->Vi(2)*veloZ);
        }
    }

    // Calculate temperature (without drift) per bin over all processes
#ifdef ENABLE_MPI
    MPI_Allreduce(ekin2_step.local.data(), ekin2_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	for (unsigned long i = 0; i < _lenVector; i++) {
        ekin2_step.global.at(i) = ekin2_step.local.at(i);
	}
#endif

    for (unsigned long i = 0; i < _lenVector; i++) {
        const unsigned int cid = i/_numBinsGlobal;
        const unsigned int dof = _simulation.getEnsemble()->getComponent(cid)->getRotationalDegreesOfFreedom();
        const unsigned long numMols = numMolecules_step.global.at(i);
        if (numMols > 0ul) {
            temperature_step_global.at(i) = ekin2_step.global.at(i) / ((3 + dof)*numMols);
        }
        cout << i << " " << cid << " " << dof << " " << temperature_step_global.at(i) << endl;
    }


    // Calculate chemical potential
    if (_sampleChemPot) {
        // Calculate number of test particles per bin
        std::vector<double> dX(_numBinsGlobal, 0.0);
        std::vector<double> dY(_numBinsGlobal, 0.0);
        std::vector<double> dZ(_numBinsGlobal, 0.0);
        std::vector<unsigned long> nX(_numBinsGlobal, 0ul);
        std::vector<unsigned long> nY(_numBinsGlobal, 0ul);
        std::vector<unsigned long> nZ(_numBinsGlobal, 0ul);
        unsigned long nTestGlobal {0ul};

        for (uint16_t i = 0; i < _numBinsGlobal; i++) {
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
            const uint16_t idxStart = std::min(_numBinsGlobal, static_cast<uint16_t>(regionLowCorner[1]/_binwidth));
            double rY = regionLowCorner[1]+0.5*dY.at(idxStart);
            while (rY < regionHighCorner[1]) {
                const uint16_t index = std::min(_numBinsGlobal, static_cast<uint16_t>(rY/_binwidth));  // Index of bin

                double rX = regionLowCorner[0]+0.5*dX.at(idxStart);
                while (rX < regionHighCorner[0]) {

                    double rZ = regionLowCorner[2]+0.5*dZ.at(idxStart);
                    while (rZ < regionHighCorner[2]) {
                        if (temperature_step_global.at(index) > 1e-9) {
                            _mTest.setr(0,rX);
                            _mTest.setr(1,rY);
                            _mTest.setr(2,rZ);
                            const double deltaUpot = particleContainer->getEnergy(_particlePairsHandler, &_mTest, *_cellProcessor);
                            double chemPot = exp(-deltaUpot/temperature_step_global.at(index));
                            if (std::isfinite(chemPot)) {
                                chemPot_step.local.at(index) += chemPot;
                                countNTest_step.local.at(index)++;
    #ifndef NDEBUG
                                std::cout << "[ExtendedProfileSampling] Rank " << domainDecomp->getRank() << " : Inserting molecule at x,y,z = "
                                        << _mTest.r(0) << " , " << _mTest.r(1) << " , " << _mTest.r(2)
                                        << " ; chemPot = " << chemPot << " ; dU = " << deltaUpot << " ; T = " << temperature_step_global.at(index) << " ; index = " << index << std::endl;
    #endif
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
            Random* rnd = new Random();

            // Share of volume of present rank from whole domain
            const float domainShare = (regionSize[0]*regionSize[1]*regionSize[2])/(_globalBoxLength[0]*_globalBoxLength[1]*_globalBoxLength[2]); 
            const unsigned long nTest = static_cast<unsigned long>(domainShare*nTestGlobal);

    #if defined(_OPENMP)
            #pragma omp parallel
    #endif
            for (unsigned long i = 0; i < nTest; i++) {
                const double rX = regionLowCorner[0] + rnd->rnd()*regionSize[0];
                const double rY = regionLowCorner[1] + rnd->rnd()*regionSize[1];
                const double rZ = regionLowCorner[2] + rnd->rnd()*regionSize[2];
                const uint16_t index = std::min(_numBinsGlobal, static_cast<uint16_t>(rY/_binwidth));  // Index of bin
                if (temperature_step_global.at(index) > 1e-9) {
                    _mTest.setr(0,rX);
                    _mTest.setr(1,rY);
                    _mTest.setr(2,rZ);
                    const double deltaUpot = particleContainer->getEnergy(_particlePairsHandler, &_mTest, *_cellProcessor);
                    double chemPot = exp(-deltaUpot/temperature_step_global.at(index));
                    if (std::isfinite(chemPot)) {
    #if defined(_OPENMP)
                        #pragma omp atomic
    #endif
                        {
                        chemPot_step.local.at(index) += chemPot;
                        countNTest_step.local.at(index)++;
    #ifndef NDEBUG
                        std::cout << "[ExtendedProfileSampling] Rank " << domainDecomp->getRank() << " : Inserting molecule at x,y,z = "
                                << _mTest.r(0) << " , " << _mTest.r(1) << " , " << _mTest.r(2)
                                << " ; chemPot = " << chemPot << " ; dU = " << deltaUpot << " ; T = " << temperature_step_global.at(index) << " ; index = " << index << std::endl;
    #endif
                        }
                    }
                }
            }
        }
    }

    // Calculate further quantities
#ifdef ENABLE_MPI
    MPI_Reduce(epot_step.local.data(), epot_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[0].local.data(), virialVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[1].local.data(), virialVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[2].local.data(), virialVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin2Vect_step[0].local.data(), ekin2Vect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin2Vect_step[1].local.data(), ekin2Vect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin2Vect_step[2].local.data(), ekin2Vect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(energyfluxVect_step[0].local.data(), energyfluxVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(energyfluxVect_step[1].local.data(), energyfluxVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(energyfluxVect_step[2].local.data(), energyfluxVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(chemPot_step.local.data(), chemPot_step.global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(countNTest_step.local.data(), countNTest_step.global.data(), _numBinsGlobal, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	for (unsigned long i = 0; i < _lenVector; i++) {
        epot_step.global.at(i) = epot_step.local.at(i);
        virialVect_step[0].global.at(i) = virialVect_step[0].local.at(i);
        virialVect_step[1].global.at(i) = virialVect_step[1].local.at(i);
        virialVect_step[2].global.at(i) = virialVect_step[2].local.at(i);
        ekin2Vect_step[0].global.at(i) = ekin2Vect_step[0].local.at(i);
        ekin2Vect_step[1].global.at(i) = ekin2Vect_step[1].local.at(i);
        ekin2Vect_step[2].global.at(i) = ekin2Vect_step[2].local.at(i);
        energyfluxVect_step[0].global.at(i) = energyfluxVect_step[0].local.at(i);
        energyfluxVect_step[1].global.at(i) = energyfluxVect_step[1].local.at(i);
        energyfluxVect_step[2].global.at(i) = energyfluxVect_step[2].local.at(i);
        chemPot_step.global.at(i) = chemPot_step.local.at(i);
        countNTest_step.global.at(i) = countNTest_step.local.at(i);
	}
#endif

    // Only root calculates further quantities and accumulates them
    if (domainDecomp->getRank() == 0) {
        for (unsigned long i = 0; i < _lenVector; i++) {
            const unsigned long numMols = numMolecules_step.global.at(i);
            if (numMols > 0) {
                const unsigned int cid = i/_numBinsGlobal;
                const unsigned int dof = _simulation.getEnsemble()->getComponent(cid)->getRotationalDegreesOfFreedom();
                const double rho = numMols / _slabVolume;
                const double ViX = virialVect_step[0].global.at(i);
                const double ViY = virialVect_step[1].global.at(i);
                const double ViZ = virialVect_step[2].global.at(i);
                const double Tx = ekin2Vect_step[0].global.at(i) / ((1 + dof/3.0)*numMols); // ??? dof
                const double Ty = ekin2Vect_step[1].global.at(i) / ((1 + dof/3.0)*numMols); // ??? dof
                const double Tz = ekin2Vect_step[2].global.at(i) / ((1 + dof/3.0)*numMols); // ??? dof
                _numMolecules_accum.at(i)            += numMolecules_step.global.at(i);
                _density_accum.at(i)                 += rho;
                _temperature_accum.at(i)             += temperature_step_global.at(i); // ??? oder (Tx+Ty+Tz)/3.0
                cout << "EPS " << temperature_step_global.at(i) << " " << (Tx+Ty+Tz)/3.0 << endl;
                _ekin_accum.at(i)                    += ekin2Trans_step.global.at(i) / numMols;
                _epot_accum.at(i)                    += epot_step.global.at(i) / numMols;
                _pressure_accum.at(i)                += rho * ( (ViX + ViY + ViZ)/(3.0*numMols) + temperature_step_global.at(i) ); // ??? Welche Temperature? Mit Drift? Statisch/dynamisch?
                _chemPot_accum.at(i)                 += chemPot_step.global.at(i);
                _countNTest_accum.at(i)              += countNTest_step.global.at(i);

                _temperatureVect_accum[0].at(i)      += Tx;
                _temperatureVect_accum[1].at(i)      += Ty;
                _temperatureVect_accum[2].at(i)      += Tz;
                _velocityVect_accum[0].at(i)         += veloDrift_step_global[0].at(i);
                _velocityVect_accum[1].at(i)         += veloDrift_step_global[1].at(i);
                _velocityVect_accum[2].at(i)         += veloDrift_step_global[2].at(i);
                _pressureVect_accum[0].at(i)         += rho * ( ViX/numMols + Tx );
                _pressureVect_accum[1].at(i)         += rho * ( ViY/numMols + Ty );
                _pressureVect_accum[2].at(i)         += rho * ( ViZ/numMols + Tz );
                _energyfluxVect_accum[0].at(i)       += energyfluxVect_step[0].global.at(i) / numMols;
                _energyfluxVect_accum[1].at(i)       += energyfluxVect_step[1].global.at(i) / numMols;
                _energyfluxVect_accum[2].at(i)       += energyfluxVect_step[2].global.at(i) / numMols;

                _countSamples.at(i)++;
            }
        }
    }


    // Write out data every _writeFrequency step
    if ( (simstep - _startSampling) % _writeFrequency != 0 ) {
        return;
    }

    if (domainDecomp->getRank() == 0) {
        unsigned long numOutputs = (_singleComp) ? 1ul : (_simulation.getDomain()->getNumberOfComponents()+1);

        // Write output file
        std::stringstream ss;
        ss << std::setw(9) << std::setfill('0') << simstep;
        const std::string fname = "ExtendedProfileSampling_TS"+ss.str()+".dat";
        std::ofstream ofs;
        ofs.open(fname, std::ios::out);
        ofs << setw(24) << "pos";                           // Bin position
        for (unsigned long cid = 0; cid < numOutputs; cid++) {
            ofs << setw(24) << "numParts["<<cid<<"]"        // Average number of molecules in bin per step
                << setw(24) << "rho["<<cid<<"]"             // Density
                << setw(24) << "T["<<cid<<"]"               // Temperature without drift (i.e. "real" temperature)
                << setw(24) << "ekin["<<cid<<"]"            // Kinetic energy
                << setw(24) << "epot["<<cid<<"]"            // Potential energy
                << setw(24) << "p["<<cid<<"]"               // Pressure
                << setw(24) << "chemPot_res["<<cid<<"]"     // Chemical potential as known as mu_tilde (equals the ms2 value)
                << setw(24) << "numTest["<<cid<<"]"         // Number of inserted test particles per sample step
                << setw(24) << "T_x["<<cid<<"]"             // Temperature in x-direction
                << setw(24) << "T_y["<<cid<<"]"             // Temperature in y-direction
                << setw(24) << "T_z["<<cid<<"]"             // Temperature in z-direction
                << setw(24) << "v_x["<<cid<<"]"             // Drift velocity in x-direction
                << setw(24) << "v_y["<<cid<<"]"             // Drift velocity in y-direction
                << setw(24) << "v_z["<<cid<<"]"             // Drift velocity in z-direction
                << setw(24) << "p_x["<<cid<<"]"             // Pressure in x-direction
                << setw(24) << "p_y["<<cid<<"]"             // Pressure in y-direction
                << setw(24) << "p_z["<<cid<<"]"             // Pressure in z-direction
                << setw(24) << "jEF_x["<<cid<<"]"           // Energy flux in x-direction
                << setw(24) << "jEF_y["<<cid<<"]"           // Energy flux in y-direction
                << setw(24) << "jEF_z["<<cid<<"]"           // Energy flux in z-direction
                << setw(24) << "numSamples["<<cid<<"]";     // Number of samples (<= _writeFrequency)
        }
        ofs << std::endl;
        for (unsigned long i = 0; i < _lenVector; i++) {
            double numMolsPerStep {0.0}; // Not an int as particles change bin during simulation
            double rho {0.0};
            double T {0.0};
            double ekin {0.0};
            double epot {0.0};
            double p {0.0};
            double chemPot_res {0.0};
            double numTest {0.0};
            double T_x {0.0};
            double T_y {0.0};
            double T_z {0.0};
            double v_x {0.0};
            double v_y {0.0};
            double v_z {0.0};
            double p_x {0.0};
            double p_y {0.0};
            double p_z {0.0};
            double jEF_x {0.0};
            double jEF_y {0.0};
            double jEF_z {0.0};
            double numSamples {0.0};
            if (_countSamples.at(i) > 0ul) {
                numMolsPerStep = static_cast<double>(_numMolecules_accum.at(i))/_countSamples.at(i);
                rho         = _density_accum.at(i)           /_countSamples.at(i);
                T           = _temperature_accum.at(i)       /_countSamples.at(i);
                ekin        = _ekin_accum.at(i)              /_countSamples.at(i);
                epot        = _epot_accum.at(i)              /_countSamples.at(i);
                p           = _pressure_accum.at(i)          /_countSamples.at(i);
                chemPot_res = _chemPot_accum.at(i)           /_countSamples.at(i);
                numTest     = static_cast<double>(_countNTest_accum.at(i))/_countSamples.at(i);
                T_x         = _temperatureVect_accum[0].at(i)/_countSamples.at(i);
                T_y         = _temperatureVect_accum[1].at(i)/_countSamples.at(i);
                T_z         = _temperatureVect_accum[2].at(i)/_countSamples.at(i);
                v_x         = _velocityVect_accum[0].at(i)   /_countSamples.at(i);
                v_y         = _velocityVect_accum[1].at(i)   /_countSamples.at(i);
                v_z         = _velocityVect_accum[2].at(i)   /_countSamples.at(i);
                p_x         = _pressureVect_accum[0].at(i)   /_countSamples.at(i);
                p_y         = _pressureVect_accum[1].at(i)   /_countSamples.at(i);
                p_z         = _pressureVect_accum[2].at(i)   /_countSamples.at(i);
                jEF_x       = _energyfluxVect_accum[0].at(i) /_countSamples.at(i);
                jEF_y       = _energyfluxVect_accum[1].at(i) /_countSamples.at(i);
                jEF_z       = _energyfluxVect_accum[2].at(i) /_countSamples.at(i);
                numSamples  = _countSamples.at(i);
            }
            if ((_chemPot_accum.at(i) > 0.0) and (_countNTest_accum.at(i) > 0ul) and (_countSamples.at(i) > 0ul)) {
                numTest     = _countNTest_accum.at(i)/_countSamples.at(i);
                chemPot_res = -log(_chemPot_accum.at(i)/_countNTest_accum.at(i)) + log(rho);
            }
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << (i+0.5)*_binwidth;
            for (unsigned long cid = 0; cid < numOutputs; cid++) {
                ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                    << numMolsPerStep
                    << rho
                    << T
                    << ekin
                    << epot
                    << p
                    << chemPot_res
                    << numTest
                    << T_x
                    << T_y
                    << T_z
                    << v_x
                    << v_y
                    << v_z
                    << p_x
                    << p_y
                    << p_z
                    << jEF_x
                    << jEF_y
                    << jEF_z
                    << numSamples;
            }
            ofs << std::endl;
        }
        ofs.close();
    }

    // Reset vectors to zero
    resetVectors();
}

// Resize vectors
void ExtendedProfileSampling::resizeVectors() {
    _numMolecules_accum.resize(_lenVector);
    _density_accum.resize(_lenVector);
    _temperature_accum.resize(_lenVector);
    _ekin_accum.resize(_lenVector);
    _epot_accum.resize(_lenVector);
    _pressure_accum.resize(_lenVector);
    _chemPot_accum.resize(_lenVector);
    _countNTest_accum.resize(_lenVector);
    
    for (unsigned short d = 0; d < 3; d++) {
        _temperatureVect_accum.at(d).resize(_lenVector);
        _velocityVect_accum.at(d).resize(_lenVector);
        _pressureVect_accum.at(d).resize(_lenVector);
        _energyfluxVect_accum.at(d).resize(_lenVector);
    }

    _countSamples.resize(_lenVector);
}

// Fill vectors with zeros
void ExtendedProfileSampling::resetVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    std::fill(_density_accum.begin(), _density_accum.end(), 0.0f);
    std::fill(_temperature_accum.begin(), _temperature_accum.end(), 0.0f);
    std::fill(_ekin_accum.begin(), _ekin_accum.end(), 0.0f);
    std::fill(_epot_accum.begin(), _epot_accum.end(), 0.0f);
    std::fill(_pressure_accum.begin(), _pressure_accum.end(), 0.0f);
    std::fill(_chemPot_accum.begin(), _chemPot_accum.end(), 0.0f);
    std::fill(_countNTest_accum.begin(), _countNTest_accum.end(), 0ul);
    
    for (unsigned short d = 0; d < 3; d++) {
        std::fill(_temperatureVect_accum.at(d).begin(), _temperatureVect_accum.at(d).end(), 0.0f);
        std::fill(_velocityVect_accum.at(d).begin(), _velocityVect_accum.at(d).end(), 0.0f);
        std::fill(_pressureVect_accum.at(d).begin(), _pressureVect_accum.at(d).end(), 0.0f);
        std::fill(_energyfluxVect_accum.at(d).begin(), _energyfluxVect_accum.at(d).end(), 0.0f);
    }

    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);
}
