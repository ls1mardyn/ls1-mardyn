/*
 * ChemPotSampling.cpp
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#include "ChemPotSampling.h"

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

#include <math.h>


ChemPotSampling::ChemPotSampling()
    :
    _binwidth(1.0f),
    _factorNumTest(4.0f),
    _startSampling(0ul),
    _writeFrequency(5000ul),
    _stopSampling(1000000000ul)
{}

ChemPotSampling::~ChemPotSampling() {}

void ChemPotSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    _numBinsGlobal = static_cast<uint16_t>(_globalBoxLength[1]/_binwidth);
    if (_globalBoxLength[1]/_binwidth != static_cast<float>(_numBinsGlobal)) {
        global_log->error() << "[ChemPotSampling] Can not divide domain without remainder! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }
    _slabVolume = domain->getGlobalLength(0)*domain->getGlobalLength(2)*_binwidth;

    _chemPotSum.local.resize(_numBinsGlobal);
    _chemPotSum.global.resize(_numBinsGlobal);

    _temperatureSumGlobal.resize(_numBinsGlobal);
    _temperatureWithDriftSumGlobal.resize(_numBinsGlobal);
    _numMoleculesSumGlobal.resize(_numBinsGlobal);

    // CHANGE OF CELL PROCESSOR DURING SIMULATION???
    _cellProcessor = _simulation.getCellProcessor();
    // CHANGE OF PP HANDLER DURING SIMULATION???
    _particlePairsHandler = new ParticlePairs2PotForceAdapter(*domain);
    // MolID is maximum possible number minus rank to prevent duplicate IDs
    // Always insert molecule of first component
    const unsigned long molID = std::numeric_limits<unsigned long>::max()-static_cast<unsigned long>(domainDecomp->getRank());
    _mTest = Molecule(molID, &_simulation.getEnsemble()->getComponents()->at(0));

    resetVectors();
}

void ChemPotSampling::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("binwidth", _binwidth);  // Default: 1.0
    xmlconfig.getNodeValue("numTest", _factorNumTest);  // Default: 4.0
    xmlconfig.getNodeValue("start", _startSampling);  // Default: 0
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);  // Default: 10000
    xmlconfig.getNodeValue("stop", _stopSampling);  // Default: 1000000000

    global_log->info() << "[ChemPotSampling] Start:Freq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    global_log->info() << "[ChemPotSampling] Binwidth: " << _binwidth << std::endl;
}

void ChemPotSampling::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* /* domain */,
                         unsigned long simstep) {

    // sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep <= _startSampling) or (simstep > _stopSampling)) {
        return;
    }
    // do not write data directly after (re)start
	if ( simstep == global_simulation->getNumInitTimesteps() ) {
        return;
    }
    
    double regionLowCorner[3], regionHighCorner[3];

    // Accumulated over one sampling step
    CommVar<std::vector<double>> ekin2;
    std::array<CommVar<std::vector<double>>,3> velocity;
    CommVar<std::vector<unsigned long>> numMols;

    ekin2.local.resize(_numBinsGlobal);
    numMols.local.resize(_numBinsGlobal);
    ekin2.global.resize(_numBinsGlobal);
    numMols.global.resize(_numBinsGlobal);

    for (int d = 0; d < 3; d++) {
        velocity[d].local.resize(_numBinsGlobal);
        velocity[d].global.resize(_numBinsGlobal);
        std::fill(velocity[d].local.begin(), velocity[d].local.end(), 0.0);
        std::fill(velocity[d].global.begin(), velocity[d].global.end(), 0.0);
    }

    std::fill(ekin2.local.begin(), ekin2.local.end(), 0.0);
    std::fill(numMols.local.begin(), numMols.local.end(), 0);
    std::fill(ekin2.global.begin(), ekin2.global.end(), 0.0);
    std::fill(numMols.global.begin(), numMols.global.end(), 0);

    for (int d = 0; d < 3; d++) {
        regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
    }

    // Sample temperature and number of molecules
    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        // sample density profile
        double ry = pit->r(1);
        uint16_t index = std::min(_numBinsGlobal, static_cast<uint16_t>(ry/_binwidth));  // Index of bin

        numMols.local.at(index) += 1;
        ekin2.local.at(index) += pit->U_trans_2() + pit->U_rot_2();

        velocity[0].local.at(index) += pit->v(0);
        velocity[1].local.at(index) += pit->v(1);
        velocity[2].local.at(index) += pit->v(2);
    }

#ifdef ENABLE_MPI
	MPI_Allreduce(numMols.local.data(), numMols.global.data(), _numBinsGlobal, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(ekin2.local.data(), ekin2.global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocity[0].local.data(), velocity[0].global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocity[1].local.data(), velocity[1].global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocity[2].local.data(), velocity[2].global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	for (unsigned int i = 0; i < _numBinsGlobal; i++) {
		numMols.global.at(i) = numMols.local.at(i);
        ekin2.global.at(i) = ekin2.local.at(i);
        velocity[0].global.at(i) = velocity[0].local.at(i);
        velocity[1].global.at(i) = velocity[1].local.at(i);
        velocity[2].global.at(i) = velocity[2].local.at(i);
	}
#endif

    std::vector<double> temperatureStep(_numBinsGlobal, 0.0);
    for (uint16_t i = 0; i < _numBinsGlobal; i++) {
        double veloDrift2 = velocity[0].global.at(i) * velocity[0].global.at(i)
                                   + velocity[1].global.at(i) * velocity[1].global.at(i)
                                   + velocity[2].global.at(i) * velocity[2].global.at(i);
        const double ekin2_T = ekin2.global.at(i) - veloDrift2/numMols.global.at(i);
        temperatureStep.at(i) = ekin2_T/(numMols.global.at(i)*3.0);
        _temperatureSumGlobal.at(i) += temperatureStep.at(i);
        _temperatureWithDriftSumGlobal.at(i) += ekin2.global.at(i)/(numMols.global.at(i)*3.0);
        _numMoleculesSumGlobal.at(i) += numMols.global.at(i);
    }

    // Insert particles in lattice structure and sample chem. pot.
    std::vector<double> dX(_numBinsGlobal, 0.0);
    std::vector<double> dY(_numBinsGlobal, 0.0);
    std::vector<double> dZ(_numBinsGlobal, 0.0);
    for (uint16_t i = 0; i < _numBinsGlobal; i++) {
        const unsigned long nTest = _factorNumTest*numMols.global.at(i);
        const unsigned long nY = std::max(1.0,std::pow((nTest*_binwidth*_binwidth)/(_globalBoxLength[0]*_globalBoxLength[0]),(1./3.)));
        dY.at(i) = _binwidth/nY;
        const unsigned long nX = std::max(1.0,std::pow((nTest*_globalBoxLength[0]*_globalBoxLength[0])/(_binwidth*_globalBoxLength[0]),(1./3.)));
        dX.at(i) = _globalBoxLength[0]/nX;
        const unsigned long nZ = nTest/(nX*nY);
        dZ.at(i) = _globalBoxLength[2]/nZ;
    }

    // Index of bin in which the left region boundary (y-dir) is in
    uint16_t idxStart = std::min(_numBinsGlobal, static_cast<uint16_t>(regionLowCorner[1]/_binwidth));
    double rY = regionLowCorner[1]+0.5*dY.at(idxStart);
    while (rY < regionHighCorner[1]) {
        uint16_t index = std::min(_numBinsGlobal, static_cast<uint16_t>(rY/_binwidth));  // Index of bin

        double rX = regionLowCorner[0]+0.5*dX.at(idxStart);
        while (rX < regionHighCorner[0]) {

            double rZ = regionLowCorner[2]+0.5*dZ.at(idxStart);
            while (rZ < regionHighCorner[2]) {
                _mTest.setr(0,rX);
                _mTest.setr(1,rY);
                _mTest.setr(2,rZ);
                double deltaUpot = particleContainer->getEnergy(_particlePairsHandler, &_mTest, *_cellProcessor);
                global_log->debug() << "Inserting molecule at x,y,z = " << _mTest.r(0) << " , " << _mTest.r(1) << " , " << _mTest.r(2)
                                    << " ; dU = " << deltaUpot << " ; index = " << index << std::endl;
                _chemPotSum.local.at(index) += exp(-deltaUpot/temperatureStep.at(index));
                rZ += dZ.at(index);
            }
            rX += dX.at(index);
        }
        rY += dY.at(index);
    }

    // Write out data
    if ( (simstep - _startSampling) % _writeFrequency != 0 ) {
        return;
    }

#ifdef ENABLE_MPI
	MPI_Reduce(_chemPotSum.local.data(), _chemPotSum.global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	for (unsigned int i = 0; i < _numBinsGlobal; i++) {
        _chemPotSum.global.at(i) = _chemPotSum.local.at(i);
	}
#endif

    if (domainDecomp->getRank() == 0) {
        // write output
        std::stringstream ss;
        ss << std::setw(9) << std::setfill('0') << simstep;
        const std::string fname = "ChemPotSampling_TS"+ss.str()+".dat";
        std::ofstream ofs;
        ofs.open(fname, std::ios::out);
        ofs << setw(24) << "pos" << setw(24) << "numMols" << setw(24) << "density" << setw(24) << "temperature" << setw(24) << "temp_with_Drift" << setw(24) << "chemPot" << std::endl;
        for (uint16_t i = 0; i < _numBinsGlobal; i++) {
            double T = _temperatureSumGlobal.at(i)/_writeFrequency;
            double chemPot = -T*log(_chemPotSum.global.at(i)/_writeFrequency);
            double numMolsPerStep = static_cast<double>(_numMoleculesSumGlobal.at(i))/_writeFrequency; // Not an int as particles change bin during simulation
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << (i+0.5)*_binwidth;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << numMolsPerStep;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << numMolsPerStep/_slabVolume;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << T;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _temperatureWithDriftSumGlobal.at(i)/_writeFrequency;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << chemPot << std::endl;
        }
        ofs.close();
    }

    // reset vectors to zero
    resetVectors();
}


// Filling vectors with zeros
void ChemPotSampling::resetVectors() {
    std::fill(_chemPotSum.local.begin(), _chemPotSum.local.end(), 0.0);
    std::fill(_chemPotSum.global.begin(), _chemPotSum.global.end(), 0.0);

    std::fill(_temperatureSumGlobal.begin(), _temperatureSumGlobal.end(), 0.0);
    std::fill(_temperatureWithDriftSumGlobal.begin(), _temperatureWithDriftSumGlobal.end(), 0.0);
    std::fill(_numMoleculesSumGlobal.begin(), _numMoleculesSumGlobal.end(), 0);
}