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


ChemPotSampling::ChemPotSampling(){}

ChemPotSampling::~ChemPotSampling(){}

void ChemPotSampling::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
    global_log -> info() << "[ChemPotSampling] Init" << std::endl;

    _globalBoxLengthY = domain->getGlobalLength(1);

    _numBinsGlobal = static_cast<uint16_t>(_globalBoxLengthY/_binwidth);
    if (_globalBoxLengthY/_binwidth != static_cast<float>(_numBinsGlobal)) {
        global_log -> error() << "[ChemPotSampling] Can not divide domain without remainder! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }
    _slabVolume = domain->getGlobalLength(0)*domain->getGlobalLength(2)*_binwidth;

    _chemPotSum.local.resize(_numBinsGlobal);
    _temperatureSum.local.resize(_numBinsGlobal);
    _numMoleculesSum.local.resize(_numBinsGlobal);

    _chemPotSum.global.resize(_numBinsGlobal);
    _temperatureSum.global.resize(_numBinsGlobal);
    _numMoleculesSum.global.resize(_numBinsGlobal);

    _cellProcessor = _simulation.getCellProcessor();
    _mTest = Molecule(std::numeric_limits<unsigned long>::max()-domainDecomp->getRank(), &_simulation.getEnsemble()->getComponents()->at(0));

    resetVectors();

}

void ChemPotSampling::readXML(XMLfileUnits& xmlconfig){

    xmlconfig.getNodeValue("binwidth", _binwidth);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("stop", _stopSampling);

    global_log -> info() << "[ChemPotSampling] Start:Freq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    global_log -> info() << "[ChemPotSampling] Binwidth: " << _binwidth << std::endl;

}

void ChemPotSampling::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                         unsigned long simstep) {
    
    double regionLowCorner[3], regionHighCorner[3];
    ParticlePairs2PotForceAdapter particlePairsHandler(*domain);

    double mass = 1.0;

    // Accumulated over one sampling step
    CommVar<std::vector<double>> ekin2;
    std::array<CommVar<std::vector<double>>, 3> velocity; // Array for x,y,z values (local and global) including vectors for bins
    CommVar<std::vector<unsigned long>> numMols;

    ekin2.local.resize(_numBinsGlobal);
    numMols.local.resize(_numBinsGlobal);
    velocity[0].local.resize(_numBinsGlobal);
    velocity[1].local.resize(_numBinsGlobal);
    velocity[2].local.resize(_numBinsGlobal);

    ekin2.global.resize(_numBinsGlobal);
    numMols.global.resize(_numBinsGlobal);
    velocity[0].global.resize(_numBinsGlobal);
    velocity[1].global.resize(_numBinsGlobal);
    velocity[2].global.resize(_numBinsGlobal);

    for (unsigned d = 0; d < 3; d++) {
        regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
    }

    // Sample temperature and number of molecules
    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        // sample density profile
        double ry = pit->r(1);
        uint16_t index = static_cast<uint16_t>(ry/_binwidth);  // Index of bin

        _numMoleculesSum.local.at(index) += 1;

        double v[3];
        v[0] = pit->v(0);
        v[1] = pit->v(1);
        v[2] = pit->v(2);

        ekin2.local.at(index) += pit->U_trans_2() + pit->U_rot_2();

        velocity[0].local.at(index) += v[0];
        velocity[1].local.at(index) += v[1];
        velocity[2].local.at(index) += v[2];

        numMols.local.at(index) += 1;
    }

#ifdef ENABLE_MPI
	MPI_Reduce(numMols.local.data(), numMols.global.data(), _numBinsGlobal, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin2.local.data(), ekin2.global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(velocity[0].local.data(), velocity[0].global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(velocity[1].local.data(), velocity[1].global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(velocity[2].local.data(), velocity[2].global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	for (unsigned int i = 0; i < _numBinsGlobal; i++) {
		numMols.global.at(i) = numMols.local.at(i);
        ekin2.global.at(i) = ekin2.local.at(i);
		velocity[0].global.at(i) = velocity[0].local.at(i);
        velocity[1].global.at(i) = velocity[1].local.at(i);
        velocity[2].global.at(i) = velocity[2].local.at(i);
	}
#endif

    std::vector<double> temperatureStep;
    temperatureStep.resize(_numBinsGlobal);
    for (uint16_t i = 0; i < _numBinsGlobal; i++) {
        double vDrift2 = (velocity[0].global.at(i) * velocity[0].global.at(i)
                       + velocity[1].global.at(i) * velocity[1].global.at(i)
                       + velocity[2].global.at(i) * velocity[2].global.at(i))/numMols.global.at(i);
        // Temperature in the present sampling step
        cout << "velo  " << velocity[0].global.at(i) << " " << velocity[1].global.at(i) << " " << velocity[2].global.at(i) << std::endl;
        temperatureStep.at(i) = (ekin2.global.at(i) - vDrift2*vDrift2*mass)/(numMols.global.at(i)*3);
        cout << "vDrift2  " << vDrift2 << " ekin2 " << ekin2.global.at(i) << " 3*numMols " << numMols.global.at(i)*3 << " T " << temperatureStep.at(i) << std::endl;
    }

cout << "temperatureStep ";
for (auto ele : temperatureStep) {
    cout << ele << " ";
}
cout << std::endl;

    // Insert particles in lattice structure and sample chem. pot.
    double dX = 1.0;
    double dZ = dX;
    double dY = 1.0;
/*
    for (double iY = regionLowCorner[1]; iY < regionHighCorner[1]; iY += dY) {
        uint16_t index = static_cast<uint16_t>(iY/_binwidth);  // Index of bin

        for (double iX = regionLowCorner[0]; iX < regionHighCorner[0]; iX += dX) {
            for (double iZ = regionLowCorner[2]; iZ < regionHighCorner[2]; iZ += dZ) {
                _mTest.setr(0,iX);
                _mTest.setr(1,iY);
                _mTest.setr(2,iZ);
                double deltaUpot = particleContainer->getEnergy(&particlePairsHandler, &_mTest, *_cellProcessor);
                // cout << "Inserting molecule at " << _mTest.r(0) << " " << _mTest.r(1) << " " << _mTest.r(2) << " ; dU = " << deltaUpot << " index " << index << std::endl;
                _chemPotSum.local.at(index) += exp(-deltaUpot/temperatureStep.at(index));
            }
        }

    }
*/
    cout << "CPS 6" << std::endl;

    // sampling starts after initial timestep and with respect to write frequency
    if ( (simstep - _startSampling) % _writeFrequency != 0 ) {
        return;
    }
	if ( (simstep <= _startSampling) or (simstep > _stopSampling)) {
        return;
    }
    // do not write data directly after (re)start
	if ( simstep == global_simulation->getNumInitTimesteps() ) {
        return;
    }

#ifdef ENABLE_MPI
	MPI_Reduce(_numMoleculesSum.local.data(), _numMoleculesSum.global.data(), _numBinsGlobal, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_temperatureSum.local.data(), _temperatureSum.global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_chemPotSum.local.data(), _chemPotSum.global.data(), _numBinsGlobal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	for (unsigned int i = 0; i < _numBinsGlobal; i++) {
		_numMoleculesSum.global.at(i) = _numMoleculesSum.local.at(i);
		_temperatureSum.global.at(i) = _temperatureSum.local.at(i);
        _chemPotSum.global.at(i) = _chemPotSum.local.at(i);
	}
#endif

cout << "CPS 7" << std::endl;

    if (domainDecomp->getRank() != 0) {
        // write output
        std::stringstream ss;
        ss << std::setw(9) << std::setfill('0') << simstep;
        const std::string fname = "ChemPotSampling_TS"+ss.str()+".dat";
        std::ofstream ofs;
        ofs.open(fname, std::ios::out);
        ofs << setw(12) << "pos" << setw(12) << "numMols" << setw(12) << "density" << setw(12) << "temperature" << setw(12) << "chemPot" << std::endl;
        for (uint16_t i = 0; i < _numBinsGlobal; i++) {
            double T = _temperatureSum.global.at(i)/_writeFrequency;
            double chemPot = -T*log(_chemPotSum.global.at(i)/_writeFrequency);
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << (i+0.5)*_binwidth;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _numMoleculesSum.global.at(i)/_writeFrequency;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _numMoleculesSum.global.at(i)/(_writeFrequency*_slabVolume);
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << T;
            ofs << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << chemPot;
        }
        ofs.close();
    }

    cout << "CPS 8" << std::endl;

    // reset vectors to zero
    resetVectors();

    cout << "CPS 9" << std::endl;
}


// Filling vectors with zeros
void ChemPotSampling::resetVectors() {
    std::fill(_chemPotSum.local.begin(), _chemPotSum.local.end(), 0);
    std::fill(_temperatureSum.local.begin(), _temperatureSum.local.end(), 0);
    std::fill(_numMoleculesSum.local.begin(), _numMoleculesSum.local.end(), 0);   

    std::fill(_chemPotSum.global.begin(), _chemPotSum.global.end(), 0);
    std::fill(_temperatureSum.global.begin(), _temperatureSum.global.end(), 0);
    std::fill(_numMoleculesSum.global.begin(), _numMoleculesSum.global.end(), 0);
}