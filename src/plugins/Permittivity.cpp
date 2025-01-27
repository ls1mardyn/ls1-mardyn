/*
 * Permittivity.h
 *
 *  Created on: November 2019
 *      Author: Joshua Marx
 */

 // DESCRIPTION: Samples the relative permittivity of Stockmayer fluids in the NVT ensemble
 // Important: Always plot the running average data to make sure convergence has been achieved. Permittivity may take a long time to converge, i.e. a few million steps with ~1000 particles. Reducing number of slabs for the thermostat can drastically improve results/convergence!
 // If a simulation is resumed from a restart file, then the existing running average file is ammended but the computation of the running averages starts anew at the time step of the restart file
#include "Permittivity.h"
void Permittivity::readXML(XMLfileUnits& xmlconfig) {
	Log::global_log->info() << "Calculation of relative permittivity enabled." << std::endl;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	xmlconfig.getNodeValue("initstatistics", _initStatistics);
	xmlconfig.getNodeValue("recordingtimesteps", _recordingTimesteps);
	xmlconfig.getNodeValue("runningaveragestep", _runningAverageSteps);
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
}

void Permittivity::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	_totalNumTimeSteps = _simulation.getNumTimesteps();
	Log::global_log->warning() << "Permittivity is being sampled. Using slab thermostat is not recommended" << std::endl;
	Log::global_log->info() << "Number of blocks for calculation of permittivity: "<<_numOutputs << std::endl;
	_readStartingStep = false;
	_currentOutputNum = 0;
	_accumulatedSteps = 1;
	_numParticlesLocal = 0;
	_totalAverageSquaredM = 0.;
	_MSquaredRAV = 0.;
	_RAVCounter = 0;

	std::string permName = _outputPrefix + ".permRAV";
	int mpi_rank = domainDecomp->getRank();
	if (mpi_rank == 0) {
		_ravStream.open(permName.c_str(), std::ios::app);
		_ravStream.precision(7);
	}
	_V = domain->getGlobalLength(0) * domain->getGlobalLength(1) * domain->getGlobalLength(2);
	_T = domain->getTargetTemperature(0);

	for (unsigned long i = 0; i<_writeFrequency;i++){
		for (unsigned int j = 0; j<3; j++){
			_localM[i][j] = 0.;
			_globalM[i][j] = 0.;
		}
	}

	for (unsigned int i = 0; i < _numOutputs; i++) {
		_outputSquaredM[i] = 0.;
		_numParticles[i] = 0;
		for (unsigned int j = 0; j < 3; j++) {
			_outputM[i][j] = 0.;
		}
	}

	for (unsigned int i = 0; i < 3; i++) {
		_totalAverageM[i] = 0.;
		_MRAV[i] = 0.;
	}

	std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
	_numComponents = components->size();
	std::vector<unsigned int> isDipole(_numComponents);

	for (unsigned int i = 0; i < _numComponents; ++i) {
		Component& ci = (*components)[i];
		isDipole[i] = ci.numDipoles();

		if (isDipole[i] == 1) {
			bool orientationIsCorrect = ci.dipole(0).e() == std::array<double, 3>{0,0,1};
			_myAbs[i] = ci.dipole(0).abs();
			if(not orientationIsCorrect){
				Log::global_log->error() << "Wrong dipole vector chosen! Please always choose [eMyx eMyy eMyz] = [0 0 1] when using the permittivity plugin" << std::endl;
			}
		}
	}
}
void Permittivity::record(ParticleContainer* particleContainer) {

	double Quaternion[4];
	double orientationVector[3];

	for(ParticleIterator tempMol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol){

		Quaternion[0] = tempMol->q().qw();
		Quaternion[1] = tempMol->q().qx();
		Quaternion[2] = tempMol->q().qy();
		Quaternion[3] = tempMol->q().qz();

		// Calculates dipole moment vector from quaternions
		orientationVector[0] = 2 * (Quaternion[1] * Quaternion[3] + Quaternion[0] * Quaternion[2]);
		orientationVector[1] = 2 * (Quaternion[2] * Quaternion[3] - Quaternion[0] * Quaternion[1]);
		orientationVector[2] = 1 - 2 * (Quaternion[1] * Quaternion[1] + Quaternion[2] * Quaternion[2]);

		_numParticlesLocal++;

		for (unsigned int i = 0; i < 3; i++) {
			// Calculates M as sum of all molecular dipole moment vectors for current time step
			_localM[_accumulatedSteps][i] += orientationVector[i] * _myAbs[tempMol->getComponentLookUpID()];
		}
	}
}

void Permittivity::writeRunningAverage(unsigned long indexM, double tempMX, double tempMY, double tempMZ, double tempMSquared) { // writes instantaneous values, temporary averages and running averages of Mx, My, Mz, <M2> and epsilon
	double timestep;
	if (_startingStep <= 1) {
		timestep = _initStatistics + _RAVCounter * _runningAverageSteps;
	}
	else {
		timestep = _startingStep + _RAVCounter * _runningAverageSteps;
	}
	double M[3];
	unsigned long RAVsteps = _RAVCounter * _runningAverageSteps;
	double Msquared = _MSquaredRAV / (double)RAVsteps;
	double MsquaredInst = 0.;

	for (unsigned int i = 0; i < 3; i++) {
		M[i] = _MRAV[i] / (double)RAVsteps;
		MsquaredInst += _globalM[indexM][i] * _globalM[indexM][i];
	}

	double perm = 1. + 4. * M_PI * Msquared / (3. * _T * _V);
	double permInst = 1. + 4. * M_PI * MsquaredInst / (3. * _T * _V);
	double permTemp = 1. + 4. * M_PI * tempMSquared / (_runningAverageSteps * 3. * _T * _V);

	_ravStream << timestep << "\t" << RAVsteps << "\t" << _globalM[indexM][0] << "\t" << _globalM[indexM][1] << "\t" << _globalM[indexM][2] << "\t" << MsquaredInst << "\t" << permInst << "\t" << tempMX / (double)(_runningAverageSteps) << "\t" << tempMY / (double)(_runningAverageSteps) << "\t" << tempMZ / (double)(_runningAverageSteps) << "\t"
	<< tempMSquared / (double)(_runningAverageSteps) << "\t" << permTemp << "\t" << M[0] << "\t" << M[1] << "\t" << M[2] << "\t" << Msquared << "\t" << perm << std::endl;
}
void Permittivity::collect(DomainDecompBase* domainDecomp) {

	//Calculate global M for each time step of the current block from local values
	domainDecomp->collCommInit(3 * (_accumulatedSteps-1));
	for (unsigned long j = 0; j < (_accumulatedSteps-1); j++){
		for (unsigned i = 0; i < 3; i++) {
			domainDecomp->collCommAppendDouble(_localM[j][i]);
		}
	}
	domainDecomp->collCommAllreduceSum();

	for (unsigned long j = 1; j < _accumulatedSteps; j++){
		for (unsigned long i = 0; i < 3; i++) {
			_globalM[j][i] = domainDecomp->collCommGetDouble();
		}
	}

	domainDecomp->collCommFinalize();

	double tempM[3];
	double tempMSquared = 0;
	for (unsigned int i=1; i < 3; i++) {
		tempM[i] = 0.;
	}

	// Calculates ensemble averages <M²> and <M> for current block
	for (unsigned long i = 1; i < _accumulatedSteps; i++) {

		tempMSquared += _globalM[i][0] * _globalM[i][0] + _globalM[i][1] * _globalM[i][1] + _globalM[i][2] * _globalM[i][2];
		_MSquaredRAV += _globalM[i][0] * _globalM[i][0] + _globalM[i][1] * _globalM[i][1] + _globalM[i][2] * _globalM[i][2];
		_outputSquaredM[_currentOutputNum] += _globalM[i][0] * _globalM[i][0] + _globalM[i][1] * _globalM[i][1] + _globalM[i][2] * _globalM[i][2];

		for (unsigned int j = 0; j < 3; j++) {

			_outputM[_currentOutputNum][j] += _globalM[i][j];
			_MRAV[j] += _globalM[i][j];
			tempM[j] += _globalM[i][j];
		}

		if (i % _runningAverageSteps == 0 && i > 0) {

			_RAVCounter++;
			writeRunningAverage(i, tempM[0], tempM[1], tempM[2], tempMSquared);
			tempMSquared = 0.;

			for (unsigned int j=0; j < 3; j++) {

				tempM[j] = 0.;
			}
		}
	}

	for (unsigned int i = 0; i < 3; i++) {
		_outputM[_currentOutputNum][i] /= ((double)_accumulatedSteps-1);
	}

	_outputSquaredM[_currentOutputNum] /= ((double)_accumulatedSteps-1);
#ifdef ENABLE_PERSISTENT
	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), _numParticlesLocal);
	collComm.persistent();
	collComm.get(_numParticles[_currentOutputNum]);
#else
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(_numParticlesLocal);
	domainDecomp->collCommAllreduceSum();
	_numParticles[_currentOutputNum] = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
#endif
	_currentOutputNum++;
}

void Permittivity::reset() {
	for (unsigned long i = 1; i<_accumulatedSteps;i++){
		for (unsigned j = 0; j<3; j++){
			_localM[i][j] = 0;
			_globalM[i][j] = 0;
		}
	}
	_accumulatedSteps = 1;
	_numParticlesLocal = 0;
}

void Permittivity::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {

	if (_readStartingStep == false) {  // this if-statement is used to differentiate between a simulation that is newly started and one that is resumed from a restart file. For restarted simulations, the running average file from
		_startingStep = simstep;       // the original simulation is then simply continued
		if (_startingStep <= 1) {
			_numOutputs = (unsigned int)ceil(((double)_totalNumTimeSteps - (double)_initStatistics)/ (double)_writeFrequency);
			_ravStream << "//Output generated by Permittivity plugin\n" << std::endl;
			_ravStream << "time steps\trecording steps\tMx_inst\tMy_inst\tMz_inst\tMsquared_inst\tperm_inst\tMx_temp\tMy_temp\tMz_temp\tMsquared_temp\tperm_temp\tMx_rav\tMy_rav\tMz_rav\tMsquared_rav\tperm_rav\t" << std::endl;
		}
		else {
			_numOutputs = (unsigned int)ceil(((double)_totalNumTimeSteps - (double)_startingStep)/ (double)_writeFrequency);
		}
		_readStartingStep = true;
		return;
	}

	if (simstep > _initStatistics && simstep % _recordingTimesteps == 0) {
		record(particleContainer);
		_accumulatedSteps++;
	}

	if (simstep > _initStatistics && simstep % _writeFrequency == 0) {
		collect(domainDecomp);
		reset();
	}

	int mpi_rank = domainDecomp->getRank();
	if (simstep == _totalNumTimeSteps) {
		if (mpi_rank == 0){
			output(domain, simstep);
		}
	}
}

void Permittivity::output(Domain* domain, unsigned long timestep) {
	double writeInit;
	double averagePermittivity = 0.;
	double averagePermittivity2 = 0.;
	double blockPermittivity = 0.;
	double blockPermittivity2 = 0.;
	double correctionFirstBlock;
	unsigned long numSteps;
	if (_initStatistics % _writeFrequency == 0 || _startingStep > 1) {
		correctionFirstBlock = 0.;
	}
	else if (_initStatistics < _writeFrequency) {
		correctionFirstBlock = (double)_initStatistics/(double)_writeFrequency;
	}
	else {
		correctionFirstBlock = (double)(_initStatistics % _writeFrequency) / (double)_writeFrequency;
	}
	if (_startingStep > 1) {
		writeInit = 0.;
	}
	else {
		writeInit = floor((double)_initStatistics / (double)_writeFrequency);
	}
	std::string prefix;
	std::ostringstream osstrm;
	osstrm << _outputPrefix;
	osstrm << std::right << _startingStep;
	prefix = osstrm.str();
	osstrm.str("");
	osstrm.clear();
	std::string permName = prefix + ".perm";
	std::ofstream writer(permName.c_str());
	writer.precision(7);

	writer << "//Output generated by Permittivity plugin\n"
		   << "//Two values for the permittivity are calculated: a block average epsilon_avg and an average for the whole simulation  epsilon_total\n"
		   << "//For epsilon2, the ensemble average <M>2, which is supposed to converge to zero, is omitted in the calculation\n"
		   << "timestep\tN_particles\tMx\tMy\tMz\t<M>squared\t<M_squared>\tepsilon\tepsilon2\n";

	for (unsigned int i = 0; i < _numOutputs; i++) {

		double squaredM = _outputM[i][0] *  _outputM[i][0] + _outputM[i][1] * _outputM[i][1] + _outputM[i][2] * _outputM[i][2];  // Calculate <M>² from <M> for each block
		blockPermittivity = 1. + 4. * M_PI /(3. * _T * _V) * (_outputSquaredM[i] - squaredM); // Calculates relative permittivity for each block
		blockPermittivity2 = 1. + 4. * M_PI /(3. * _T * _V) * (_outputSquaredM[i]);
		if ( i == 0) {
			averagePermittivity += blockPermittivity * (1. - correctionFirstBlock); // Sums up permittivities for total block average
			averagePermittivity2 += blockPermittivity2 * (1. - correctionFirstBlock);
		}
		else {
			averagePermittivity += blockPermittivity; // Sums up permittivities for total block average
			averagePermittivity2 += blockPermittivity2;
		}
		numSteps = (i + 1 + writeInit) * _writeFrequency + _startingStep;
		writer << numSteps << "\t" << _numParticles[i] << "\t" << _outputM[i][0] << "\t" << _outputM[i][1] << "\t" << _outputM[i][2] << "\t" << squaredM <<"\t" << _outputSquaredM[i] << "\t" << blockPermittivity << "\t" << blockPermittivity2 << "\n";

		if (i == 0) {
			// Sums up block values for <M²> and <M> to calculate total simulation averages
			_totalAverageSquaredM += _outputSquaredM[i] * (1. - correctionFirstBlock);
			for (unsigned int j = 0; j < 3; j++){
				_totalAverageM[j] += _outputM[i][j] * (1. - correctionFirstBlock);
			}
		}
		else {
			_totalAverageSquaredM += _outputSquaredM[i];
			for (unsigned int j = 0; j < 3; j++){
				_totalAverageM[j] += _outputM[i][j];
			}
		}
	}
	// Permittivity as total simulation average (i.e. NOT block average)
	double permittivityTotal = 1. + 4. * M_PI /(3. * _T * _V) * (_totalAverageSquaredM / ((double)_numOutputs - correctionFirstBlock)- (_totalAverageM[0] * _totalAverageM[0] + _totalAverageM[1] * _totalAverageM[1] + _totalAverageM[2] * _totalAverageM[2]) / (((double)_numOutputs - correctionFirstBlock) * ((double)_numOutputs - correctionFirstBlock)));
	double permittivityTotal2 = 1. + 4. * M_PI /(3. * _T * _V) * (_totalAverageSquaredM / ((double)_numOutputs - correctionFirstBlock));
	writer << "\n\n epsilon_total = " << permittivityTotal <<"\tepsilon_total2 = " << permittivityTotal2<<"\t epsilon_avg = "<< averagePermittivity / ((double)_numOutputs - correctionFirstBlock) <<"\tepsilon_avg2 = " << averagePermittivity2 / ((double)_numOutputs - correctionFirstBlock)<<"\tT = " << _T << "\tV = " << _V << "\t<M2> = " << _totalAverageSquaredM / ((double)_numOutputs - correctionFirstBlock) <<"\t<M>2 = " <<  (_totalAverageM[0] * _totalAverageM[0] + _totalAverageM[1] * _totalAverageM[1] + _totalAverageM[2] * _totalAverageM[2]) / (((double)_numOutputs - correctionFirstBlock) * ((double)_numOutputs - correctionFirstBlock)) << std::endl;
	writer.close();
	_ravStream.close();
}

