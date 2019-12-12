
#include "Permittivity.h"
void Permittivity::readXML(XMLfileUnits& xmlconfig) {
	global_log->info() << "Calculation of relative permittivity enabled." << std::endl;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	xmlconfig.getNodeValue("initstatistics", _initStatistics);
	xmlconfig.getNodeValue("recordingtimesteps", _recordingTimesteps);
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
}

void Permittivity::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	_totalNumTimeSteps = _simulation.getNumTimesteps();
	_numOutputs = floor(_totalNumTimeSteps / _writeFrequency);
	_currentOutputNum = 0;
	_accumulatedSteps = 0;
	_numParticlesLocal = 0;
	_totalAverageSquaredM = 0.;


	for (unsigned long i = 0; i<_writeFrequency;i++){
		for (unsigned int j = 0; j<3; j++){
			_localM[i][j] = 0;
			_globalM[i][j] = 0;
		}
	}

	for (unsigned int i = 0; i < _numOutputs; i++) {
		_outputSquaredM[i] = 0.;
		_permittivity[i] = 0.;
		_numParticles[i] = 0;
		for (unsigned int j = 0; j < 3; j++) {
			_outputM[i][j] = 0.;
		}
	}

	for (unsigned int i = 0; i < 3; i++) {
		_totalAverageM[i] = 0.;
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
				global_log->error() << "Wrong dipole vector chosen! Please always choose [eMyx eMyy eMyz] = [0 0 1] when using the permittivity plugin" << endl;
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

void Permittivity::collect(DomainDecompBase* domainDecomp) {

	//Calculate global M for each time step of the current block from local values
	domainDecomp->collCommInit(3 * _accumulatedSteps);
	for (unsigned long j = 0; j < _accumulatedSteps; j++){
		for (unsigned i = 0; i < 3; i++) {
			domainDecomp->collCommAppendDouble(_localM[j][i]);
		}
	}
	domainDecomp->collCommAllreduceSum();

	for (unsigned long j = 0; j < _accumulatedSteps; j++){
		for (unsigned long i = 0; i < 3; i++) {
			_globalM[j][i] = domainDecomp->collCommGetDouble();
		}
	}
		
	domainDecomp->collCommFinalize();
	
	// Calculates ensemble averages <M²> and <M> for current block
	for (unsigned long i = 0; i < _accumulatedSteps; i++) {
		_outputSquaredM[_currentOutputNum] += _globalM[i][0] * _globalM[i][0] + _globalM[i][1] * _globalM[i][1] + _globalM[i][2] * _globalM[i][2];
		for (unsigned int j = 0; j < 3; j++) {
			_outputM[_currentOutputNum][j] += _globalM[i][j];
		}
	}
		
		for (unsigned int i = 0; i < 3; i++) {
			_outputM[_currentOutputNum][i] /= (double)_accumulatedSteps;
		}
		_outputSquaredM[_currentOutputNum] /= (double)_accumulatedSteps;

		domainDecomp->collCommInit(1);
		domainDecomp->collCommAppendUnsLong(_numParticlesLocal);
		domainDecomp->collCommAllreduceSum();
		_numParticles[_currentOutputNum] = domainDecomp->collCommGetUnsLong();
		domainDecomp->collCommFinalize();	
		_currentOutputNum++;
}

void Permittivity::reset() {
	for (unsigned long i = 0; i<_accumulatedSteps;i++){
		for (unsigned j = 0; j<3; j++){
			_localM[i][j] = 0;
			_globalM[i][j] = 0;
		}
	}
	_accumulatedSteps = 0;
	_numParticlesLocal = 0;

}

void Permittivity::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	if (simstep > _initStatistics && simstep % _recordingTimesteps == 0) {
		record(particleContainer);
		_accumulatedSteps++;
	}
	
	if (simstep > _initStatistics && simstep % _writeFrequency == 0 && simstep > 0) {
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
	
	std::array<double, 3> simBoxSize = {domain->getGlobalLength(0), domain->getGlobalLength(1),
										domain->getGlobalLength(2)};
	double V = simBoxSize[0] * simBoxSize[1] * simBoxSize[2]; // Simulation box volume
	double T = domain->getTargetTemperature(0); // Temperature
	double averagePermittivity = 0.; 
	string permName = _outputPrefix + ".perm";
	ofstream writer(permName.c_str());
	writer.precision(7);
	
	writer << "//Output generated by Permittivity plugin\n"
		   << "//Two values for the permittivity are calculated: a block average epsilon_avg and an average for the whole simulation epsilon_total\n"
		   << "timestep\tN_particles\tMx\tMy\tMz\t<M>squared\t<M_squared>\tepsilon\n";
		   
	for (unsigned int i = 0; i < _numOutputs; i++) {

		double squaredM = _outputM[i][0] *  _outputM[i][0] + _outputM[i][1] * _outputM[i][1] + _outputM[i][2] * _outputM[i][2];  // Calculate <M>² from <M> for each block
		_permittivity[i] = 1. + 4. * M_PI /(3. * T * V) * (_outputSquaredM[i] - squaredM); // Calculates relative permittivity for each block 
		averagePermittivity += _permittivity[i]; // Sums up permittivities for total block average 
		writer << (i + 1) * _writeFrequency << "\t" << _numParticles[i] << "\t" << _outputM[i][0] << "\t" << _outputM[i][1] << "\t" << _outputM[i][2] << "\t" << squaredM <<"\t" << _outputSquaredM[i] << "\t" << _permittivity[i] << "\n";
		
		// Sums up block values for <M²> and <M> to calculate total simulation averages
		_totalAverageSquaredM += _outputSquaredM[i];
		for (unsigned int j = 0; j < 3; j++){ 
			_totalAverageM[j] += _outputM[i][j];
		}
	}
	// Permittivity as total simulation average (i.e. NOT block average)
	double permittivityTotal = 1. + 4. * M_PI /(3. * T * V) * (_totalAverageSquaredM / (double)_numOutputs - (_totalAverageM[0] * _totalAverageM[0] + _totalAverageM[1] * _totalAverageM[1] + _totalAverageM[2] * _totalAverageM[2]) / ((double)_numOutputs * (double)_numOutputs));
	
	writer << "\n\n epsilon_total = " << permittivityTotal <<"\t epsilon_avg = "<< averagePermittivity / (double)_numOutputs <<"\tT = " << T << "\tV = " << V << "\t<M2> = " << _totalAverageSquaredM / (double)_numOutputs <<"\t<M>2 = " <<  (_totalAverageM[0] * _totalAverageM[0] + _totalAverageM[1] * _totalAverageM[1] + _totalAverageM[2] * _totalAverageM[2]) / ((double)_numOutputs * (double)_numOutputs) << endl;
	writer.close();
}
	
