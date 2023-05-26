// Created by Joshua Marx 08/2019

#include "ODF.h"
#include "WrapOpenMP.h"

void ODF::readXML(XMLfileUnits& xmlconfig) {
	Log::global_log->debug() << "[ODF] enabled. Dipole orientations must be set to [0 0 1]!" << std::endl;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "[ODF] Write frequency: " << _writeFrequency << std::endl;
	xmlconfig.getNodeValue("initstatistics", _initStatistics);
	Log::global_log->info() << "[ODF] Init Statistics: " << _initStatistics << std::endl;
	xmlconfig.getNodeValue("recordingtimesteps", _recordingTimesteps);
	Log::global_log->info() << "[ODF] Recording Timesteps: " << _recordingTimesteps << std::endl;
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[ODF] Output prefix: " << _outputPrefix << std::endl;
	xmlconfig.getNodeValue("phi1increments", _phi1Increments);
	Log::global_log->info() << "[ODF] Phi1 increments: " << _phi1Increments << std::endl;
	xmlconfig.getNodeValue("phi2increments", _phi2Increments);
	Log::global_log->info() << "[ODF] Phi2 increments: " << _phi2Increments << std::endl;
	xmlconfig.getNodeValue("gammaincrements", _gammaIncrements);
	Log::global_log->info() << "[ODF] Gamma increments: " << _gammaIncrements << std::endl;
	xmlconfig.getNodeValue("shellcutoff", _shellCutOff);
	Log::global_log->info() << "[ODF] Shell cutoff: " << _shellCutOff << std::endl;
}

void ODF::init(ParticleContainer* particleContainer, DomainDecompBase* /*domainDecomp*/, Domain* domain) {
	std::array<double, 3> simBoxSize = {domain->getGlobalLength(0), domain->getGlobalLength(1),
										domain->getGlobalLength(2)};
	_cellProcessor.reset(new ODFCellProcessor(particleContainer->getCutoff(), this, simBoxSize));

	std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
	_numComponents = components->size();
	std::vector<unsigned int> isDipole(_numComponents);
	unsigned int numPairs = 0;
	_readStartingStep = false;

	for (unsigned int i = 0; i < _numComponents; ++i) {
		Component& ci = (*components)[i];
		isDipole[i] = ci.numDipoles();

		if (isDipole[i] == 1) {
			bool orientationIsCorrect = ci.dipole(0).e() == std::array<double, 3>{0,0,1};
			if(orientationIsCorrect == false){
				Log::global_log->error() << "Wrong dipole vector chosen! Please always choose [eMyx eMyy eMyz] = [0 0 1] when using the ODF plugin" << std::endl;
			}
			numPairs++;
		}
	}

	_numPairs = numPairs * numPairs;
	_numElements = _phi1Increments * _phi2Increments * _gammaIncrements + 1;
	Log::global_log->info() << "ODF arrays contains " << _numElements << " elements each for " << _numPairs << "pairings"
					   << std::endl;
	_ODF11.resize(_numElements);
	_ODF12.resize(_numElements);
	_ODF21.resize(_numElements);
	_ODF22.resize(_numElements);

	using vecType2D = decltype(_threadLocalODF11);
	auto resize2D = [](vecType2D& vec, size_t newSizeOuter, size_t newSizeInner) {
		vec.resize(newSizeOuter);
		using vecType = decltype(vec[0]);
		std::for_each(vec.begin(), vec.end(), [newSizeInner](vecType& v) { v.resize(newSizeInner); });
	};

	resize2D(_threadLocalODF11, mardyn_get_max_threads(), _numElements);
	resize2D(_threadLocalODF12, mardyn_get_max_threads(), _numElements);
	resize2D(_threadLocalODF21, mardyn_get_max_threads(), _numElements);
	resize2D(_threadLocalODF22, mardyn_get_max_threads(), _numElements);

	if (_numPairs < 1) {
		Log::global_log->error() << "No components with dipoles. ODF's not being calculated!" << std::endl;
	} else if (_numPairs > 4) {
		Log::global_log->error()
			<< "Number of pairings for ODF calculation too high. Current maximum number of ODF pairings is 4." << std::endl;
	}

	reset();
}

void ODF::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep > _initStatistics && simstep % _recordingTimesteps == 0) {
		particleContainer->traverseCells(*_cellProcessor);
	}
}

void ODF::endStep(ParticleContainer* /*particleContainer*/, DomainDecompBase* domainDecomp, Domain* domain,
				  unsigned long simstep) {
	if (_readStartingStep == false) { // this makes sure, that no output is written at the first time step if the simulation is resumed from a restart file. Otherwise the last output from the previous simulation would be overwritten by a useless file
		_readStartingStep = true;
		return;
	}

	int mpi_rank = domainDecomp->getRank();
	if (simstep > _initStatistics && simstep % _writeFrequency == 0) {
		collect(domainDecomp);
		if (mpi_rank == 0){
			output(domain, simstep);
		}
		reset();
	}
}

void ODF::reset() {
	Log::global_log->info() << "[ODF] resetting data sets" << std::endl;

	//	// C++ 14:
	//	auto fillZero = [](auto& vec) {std::fill(vec.begin(), vec.end(), 0);};
	using vecType = decltype(_ODF11);
	auto fillZero = [](vecType& vec) { std::fill(vec.begin(), vec.end(), 0); };

	fillZero(_ODF11);
	fillZero(_ODF12);
	fillZero(_ODF21);
	fillZero(_ODF22);

	std::for_each(_threadLocalODF11.begin(), _threadLocalODF11.end(), fillZero);
	std::for_each(_threadLocalODF12.begin(), _threadLocalODF12.end(), fillZero);
	std::for_each(_threadLocalODF21.begin(), _threadLocalODF21.end(), fillZero);
	std::for_each(_threadLocalODF22.begin(), _threadLocalODF22.end(), fillZero);
}

void ODF::calculateOrientation(const std::array<double, 3>& simBoxSize, const Molecule& mol1, const Molecule& mol2,
							   const std::array<double, 3>& orientationVector1) {

	// TODO Implement rotation matrices to calculate orientations for dipole direction unit vectors other than [0 0 1];
	
	double Quaternion2[4], orientationVector2[3], distanceVector12[3], moleculeDistance1D[3], auxiliaryVector1[3], auxiliaryVector2[3], projectionVector1[3], projectionVector2[3]; 
	auto cid = mol1.getComponentLookUpID();
	double cosPhi1, cosPhi2, cosGamma12, Gamma12, absoluteProjection1, absoluteProjection2, absoluteDistanceVector12/*, shellCutoff = _shellCutOff[cid]*/;
	double roundingThreshold = 0.0001;
	unsigned long indexPhi1, indexPhi2, indexGamma12, elementID, elementIDreversed;
	unsigned maximumIncrements;
	bool assignPhi1, assignPhi2, assignGamma12;  

	
	double distanceSquared = 0.;
	// Determine modular distance between molecules
	for (int i = 0; i < 3; i++) {
		moleculeDistance1D[i] = mol1.r(i) - mol2.r(i);
		if (abs(moleculeDistance1D[i]) > 0.5 * simBoxSize[i]) {
			moleculeDistance1D[i] = simBoxSize[i] - abs(moleculeDistance1D[i]);
		}
		distanceSquared += moleculeDistance1D[i] * moleculeDistance1D[i];
	}
	
	/*if (_mixingRule == 1) {
		shellCutoff = 1. / 2. * _shellCutOff[cid] + _shellCutOff[mol2.getComponentLookUpID()];
	}*/

	if (distanceSquared < _shellCutOff * _shellCutOff && mol1.getID() != mol2.getID()) {
		absoluteDistanceVector12 = 0.;

		// reading second molecule's quaternions

		Quaternion2[0] = mol2.q().qw();
		Quaternion2[1] = mol2.q().qx();
		Quaternion2[2] = mol2.q().qy();
		Quaternion2[3] = mol2.q().qz();

		cosPhi1 = 0.;
		cosPhi2 = 0.;
		cosGamma12 = 0.;
		absoluteProjection1 = 0.;
		absoluteProjection2 = 0.;
		indexPhi1 = 0;
		indexPhi2 = 0;
		indexGamma12 = 0;
		assignPhi1 = false;
		assignPhi2 = false;
		assignGamma12 = false;

		// calculate distance vector between molecules
		for (unsigned i = 0; i < 3; i++) {
			distanceVector12[i] = mol2.r(i) - mol1.r(i);
			if (distanceVector12[i] > 0.5 * simBoxSize[i]) {
				distanceVector12[i] = -(simBoxSize[i] - distanceVector12[i]);
			} else if (distanceVector12[i] < -0.5 * simBoxSize[i]) {
				distanceVector12[i] = -(distanceVector12[i] - simBoxSize[i]);
			}
			absoluteDistanceVector12 += distanceVector12[i] * distanceVector12[i];
		}

		absoluteDistanceVector12 = sqrt(absoluteDistanceVector12);
		// norm the distance vector
		for (double& i : distanceVector12) {
			i /= absoluteDistanceVector12;
		}

		// calculate dipolar orientation vector from quaternion
		orientationVector2[0] = 2 * (Quaternion2[1] * Quaternion2[3] + Quaternion2[0] * Quaternion2[2]);
		orientationVector2[1] = 2 * (Quaternion2[2] * Quaternion2[3] - Quaternion2[0] * Quaternion2[1]);
		orientationVector2[2] = 1 - 2 * (Quaternion2[1] * Quaternion2[1] + Quaternion2[2] * Quaternion2[2]);

		// calculate projection of the vectors onto plane perpendicular to the distance vector with cross product for calculation of the torque angle gamma
		
		auxiliaryVector1[0] = orientationVector1[1] * distanceVector12[2] - orientationVector1[2] * distanceVector12[1];
		auxiliaryVector1[1] = orientationVector1[2] * distanceVector12[0] - orientationVector1[0] * distanceVector12[2];
		auxiliaryVector1[2] = orientationVector1[0] * distanceVector12[1] - orientationVector1[1] * distanceVector12[0];

		projectionVector1[0] = distanceVector12[1] * auxiliaryVector1[2] - distanceVector12[2] * auxiliaryVector1[1];
		projectionVector1[1] = distanceVector12[2] * auxiliaryVector1[0] - distanceVector12[0] * auxiliaryVector1[2];
		projectionVector1[2] = distanceVector12[0] * auxiliaryVector1[1] - distanceVector12[1] * auxiliaryVector1[0];

		auxiliaryVector2[0] = orientationVector2[1] * distanceVector12[2] - orientationVector2[2] * distanceVector12[1];
		auxiliaryVector2[1] = orientationVector2[2] * distanceVector12[0] - orientationVector2[0] * distanceVector12[2];
		auxiliaryVector2[2] = orientationVector2[0] * distanceVector12[1] - orientationVector2[1] * distanceVector12[0];

		projectionVector2[0] = distanceVector12[1] * auxiliaryVector2[2] - distanceVector12[2] * auxiliaryVector2[1];
		projectionVector2[1] = distanceVector12[2] * auxiliaryVector2[0] - distanceVector12[0] * auxiliaryVector2[2];
		projectionVector2[2] = distanceVector12[0] * auxiliaryVector2[1] - distanceVector12[1] * auxiliaryVector2[0];

		// calculate cos(phi) and norm of projection vector
		for (unsigned i = 0; i < 3; i++) {
			cosPhi1 += distanceVector12[i] * orientationVector1[i];
			cosPhi2 -= distanceVector12[i] * orientationVector2[i];
			absoluteProjection1 += projectionVector1[i] * projectionVector1[i];
			absoluteProjection2 += projectionVector2[i] * projectionVector2[i];
		}

		absoluteProjection1 = sqrt(absoluteProjection1);
		absoluteProjection2 = sqrt(absoluteProjection2);

		// calculate cos(gamma) as dot product of projections
		for (unsigned i = 0; i < 3; i++) {
			projectionVector1[i] /= absoluteProjection1;
			projectionVector2[i] /= absoluteProjection2;
			cosGamma12 += projectionVector1[i] * projectionVector2[i];
		}

		// precaution to prevent numerically intractable values (e.g. VERY close to zero but not zero) and
		// values just barely out of boundaries -1/1, by rounding to 0,-1 or 1 respectively

		if (abs(cosPhi1) < roundingThreshold || abs(abs(cosPhi1) - 1) < roundingThreshold) {
			cosPhi1 = round(cosPhi1);
		}
		if (abs(cosPhi2) < roundingThreshold || abs(abs(cosPhi2) - 1) < roundingThreshold) {
			cosPhi2 = round(cosPhi2);
		}
		if (abs(cosGamma12) < roundingThreshold || abs(abs(cosGamma12) - 1) < roundingThreshold) {
			cosGamma12 = round(cosGamma12);
		}

		Gamma12 = acos(cosGamma12);
		
		// determine array element
		// NOTE: element 0 of array ODF is unused

		maximumIncrements = std::max(_phi1Increments, _phi2Increments);
		maximumIncrements = std::max(maximumIncrements, _gammaIncrements);
		
		// calculate indices for phi1, phi2 and gamma12 for bin assignment
		for (unsigned i = 0; i < maximumIncrements; i++) {
			if (1. - i * 2. / (double)_phi1Increments >= cosPhi1 &&
				cosPhi1 > 1. - (i + 1) * 2. / (double)_phi1Increments) {
				indexPhi1 = i;
				assignPhi1 = true;
			}

			if (1. - i * 2. / (double)_phi2Increments >= cosPhi2 &&
				cosPhi2 > 1. - (i + 1) * 2. / (double)_phi2Increments) {
				indexPhi2 = i;
				assignPhi2 = true;
			}

			if (i * M_PI / (double)_gammaIncrements <= Gamma12 && Gamma12 < (i + 1) * M_PI / (double)_gammaIncrements) {
				indexGamma12 = i + 1;
				assignGamma12 = true;
			}
		}

		if (indexGamma12 == _gammaIncrements + 1) {
			indexGamma12 = _gammaIncrements;
		}

		// manually assign bin for cos(...) == M_PI/-1, because loop only includes values < pi
		if (assignPhi1 == 0 && cosPhi1 == -1.) {
			indexPhi1 = _phi1Increments - 1;
			assignPhi1 = true;
		}

		if (assignPhi2 == 0 && cosPhi2 == -1.) {
			indexPhi2 = _phi2Increments - 1;
			assignPhi2 = true;
		}

		if (assignGamma12 == 0 && Gamma12 == M_PI) {
			indexGamma12 = _gammaIncrements;
			assignGamma12 = true;
		}

		// notification if anything goes wrong during calculataion
		if (assignPhi1 == 0 || assignPhi2 == 0 || assignGamma12 == 0) {
			Log::global_log->warning() << "Array element in ODF calculation not properly assigned!" << std::endl;
			Log::global_log->warning() << "Mol-ID 1 = " << mol1.getID() << "  Mol-ID 2 = " << mol2.getID() << std::endl;
			Log::global_log->warning() << "orientationVector1=" << orientationVector1[0] << " " << orientationVector1[1] << " " << orientationVector1[2] << " " << std::endl;
			Log::global_log->warning() << "orientationVector2=" << orientationVector2[0] << " " << orientationVector2[1] << " " << orientationVector2[2] << " " << std::endl;
			Log::global_log->warning() << "distanceVector12=" << distanceVector12[0] << " " << distanceVector12[1] << " " << distanceVector12[2] << " " << std::endl;
			Log::global_log->warning() << "[cosphi1 cosphi2 cosgamma12] = [" << cosPhi1 << " " << cosPhi2 << " "
								  << cosGamma12 << "]" << std::endl;
			Log::global_log->warning() << "indices are " << indexPhi1 << " " << indexPhi2 << " " << indexGamma12 << std::endl;
		}
		
		// assignment of bin ID
		elementID = indexPhi1 * _phi2Increments * _gammaIncrements + (indexPhi2 * _gammaIncrements) + indexGamma12;
		elementIDreversed = indexPhi2 * _phi2Increments * _gammaIncrements + (indexPhi1 * _gammaIncrements) + indexGamma12; 
		//the ODFcellProcessor calculates every particle interaction only once. Therefore the reverse interaction is considered here as well

		// determine component pairing and add to bin

		if (cid == 0 && mol2.getComponentLookUpID() == 0) {
			_threadLocalODF11[mardyn_get_thread_num()][elementID]++;
			_threadLocalODF11[mardyn_get_thread_num()][elementIDreversed]++;
		}

		else if (cid == 0 && mol2.getComponentLookUpID() == 1) {
			_threadLocalODF12[mardyn_get_thread_num()][elementID]++;
			_threadLocalODF21[mardyn_get_thread_num()][elementIDreversed]++;
		}

		else if (cid == 1 && mol2.getComponentLookUpID() == 1) {
			_threadLocalODF22[mardyn_get_thread_num()][elementID]++;
			_threadLocalODF22[mardyn_get_thread_num()][elementIDreversed]++;
		}

		else {
			_threadLocalODF21[mardyn_get_thread_num()][elementID]++;
			_threadLocalODF12[mardyn_get_thread_num()][elementIDreversed]++;
		}
	}
}

void ODF::collect(DomainDecompBase* domainDecomp) {
	// accumulate thread buffers
	std::vector<unsigned long> localODF11(_threadLocalODF11[0].size(), 0ul);
	std::vector<unsigned long> localODF12(_threadLocalODF12[0].size(), 0ul);
	std::vector<unsigned long> localODF21(_threadLocalODF21[0].size(), 0ul);
	std::vector<unsigned long> localODF22(_threadLocalODF22[0].size(), 0ul);
	for (size_t t = 0; t < static_cast<size_t>(mardyn_get_max_threads()); ++t) {
		using plusType = unsigned long;
		std::transform(localODF11.begin(), localODF11.end(), _threadLocalODF11[t].begin(), localODF11.begin(),
					   std::plus<plusType>());
		std::transform(localODF12.begin(), localODF12.end(), _threadLocalODF12[t].begin(), localODF12.begin(),
					   std::plus<plusType>());
		std::transform(localODF21.begin(), localODF21.end(), _threadLocalODF21[t].begin(), localODF21.begin(),
					   std::plus<plusType>());
		std::transform(localODF22.begin(), localODF22.end(), _threadLocalODF22[t].begin(), localODF22.begin(),
					   std::plus<plusType>());
	}

	if (_numPairs == 1) {
		domainDecomp->collCommInit(_numElements);

		for (unsigned long i = 0; i < _numElements; i++) {
			domainDecomp->collCommAppendUnsLong(localODF11[i]);
		}
		domainDecomp->collCommAllreduceSum();

		for (unsigned long i = 0; i < _numElements; i++) {
			_ODF11[i] = domainDecomp->collCommGetUnsLong();
		}
		domainDecomp->collCommFinalize();
	}

	else {
		domainDecomp->collCommInit(_numElements * 4);

		for (unsigned long i = 0; i < _numElements; i++) {
			domainDecomp->collCommAppendUnsLong(localODF11[i]);
			domainDecomp->collCommAppendUnsLong(localODF12[i]);
			domainDecomp->collCommAppendUnsLong(localODF22[i]);
			domainDecomp->collCommAppendUnsLong(localODF21[i]);
		}
		domainDecomp->collCommAllreduceSum();

		for (unsigned long i = 0; i < _numElements; i++) {
			_ODF11[i] = domainDecomp->collCommGetUnsLong();
			_ODF12[i] = domainDecomp->collCommGetUnsLong();
			_ODF22[i] = domainDecomp->collCommGetUnsLong();
			_ODF21[i] = domainDecomp->collCommGetUnsLong();
		}
		domainDecomp->collCommFinalize();
	}
}

void ODF::output(Domain* /*domain*/, long unsigned timestep) {
	Log::global_log->info() << "[ODF] writing output" << std::endl;
	// Setup outfile
	constexpr double piHalf = 0.5 * M_PI;
	double cosPhi1 = 1. + 1. / (double)_phi1Increments;
	double cosPhi2 = 1. - 2. / _phi2Increments;
	double Gamma12 = 0.;
	double averageODF11 = 0.;
	std::string prefix;
	std::ostringstream osstrm;
	osstrm << _outputPrefix;
	osstrm.fill('0');
	osstrm.width(7);
	osstrm << std::right << timestep;
	prefix = osstrm.str();
	osstrm.str("");
	osstrm.clear();

	if (_numPairs == 1) {
		std::string ODF11name = prefix + ".ODF11";
		std::ofstream outfile(ODF11name.c_str());
		outfile.precision(6);
		
		for (unsigned long i = 1; i < _numElements; i++) {
			averageODF11 += (double)_ODF11[i];
		}
		
		averageODF11 /= ((double)_numElements - 1.);
		
		outfile << "//Output generated by ODF plugin\n"
				<< "//The output is to be interpreted as follows: an orientation with a value of ODF_normed = 1 has the same probability as random orientation. ODF_normed = 1.4 means 40% higher probability than random orientation, 0.4 means 60% lower probability than random orientation etc.\n"
				<< "//Angular distribution at time step = " << timestep << " for component pairing pairing 11\n";
		outfile << "cosPhi1\tcosPhi2\tGamma12\tODF_normed\n";

		for (unsigned long i = 0; i < _numElements - 1; i++) {
			Gamma12 += M_PI / (double)_gammaIncrements;
			if (i % _gammaIncrements == 0) {
				cosPhi2 -= 2. / (double)_phi2Increments;
				Gamma12 = piHalf / (double)_gammaIncrements;
			}
			if (i % (_gammaIncrements * _phi2Increments) == 0) {
				cosPhi1 -= 2. / (double)_phi1Increments;
				cosPhi2 = 1. - 1. / (double)_phi2Increments;
			}
			outfile << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << (double)_ODF11[i + 1] / averageODF11 << "\n";
		}
		outfile.close();
	} else {
		double averageODF12 = 0.;
		double averageODF21 = 0.;
		double averageODF22 = 0.;
		
		for (unsigned long i = 1; i < _numElements; i++) {
			averageODF11 += (double)_ODF11[i];
			averageODF12 += (double)_ODF12[i];
			averageODF21 += (double)_ODF21[i];
			averageODF22 += (double)_ODF22[i];
		}
		
		averageODF11 /= ((double)_numElements - 1.);
		averageODF12 /= ((double)_numElements - 1.);
		averageODF21 /= ((double)_numElements - 1.);
		averageODF22 /= ((double)_numElements - 1.);
		
		std::string ODF11name = prefix + ".ODF11";
		std::string ODF12name = prefix + ".ODF12";
		std::string ODF22name = prefix + ".ODF22";
		std::string ODF21name = prefix + ".ODF21";

		std::ofstream ODF11(ODF11name.c_str());
		std::ofstream ODF12(ODF12name.c_str());
		std::ofstream ODF22(ODF22name.c_str());
		std::ofstream ODF21(ODF21name.c_str());
		ODF11.precision(5);
		ODF12.precision(5);
		ODF22.precision(5);
		ODF21.precision(5);
		
		

		ODF11 << "//Output generated by ODF plugin\n"
			  << "//The output is to be interpreted as follows: an orientation with a value of ODF_normed = 1 has the same probability as random orientation. ODF_normed = 1.4 means 40% higher probability than random orientation, 0.4 means 60% lower probability than random orientation etc.\n"
			  << "//Angular distribution at time step = " << timestep << " for component pairing pairing 11.\n";
		ODF11 << "cosPhi1\tcosPhi2\tcosGamma12\tODF_normed\n";
		ODF12 << "//Output generated by ODF plugin\n"
			  << "//The output is to be interpreted as follows: an orientation with a value of ODF_normed = 1 has the same probability as random orientation. ODF_normed = 1.4 means 40% higher probability than random orientation, 0.4 means 60% lower probability than random orientation etc.\n"
		      << "//Angular distribution at time step = " << timestep << " for component pairing pairing 12\n";
		ODF12 << "cosPhi1\tcosPhi2\tcosGamma12\tODF_normed\n";
		ODF22 << "//Output generated by ODF plugin\n"
			  << "//The output is to be interpreted as follows: an orientation with a value of ODF_normed = 1 has the same probability as random orientation. ODF_normed = 1.4 means 40% higher probability than random orientation, 0.4 means 60% lower probability than random orientation etc.\n"
		      << "//Angular distribution at time step = " << timestep << " for component pairing pairing 22\n";
		ODF22 << "phi1\tphi2\tgamm12\tODF_normed\n";
		ODF21 << "//Output generated by ODF plugin\n"
			  << "//The output is to be interpreted as follows: an orientation with a value of ODF_normed = 1 has the same probability as random orientation. ODF_normed = 1.4 means 40% higher probability than random orientation, 0.4 means 60% lower probability than random orientation etc.\n"
			  << "//Angular distribution at time step = " << timestep << " for component pairing pairing 21\n";
		ODF21 << "cosPhi1\tcosPhi2\tcosGamma12\tODF_normed\n";
			
				

		for (unsigned long i = 0; i < _numElements - 1; i++) {
			Gamma12 += M_PI / (double)_gammaIncrements;
			if (i % _gammaIncrements == 0) {
				cosPhi2 -= 2. / (double)_phi2Increments;
				Gamma12 = piHalf / (double)_gammaIncrements;
			}
			if (i % (_gammaIncrements * _phi2Increments) == 0) {
				cosPhi1 -= 2. / (double)_phi1Increments;
				cosPhi2 = 1. - 1. / (double)_phi2Increments;
			}

			ODF11 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << (double)_ODF11[i + 1] / averageODF11 << "\n";
			ODF12 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << (double)_ODF12[i + 1] / averageODF12 << "\n";
			ODF22 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << (double)_ODF22[i + 1] / averageODF22 << "\n";
			ODF21 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << (double)_ODF21[i + 1] / averageODF21 << "\n";
		}
		ODF11.close();
		ODF12.close();
		ODF22.close();
		ODF21.close();
	}
}
