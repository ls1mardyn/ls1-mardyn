
#include "Domain.h"
#include "longRange/Spherical.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"

#include <vector>
#include <cmath>
#include <algorithm>

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include "utils/FileUtils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;  // TODO: remove (Fabio says this is bad practice!)
using Log::global_log;

Spherical::Spherical(double /*cutoffT*/, double cutoffLJ, Domain* domain, DomainDecompBase* domainDecomposition,
					 ParticleContainer* particleContainer, Simulation* simulation) {
	global_log->info() << "Long Range Correction for spherical interfaces is used" << std::endl;
	// default values
	rc = cutoffLJ;
	_domain = domain;
	_domainDecomposition = domainDecomposition;
	_particleContainer = particleContainer;
	NShells = 70;
	NSMean = 50;
	nSamples = 0;
	numComp = 1;
	globalNumMols = 31000;
	global_log->info() << "[Long Range Correction] golbalNumMols set to: " << globalNumMols << " during initialization"
					   << std::endl;

	rcmax = 8.0;
	UpotKorrLJ = 0.0;
	VirialKorrLJ = 0.0;
	droplet = true;
	disableLRC = false;
	_outputPrefix = "mardyn";
	_T = 0.0;

	calcFreq = 100;
	writeFreq = 10000;
	// making sure that writeFreq is a multiple of calcFreq:
	writeFreq = (writeFreq / calcFreq) * calcFreq;
	if (writeFreq == 0) {
		writeFreq = calcFreq;
	}
}

Spherical::~Spherical() {}

void Spherical::init() {
	global_log->info() << "[Long Range Correction] Initializing. Is this function called twice?!" << std::endl;

	std::cout << "[Long Range Correction] floatFix rho" << std::endl;

	globalNumMols = _domain->getglobalNumMolecules(true, _particleContainer, _domainDecomposition);
	global_log->info() << "[Long Range Correction] global number of molecules: " << globalNumMols << std::endl;

	vector<Component>& components = *_simulation.getEnsemble()->getComponents();
	numComp = components.size();
	global_log->info() << "[Long Range Correction] Anzahl Comps: " << numComp << std::endl;
	resizeExactly(numLJ, numComp);
	numLJSum = 0;
	resizeExactly(numLJSum2, numComp);
	for (unsigned i = 0; i < numComp; i++) {
		numLJSum2[i] = 0;
	}
	for (unsigned i = 0; i < numComp; i++) {
		numLJ[i] = components[i].numLJcenters();
		numLJSum += numLJ[i];
		for (unsigned j = i + 1; j < numComp; j++) {
			numLJSum2[j] += numLJ[i];
		}
	}

	resizeExactly(eLong, numLJSum);

	resizeExactly(ksi, globalNumMols);
	resizeExactly(FcorrX, globalNumMols);         // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrY, globalNumMols);         // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrZ, globalNumMols);         // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrX_global, globalNumMols);  // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrY_global, globalNumMols);  // requires refreshIDs=1 in simulation options
	resizeExactly(FcorrZ_global, globalNumMols);  // requires refreshIDs=1 in simulation options
	resizeExactly(lowerS, globalNumMols);         // requires refreshIDs=1 in simulation options
	resizeExactly(interS, globalNumMols);         // requires refreshIDs=1 in simulation options
	resizeExactly(upperS, globalNumMols);         // requires refreshIDs=1 in simulation options

	std::fill(ksi.begin(), ksi.end(), 0.0);
	std::fill(FcorrX.begin(), FcorrX.end(), 0.0);
	std::fill(FcorrY.begin(), FcorrY.end(), 0.0);
	std::fill(FcorrZ.begin(), FcorrZ.end(), 0.0);
	std::fill(FcorrX_global.begin(), FcorrX_global.end(), 0.0);
	std::fill(FcorrY_global.begin(), FcorrY_global.end(), 0.0);
	std::fill(FcorrZ_global.begin(), FcorrZ_global.end(), 0.0);
	std::fill(lowerS.begin(), lowerS.end(), 0.0);
	std::fill(interS.begin(), interS.end(), 0.0);
	std::fill(upperS.begin(), upperS.end(), 0.0);

	resizeExactly(RShells, NShells);
	resizeExactly(RShells2, NShells);
	resizeExactly(RShells3, NShells);
	resizeExactly(VShells, NShells);
	resizeExactly(rhoShellsTemp, NShells);
	resizeExactly(rhoShellsTemp_global, NShells);
	resizeExactly(rhoShellsMean, NShells * NSMean);
	resizeExactly(rhoShellsAvg, NShells);
	resizeExactly(rhoShellsAvg_global, NShells);
	resizeExactly(TShellsAvg_global, NShells);
	resizeExactly(TShellsTemp, NShells);
	resizeExactly(TShellsAvg, NShells);
	resizeExactly(rhoShells, NShells);
	resizeExactly(rhoShells_global, NShells);
	resizeExactly(rhoShellsT, NShells);
	resizeExactly(PartShells, globalNumMols);
	resizeExactly(UShells_Mean, NShells);
	resizeExactly(UShells_Mean_global, NShells);
	resizeExactly(FShells_Mean, NShells);
	resizeExactly(FShells_Mean_global, NShells);
	resizeExactly(PNShells_Mean, NShells);
	resizeExactly(PNShells_Mean_global, NShells);
	resizeExactly(PTShells_Mean, NShells);
	resizeExactly(PTShells_Mean_global, NShells);
	resizeExactly(PNShellsCorr_Avg, NShells);
	resizeExactly(PNShellsCorr_Avg_global, NShells);
	resizeExactly(PTShellsCorr_Avg, NShells);
	resizeExactly(PTShellsCorr_Avg_global, NShells);
	resizeExactly(VirShells_N, NShells);
	resizeExactly(VirShells_N_global, NShells);
	resizeExactly(VirShells_T, NShells);
	resizeExactly(VirShells_T_global, NShells);

	std::fill(RShells.begin(), RShells.end(), 0.0);
	std::fill(RShells2.begin(), RShells2.end(), 0.0);
	std::fill(RShells3.begin(), RShells3.end(), 0.0);
	std::fill(VShells.begin(), VShells.end(), 0.0);
	std::fill(rhoShellsTemp.begin(), rhoShellsTemp.end(), 0.0);
	std::fill(rhoShellsTemp_global.begin(), rhoShellsTemp_global.end(), 0.0);
	std::fill(rhoShellsMean.begin(), rhoShellsMean.end(), 0.0);
	std::fill(rhoShellsAvg.begin(), rhoShellsAvg.end(), 0.0);
	std::fill(rhoShellsAvg_global.begin(), rhoShellsAvg_global.end(), 0.0);
	std::fill(TShellsAvg_global.begin(), TShellsAvg_global.end(), 0.0);
	std::fill(TShellsTemp.begin(), TShellsTemp.end(), 0.0);
	std::fill(TShellsAvg.begin(), TShellsAvg.end(), 0.0);
	std::fill(rhoShells.begin(), rhoShells.end(), 0.0);
	std::fill(rhoShells_global.begin(), rhoShells_global.end(), 0.0);
	std::fill(rhoShellsT.begin(), rhoShellsT.end(), 0.0);
	std::fill(PartShells.begin(), PartShells.end(), 0.0);
	std::fill(UShells_Mean.begin(), UShells_Mean.end(), 0.0);
	std::fill(UShells_Mean_global.begin(), UShells_Mean_global.end(), 0.0);
	std::fill(FShells_Mean.begin(), FShells_Mean.end(), 0.0);
	std::fill(FShells_Mean_global.begin(), FShells_Mean_global.end(), 0.0);
	std::fill(PNShells_Mean.begin(), PNShells_Mean.end(), 0.0);
	std::fill(PNShells_Mean_global.begin(), PNShells_Mean_global.end(), 0.0);
	std::fill(PTShells_Mean.begin(), PTShells_Mean.end(), 0.0);
	std::fill(PTShells_Mean_global.begin(), PTShells_Mean_global.end(), 0.0);
	std::fill(PNShellsCorr_Avg.begin(), PNShellsCorr_Avg.end(), 0.0);
	std::fill(PNShellsCorr_Avg_global.begin(), PNShellsCorr_Avg_global.end(), 0.0);
	std::fill(PTShellsCorr_Avg.begin(), PTShellsCorr_Avg.end(), 0.0);
	std::fill(PTShellsCorr_Avg_global.begin(), PTShellsCorr_Avg_global.end(), 0.0);
	std::fill(VirShells_N.begin(), VirShells_N.end(), 0.0);
	std::fill(VirShells_N_global.begin(), VirShells_N_global.end(), 0.0);
	std::fill(VirShells_T.begin(), VirShells_T.end(), 0.0);
	std::fill(VirShells_T_global.begin(), VirShells_T_global.end(), 0.0);

	// Determination of the elongation of the Lennard-Jones sites
	unsigned int counter = 0;
	for (unsigned i = 0; i < numComp; i++) {
		for (unsigned j = 0; j < components[i].numLJcenters(); j++) {
			const LJcenter& ljcenteri = static_cast<const LJcenter&>(components[i].ljcenter(j));
			double dX[3];
			dX[0] = ljcenteri.rx();
			dX[1] = ljcenteri.ry();
			dX[2] = ljcenteri.rz();
			for (unsigned d = 0; d < 3; d++) {
				dX[d] *= dX[d];
			}
			eLong[counter] = sqrt(dX[0] + dX[1] + dX[2]);
			counter++;
		}
	}

	boxlength[0] = _domain->getGlobalLength(0);
	boxlength[1] = _domain->getGlobalLength(1);
	boxlength[2] = _domain->getGlobalLength(2);
	for (unsigned short d = 0; d < 3; ++d) {
		systemcenter[d] = 0.5 * boxlength[d];
	}

	if ((boxlength[0] != boxlength[1]) or (boxlength[0] != boxlength[2])) {
		global_log->warning() << "[Long Range Correction] Domain size should be the same in all directions!" << std::endl;
	}

	_drShells = 0.5 * std::min({boxlength[0], boxlength[1], boxlength[2]}) / (NShells + 1);

	_deltaShells = rc / _drShells;

	double drShells05 = 0.5 * _drShells;
	double maxShells = NShells + _deltaShells;

	for (unsigned int i = 0; i < NShells - 1; i++) {
		RShells[i] = (i + 1) * _drShells;
		RShells2[i] = RShells[i] * RShells[i];
		RShells3[i] = RShells2[i] * RShells[i];
		VShells[i] = (4. / 3.) * M_PI * (pow((RShells[i] + drShells05), 3) - pow((RShells[i] - drShells05), 3));
	}
	RShells[NShells - 1] = NShells * _drShells;
	RShells2[NShells - 1] = RShells[NShells - 1] * RShells[NShells - 1];
	RShells3[NShells - 1] = RShells2[NShells - 1] * RShells[NShells - 1];
	VShells[NShells - 1] =
		(boxlength[0] * boxlength[1] * boxlength[2]) - (4. / 3.) * M_PI * pow((RShells[NShells - 1] - drShells05), 3);


	// NOTE: procedure to calculate the instantaneous inside and outside densities is risky here, because it makes a
	// priori asumptions about where to find the phase boundary! Better Ideas are welcome.
	bulkBoundaries.inside_from = NShells / 10; //10; 
	bulkBoundaries.inside_to = NShells / 5 + 5; // 35;
	bulkBoundaries.outside_from = NShells - 2 - NShells / 10 - 5;  // 72;
	bulkBoundaries.outside_to = NShells - 1;

	// systemcenter is at half of domain size (densities etc are averaged spherically around systemcenter)
	double inside_from_coord = (bulkBoundaries.inside_from*systemcenter[0])/NShells;
	double inside_to_coord = (bulkBoundaries.inside_to*systemcenter[0])/NShells;
	double outside_from_coord = (bulkBoundaries.outside_from*systemcenter[0])/NShells;
	double outside_to_coord = (bulkBoundaries.outside_to*systemcenter[0])/NShells;

	global_log->warning() << "[Long Range Correction] Calculating rho_inside from " << inside_from_coord << " (Shell " << bulkBoundaries.inside_from << ")"
																	   << " to " << inside_to_coord << " (Shell " << bulkBoundaries.inside_to << ")"
														 << " and rho_out from " << outside_from_coord << " (Shell " << bulkBoundaries.outside_from << ")"
																	   << " to " << outside_to_coord << " (Shell " << bulkBoundaries.outside_to << ")" << std::endl;


	// Set names for output files and write header
	// if (_domainDecomposition->getRank() == 0) {
		// filenameTanhParams = _outputPrefix + "_LRCspherical_tanh.dat";
		// ofstream outfilestreamTanhParams(filenameTanhParams, ios::out);
		// outfilestreamTanhParams << std::setw(24) << "simstep";
		// outfilestreamTanhParams << std::setw(24) << "rhov";
		// outfilestreamTanhParams << std::setw(24) << "rhol";
		// outfilestreamTanhParams << std::setw(24) << "D0";
		// outfilestreamTanhParams << std::setw(24) << "R0";
		// outfilestreamTanhParams << std::endl;
		// outfilestreamTanhParams.close();

		// filenameThermData = _outputPrefix + "_ThermData.csv";
		// ofstream outfilestreamThermData(filenameThermData, ios::out);
		// outfilestreamThermData << std::setw(24) << "simstep";
		// // outfilestreamThermData << std::setw(24) << "gamma (iterative)";
		// outfilestreamThermData << std::setw(24) << "gamma[n-2]";
		// outfilestreamThermData << std::setw(24) << "gamma[n-1]";
		// outfilestreamThermData << std::setw(24) << "gamma";
		// outfilestreamThermData << std::setw(24) << "R_gamma";
		// outfilestreamThermData << std::setw(24) << "R_e";
		// outfilestreamThermData << std::setw(24) << "delta";
		// outfilestreamThermData << std::setw(24) << "rhoInside";
		// outfilestreamThermData << std::setw(24) << "rhoOutside ";
		// outfilestreamThermData << std::setw(24) << "pInside";
		// outfilestreamThermData << std::setw(24) << "pOutside";
		// outfilestreamThermData << std::setw(24) << "inside_from"
		// 					   << "";
		// outfilestreamThermData << std::setw(24) << "inside_to"
		// 					   << "";
		// outfilestreamThermData << std::setw(24) << "outside_from"
		// 					   << "";
		// outfilestreamThermData << std::setw(24) << "outside_to"
		// 					   << "";
		// outfilestreamThermData << std::setw(24) << "nSamples"
		// 					   << "";
		// outfilestreamThermData << std::endl;
		// outfilestreamThermData.close();

		// filenameGlobalCorrs = _outputPrefix + "_globalCorrections_0.csv";
	// }
}

void Spherical::readXML(XMLfileUnits& xmlconfig) {
	global_log->info() << "[Long Range Correction] reading XML for paramters of SphericalLRC." << std::endl;
	xmlconfig.getNodeValue("shells", NShells);
	xmlconfig.getNodeValue("droplet", droplet);
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	xmlconfig.getNodeValue("temperature", _T);
	xmlconfig.getNodeValue("disableLRC", disableLRC);
	xmlconfig.getNodeValue("calculationFreq", calcFreq);
	xmlconfig.getNodeValue("writeFreq", writeFreq);

	// making sure that writeFreq is a multiple of calcFreq:
	writeFreq = (writeFreq / calcFreq) * calcFreq;
	if (writeFreq == 0) {
		writeFreq = calcFreq;
	}

	global_log->info() << "[Long Range Correction] Using " << NShells << " shells for profiles to calculate LRC."
					   << std::endl;
	if (droplet) {
		global_log->info() << "[Long Range Correction] System contains a droplet. COMaligner plugin is required."
						   << std::endl;
	} else {
		global_log->info() << "[Long Range Correction] System contains a bubble. COMalignerBubble plugin is required."
						   << std::endl;
	}

	// init
	this->init();
}

void Spherical::calculateLongRange() {
	// global_log->info() << "[Long Range Correction] calculateLongRange has been called" << std::endl;

	int rank = _domainDecomposition->getRank();
	uint64_t simstep = _simulation.getSimulationStep();
	// global_log->info() << "[Long Range Correction] simstep = " << simstep << std::endl;

	std::fill(rhoShellsTemp.begin(), rhoShellsTemp.end(), 0.0);
	std::fill(TShellsTemp.begin(), TShellsTemp.end(), 0.0);

	for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid();
		 ++tempMol) {
		double v2i = tempMol->v(0) * tempMol->v(0) + tempMol->v(1) * tempMol->v(1) +
					 tempMol->v(2) * tempMol->v(2);  // velocity^2 of tempMol

		unsigned long molID = tempMol->getID();
		for (unsigned short d = 0; d < 3; ++d) {
			FcorrX[molID] = systemcenter[0] - tempMol->r(0);
			FcorrY[molID] = systemcenter[1] - tempMol->r(1);
			FcorrZ[molID] = systemcenter[2] - tempMol->r(2);
		}
		ksi[molID] = std::sqrt((tempMol->r(0) - systemcenter[0]) * (tempMol->r(0) - systemcenter[0]) +
							   (tempMol->r(1) - systemcenter[1]) * (tempMol->r(1) - systemcenter[1]) +
							   (tempMol->r(2) - systemcenter[2]) * (tempMol->r(2) - systemcenter[2]));
		double realk = ksi[molID] / _drShells;
		lowerS[molID] = std::min(static_cast<double>(std::floor((realk - _deltaShells))), static_cast<double>(NShells));
		interS[molID] = std::max(static_cast<double>(std::ceil(abs(realk - _deltaShells))), 1.0);
		upperS[molID] =
			std::min(static_cast<double>(std::ceil(realk + _deltaShells)), static_cast<double>(NShells + 1));
		unsigned long k = std::round(realk);
		if (k > NShells - 1) {
			rhoShellsTemp[NShells - 1] += 1.;
			TShellsTemp[NShells - 1] += v2i;
			PartShells[molID] = NShells - 1;
		} else if (k == 0) {
			// rhoShellsTemp[k] += 1.; besser?
			PartShells[molID] = k;
		} else {
			rhoShellsTemp[k - 1] += 1.;
			TShellsTemp[k - 1] += v2i;
			PartShells[molID] = k - 1;
		}
	}

	for (unsigned int i = 0; i < NShells; i++) {
		// TShellsTemp[i] /= 3.*rhoShellsTemp[i] -3; //up to this point, rhoShellsTemp[i] is just # of Particles in Cell
		// i
		TShellsAvg[i] = (TShellsAvg[i] * simstep + TShellsTemp[i]) / (simstep + 1);  // not used

		rhoShellsTemp[i] = rhoShellsTemp[i] / VShells[i];
		rhoShellsAvg[i] += rhoShellsTemp[i];
	}
	// Distribution of the Density Profile and local Temperature to every node
	_domainDecomposition->collCommInit(2 * NShells);
	for (unsigned i = 0; i < NShells; i++) {
		_domainDecomposition->collCommAppendDouble(rhoShellsAvg[i]);
		_domainDecomposition->collCommAppendDouble(TShellsAvg[i]);
	}
	_domainDecomposition->collCommAllreduceSum();
	for (unsigned i = 0; i < NShells; i++) {
		rhoShellsAvg_global[i] = _domainDecomposition->collCommGetDouble();
		TShellsAvg_global[i] =
			_domainDecomposition->collCommGetDouble() / (3 * rhoShellsAvg_global[i] / (simstep + 1) * VShells[i]);
		
		rhoShellsAvg_global[i] /= (simstep + 1);
	}
	_domainDecomposition->collCommFinalize();


	unsigned long MeanIndex = (static_cast<int>(std::floor((simstep - 1) / calcFreq))) % NSMean;  // 1000
	if (((simstep - 1) % calcFreq) == 0) {                                                      // 1000
		std::fill(rhoShellsMean.begin() + NShells * MeanIndex, rhoShellsMean.begin() + NShells * (MeanIndex + 1), 0.0);
	}
	// To get some results before NSMean*calcFreq, initialize rhoShellsMean with density at simstep = 0
	if (simstep == 0) {
		for (unsigned int j = 0; j < NSMean; j++) {
			for (unsigned int i = 0; i < NShells; i++) {
				rhoShellsMean[j * NShells + i] += calcFreq * rhoShellsTemp[i];  // 1000
				if (i==2932) { std::cout << "Init MeanIndex " << MeanIndex << " " << rhoShellsTemp[i] << std::endl; }
			}
		}
	} else {
		for (unsigned int i = 0; i < NShells; i++) {
			if (i==2932) { std::cout << "MeanIndex " << simstep << " " << MeanIndex << " " << rhoShellsMean[NShells * MeanIndex + i] << std::endl; }
			rhoShellsMean[NShells * MeanIndex + i] += rhoShellsTemp[i];
			if (i==2932) { std::cout << "MeanIndex " << simstep << " " << MeanIndex << " " << rhoShellsMean[NShells * MeanIndex + i] << " added " << rhoShellsTemp[i] << std::endl; }
		}
	}

	// mittlere Dichte
	if ((simstep) % calcFreq == 0) {  // 1000
		std::fill(rhoShells.begin(), rhoShells.end(), 0.0);
		std::fill(rhoShells_global.begin(), rhoShells_global.end(), 0.0);
		std::fill(rhoShellsTemp_global.begin(), rhoShellsTemp_global.end(), 0.0);

		for (unsigned int j = 0; j < NShells; j++) {
			for (unsigned int i = 0; i < NSMean; i++) {
				rhoShells[j] += (1./calcFreq) / NSMean * rhoShellsMean[i * NShells + j];  // 1000 --- 0.001
				if (j==2932) { std::cout << "rhoShellsMean " << rank << " " << simstep << " " << i << " " << (1./calcFreq) / NSMean * rhoShellsMean[i * NShells + j] << std::endl; }
			}
			if (j==2932) { std::cout << "rhoShells " << rank << " " << simstep << " " << j << " " << rhoShells[j] << std::endl; }
		}

		// Distribution of the Density Profile to every node
		_domainDecomposition->collCommInit(NShells);
		for (unsigned i = 0; i < NShells; i++) {
			_domainDecomposition->collCommAppendDouble(rhoShells[i]);
		}
		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i = 0; i < NShells; i++) {
			rhoShells_global[i] = _domainDecomposition->collCommGetDouble();
			if (i==2932) { std::cout << "rhoShells_global[2932] " << rank << " " << simstep << " " << i << " " << rhoShells_global[i] << " " << rhoShellsAvg_global[i] << std::endl; }
		}
		_domainDecomposition->collCommFinalize();

		// Distribution of the Density Profile to every node
		_domainDecomposition->collCommInit(NShells);
		for (unsigned i = 0; i < NShells; i++) {
			_domainDecomposition->collCommAppendDouble(rhoShellsTemp[i]);
		}
		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i = 0; i < NShells; i++) {
			rhoShellsTemp_global[i] = _domainDecomposition->collCommGetDouble();
		}
		_domainDecomposition->collCommFinalize();

		// 	// rhol and rhov  ----- this calculation scheme does not work! the instantaneous profile is far to jagged
		// for the "abs(rhol-rhoShells_global[i]) < 0.1*rhol" condition. often return 0 density for rhol in bubbles with
		// large NShells 	double rhol = 0.0;  // actually "rho_inside" 	for (unsigned int i=4; i<14; i++) { 		rhol +=
		// 0.1*rhoShells_global[i];
		// //		std::cout << "rhol" << i << " " << rhol << std::endl;
		// 	}

		// 	for (unsigned int i=14; i<NShells; i++) {
		// 		if ( abs(rhol-rhoShells_global[i]) < 0.1*rhol ) {
		// 			rhol = ((i-1)*rhol + rhoShells_global[i])/i;
		// //			std::cout << "rhol" << i << " " << rhol << std::endl;

		// 		}
		// 	}

		// 	double rhov = 0.0;   // actually "rho_outside"
		// 	for (unsigned int i=NShells-11; i<(NShells-1); i++) {
		// 		rhov += 0.1*rhoShells_global[i];
		// 	}

		// 	for (unsigned int i=1; i<(NShells-20); i++) {
		// 		unsigned int j = NShells-11-i;
		// 		if ( abs(rhov-rhoShells_global[j]) < 0.2*rhov ) {
		// 			rhov = ((11+i-1)*rhov + rhoShells_global[j])/(11+i);
		// 		}
		// 	}

		// NOTE: procedure to calculate the instantaneous inside and outside densities is risky here, because it makes a
		// priori asumptions about where to find the phase boundary! Better Ideas are welcome.
		double rho_in = 0.;
		for (unsigned int i = bulkBoundaries.inside_from; i < bulkBoundaries.inside_to; i++) {
			rho_in += rhoShells_global[i];
		}
		rho_in /= (bulkBoundaries.inside_to - bulkBoundaries.inside_from);

		double rho_out = 0.;
		for (unsigned int i = bulkBoundaries.outside_from; i < bulkBoundaries.outside_to; i++) {
			rho_out += rhoShells_global[i];
		}
		rho_out /= (bulkBoundaries.outside_to - bulkBoundaries.outside_from);

		global_log->info() << "[Long Range Correction] Averaged rho: rho_in = " << rho_in << " rho_out = " << rho_out << std::endl;

		// D0 with 1090 (Baidakov et al.)
		double Dmin = 0.0;
		double Dmax = 0.0;

		if (droplet) {
			double r10 = rho_out + 0.1 * (rho_in - rho_out);
			double r90 = rho_out + 0.9 * (rho_in - rho_out);
			for (unsigned int i = 1; i < (NShells - 10); i++) {  // TODO:WHY IS THERE A HARDCODED 10 HERE?!
				if (rhoShells_global[i] > r90) {
					Dmin = RShells[i];
				}
			}
			for (unsigned int i = 1; i < (NShells - 10); i++) {  // some value/limitation of the search space could be needed for the largest shells in case of droplet (low density)
				// This could be written much better with a "backwards" loop and/or a break condition
				unsigned int j = NShells - i;
				if (rhoShells_global[j] < r10) {
					Dmax = RShells[j];
				}
			}
		} else {
			double r10 = rho_in + 0.1 * (rho_out - rho_in);
			double r90 = rho_in + 0.9 * (rho_out - rho_in);
			for (unsigned int i = 1; i < (NShells - 10); i++) {  // some value/limitation of the search space could be needed for the smallest shells in case of bubble (low density)
				if (rhoShells_global[i] < r10) {
					Dmin = RShells[i];
				}
			}
			for (unsigned int i = 1; i < (NShells - 10); i++) {
				unsigned int j = NShells - i;
				if (rhoShells_global[j] > r90) {
					Dmax = RShells[j];
				}
			}
		}

		double D0 = (Dmax - Dmin);
		double R0 = Dmin + 0.5 * D0;

		// density profile
		for (unsigned int i = 0; i < NShells; i++) {
			rhoShellsT[i] = RhoP(RShells[i], rho_out, rho_in, D0, R0);
		}

		// if (simstep % writeFreq == 0) {
		// 	if (rank == 0) {
		// 		ofstream outfilestreamTanhParams(filenameTanhParams, ios::app);
		// 		outfilestreamTanhParams << std::setw(24) << std::setprecision(std::numeric_limits<double>::digits10)
		// 								<< simstep;
		// 		outfilestreamTanhParams << std::setw(24) << std::scientific
		// 								<< std::setprecision(std::numeric_limits<double>::digits10) << rho_out;
		// 		outfilestreamTanhParams << std::setw(24) << std::scientific
		// 								<< std::setprecision(std::numeric_limits<double>::digits10) << rho_in;
		// 		outfilestreamTanhParams << std::setw(24) << std::scientific
		// 								<< std::setprecision(std::numeric_limits<double>::digits10) << D0;
		// 		outfilestreamTanhParams << std::setw(24) << std::scientific
		// 								<< std::setprecision(std::numeric_limits<double>::digits10) << R0;
		// 		outfilestreamTanhParams << std::endl;
		// 		outfilestreamTanhParams.close();
		// 	}
		// }

		// shift for correction, since outer bulk only needs homogeneous correction
		if (droplet) {
			for (unsigned int i = 0; i < NShells; i++) {
				if (rhoShellsT[i] >= 1.02 * rho_out) {
					rhoShellsT[i] -= rho_out;
				} else {
					rhoShellsT[i] = 0.0;
				}
			}
		} else {
			for (unsigned int i = 0; i < NShells; i++) {
				rhoShellsT[i] -= rho_out;
				// if (rhoShellsT[i] <= 0.998 * rho_out) {
				// 	rhoShellsT[i] -= rho_out;
				// } else {
				// 	rhoShellsT[i] = 0.0;
				// }
			}
		}

		// U Correction of homogeneous system
		// TODO: multi component
		// rhovtot = 0.0
		// DO i=1,AnzKomp,1
		//   rhovtot = rhovtot + rhov(i)
		// }
		UpotKorrLJ = 0.0;
		VirialKorrLJ = 0.0;
		for (unsigned ci = 0; ci < numComp; ++ci) {
			for (unsigned cj = 0; cj < numComp; ++cj) {
				ParaStrm& params = _domain->getComp2Params()(ci, cj);
				params.reset_read();
				for (unsigned si = 0; si < numLJ[ci]; ++si) {  // Long Range Correction for Lennard-Jones sites
					for (unsigned sj = 0; sj < numLJ[cj]; ++sj) {
						double eps24;
						double sig2;
						double shift6;
						params >> eps24;
						params >> sig2;
						params >> shift6;
						double tau1 = eLong[numLJSum2[ci] + si];
						double tau2 = eLong[numLJSum2[cj] + sj];
						// double tau1 = 0.0;
						// double tau2 = 0.0;
						if ((tau1 == 0.0) && (tau2 == 0.0)) {
							UpotKorrLJ =
								UpotKorrLJ + rho_out * (eps24 / 6.0) * (TICCu(-6, rc, sig2) - TICCu(-3, rc, sig2));
							VirialKorrLJ =
								VirialKorrLJ + rho_out * (eps24 / 6.0) * (TICCp(-6, rc, sig2) - TICCp(-3, rc, sig2));
						} else if ((tau1 == 0.0) || (tau2 == 0.0)) {
							double tau = max(tau1, tau2);  // Filtert den Wert Null raus
							UpotKorrLJ = UpotKorrLJ + rho_out * (eps24 / 6.0) *
														  (TICSu(-6, rc, sig2, tau) - TICSu(-3, rc, sig2, tau));
							VirialKorrLJ = VirialKorrLJ + rho_out * (eps24 / 6.0) *
															  (TICSp(-6, rc, sig2, tau) - TICSp(-3, rc, sig2, tau));
						} else if ((tau1 != 0.0) && (tau2 != 0.0)) {
							UpotKorrLJ =
								UpotKorrLJ + rho_out * (eps24 / 6.0) *
												 (TISSu(-6, rc, sig2, tau1, tau2) - TISSu(-3, rc, sig2, tau1, tau2));
							VirialKorrLJ =
								VirialKorrLJ + rho_out * (eps24 / 6.0) *
												   (TISSp(-6, rc, sig2, tau1, tau2) - TISSp(-3, rc, sig2, tau1, tau2));
						}
					}
				}
				UpotKorrLJ *= 2.0 * M_PI;
				VirialKorrLJ *= 2.0 * M_PI / 3.0;
			}
		}

		// Alternative way for one component
		double UCORR =  rho_out*(8./3.)*M_PI*(1./(3.*std::pow(rc,9))-1./std::pow(rc,3));
    	double PCORR = rho_out*(16./3.)*M_PI*(2./(3.*std::pow(rc,9))-1./std::pow(rc,3));

		global_log->info() << "[Long Range Correction] Homogeneous term: rho_out = " << rho_out << " UpotKorrLJ = " << UpotKorrLJ << " ; VirialKorrLJ = " << VirialKorrLJ << std::endl;
		global_log->info() << "[Long Range Correction] Alt. homog. term: rho_out = " << rho_out << " UpotKorrLJ = " << UCORR      << " ; VirialKorrLJ = " << PCORR << std::endl;

		// Korrektur je Schale
		std::fill(UShells_Mean.begin(), UShells_Mean.end(), 0.0);
		std::fill(UShells_Mean_global.begin(), UShells_Mean_global.end(), 0.0);
		std::fill(FShells_Mean.begin(), FShells_Mean.end(), 0.0);
		std::fill(FShells_Mean_global.begin(), FShells_Mean_global.end(), 0.0);
		std::fill(PNShells_Mean.begin(), PNShells_Mean.end(), 0.0);
		std::fill(PNShells_Mean_global.begin(), PNShells_Mean_global.end(), 0.0);
		std::fill(PTShells_Mean.begin(), PTShells_Mean.end(), 0.0);
		std::fill(PTShells_Mean_global.begin(), PTShells_Mean_global.end(), 0.0);

		double rlow, rlow2, rlowInv, rlowInv2, rdash2, rdashInv, UCorrTemp, rdash, rdashInv2, FCorrTemp, PNCorrTemp,
			PTCorrTemp;
		double UCorrShells = 0.0;
		double FCorrShells = 0.0;
		double PNCorrShells = 0.0;
		double PTCorrShells = 0.0;
		double ksi2 = 0.0;

		for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid();
			 ++tempMol) {
			unsigned long molID = tempMol->getID();
			UCorrShells = 0.0;
			FCorrShells = 0.0;
			PNCorrShells = 0.0;
			PTCorrShells = 0.0;
			ksi2 = ksi[molID] * ksi[molID];
			for (unsigned ci = 0; ci < numComp; ++ci) {
				for (unsigned cj = 0; cj < numComp; ++cj) {
					ParaStrm& params = _domain->getComp2Params()(ci, cj);
					params.reset_read();
					for (unsigned si = 0; si < numLJ[ci]; ++si) {  // Long Range Correction for Lennard-Jones sites
						for (unsigned sj = 0; sj < numLJ[cj]; ++sj) {
							double eps24;
							double sigma;
							double sig2;
							double shift6;
							double eps;
							params >> eps24;
							params >> sig2;
							params >> shift6;
							// sigma = sqrt(sig2);
							// eps = eps24 / 24;
							double tau1 = eLong[numLJSum2[ci] + si];
							double tau2 = eLong[numLJSum2[cj] + sj];
							// double tau1 = 0.5;
							// double tau2 = 0;
							double sigma6 = sig2 * sig2 * sig2;
							double factorU = -M_PI * _drShells * eps24 * sigma6 / (6. * ksi[molID]);
							double factorF = -factorU / ksi[molID];
							double factorP = 0.5 * factorF / ksi[molID];
							if ((tau1 == 0.0) && (tau2 == 0.0)) {  // Center-Center
								for (unsigned long j = 1; j < (lowerS[molID] + 1);
									 j++) {  // Loop over Shells with smaller Radius
									if (rhoShellsT[j - 1] != 0.0) {
										rlow = ksi[molID] - RShells[j - 1];
										rlowInv = 1. / rlow;
										rlowInv2 = rlowInv * rlowInv;
										rdashInv = 1. / (RShells[j - 1] + ksi[molID]);
										UCorrTemp = sigma6 * 0.2 * (pow(rdashInv, 10) - pow(rlowInv, 10)) -
													0.5 * (pow(rdashInv, 4) - pow(rlowInv, 4));
										if (rlow < rcmax) {
											rdash = min(rcmax, (RShells[j - 1] + ksi[molID]));
											rdashInv = 1. / rdash;
											rdashInv2 = rdashInv * rdashInv;
											FCorrTemp =
												sigma6 * ((6. / 5. * rdash * rdash + ksi2 - RShells2[j - 1]) *
															  pow(rdashInv2, 6) -
														  (6. / 5. * rlow * rlow + ksi2 - RShells2[j - 1]) *
															  pow(rlowInv2, 6)) -
												(1.5 * rdash * rdash + ksi2 - RShells2[j - 1]) * pow(rdashInv2, 3) +
												(1.5 * rlow * rlow + ksi2 - RShells2[j - 1]) * pow(rlowInv2, 3);
											FCorrShells =
												FCorrShells + FCorrTemp * factorF * rhoShellsT[j - 1] * RShells[j - 1];
											PNCorrTemp =  // 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
												-3. * (rdashInv2 - rlowInv2) +
												2. * (ksi2 - RShells2[j - 1]) *
													(  // 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
														-1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2))) +
												pow((ksi2 - RShells2[j - 1]), 2) *
													(  // sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6))
														-(pow(rdashInv2, 3) - pow(rlowInv2, 3)));
											PNCorrShells = PNCorrShells +
														   PNCorrTemp * factorP * rhoShellsT[j - 1] * RShells[j - 1];
											PTCorrTemp = 6. / 5. * sigma6 * (pow(rdashInv2, 5) - pow(rlowInv2, 5)) -
														 1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2));
											PTCorrShells = PTCorrShells + 4. * ksi2 * PTCorrTemp * factorP *
																			  rhoShellsT[j - 1] * RShells[j - 1];
										}
										UCorrShells =
											UCorrShells + UCorrTemp * factorU * rhoShellsT[j - 1] * RShells[j - 1];
									}
								}
								for (unsigned long j = upperS[molID]; j < (NShells + 1); j++) {
									if (rhoShellsT[j - 1] != 0.0) {
										rlow = RShells[j - 1] - ksi[molID];
										rlowInv = 1 / rlow;
										rlowInv2 = rlowInv * rlowInv;
										rdashInv = 1 / (RShells[j - 1] + ksi[molID]);
										UCorrTemp = sigma6 * 0.2 * (pow(rdashInv, 10) - pow(rlowInv, 10)) -
													0.5 * (pow(rdashInv, 4) - pow(rlowInv, 4));
										if (rlow < rcmax) {
											rdash = min(rcmax, (RShells[j - 1] + ksi[molID]));
											rdashInv = 1 / rdash;
											rdashInv2 = rdashInv * rdashInv;
											FCorrTemp =
												sigma6 * ((6. / 5. * rdash * rdash + ksi2 - RShells2[j - 1]) *
															  pow(rdashInv2, 6) -
														  (6. / 5. * rlow * rlow + ksi2 - RShells2[j - 1]) *
															  pow(rlowInv2, 6)) -
												(1.5 * rdash * rdash + ksi2 - RShells2[j - 1]) * pow(rdashInv2, 3) +
												(1.5 * rlow * rlow + ksi2 - RShells2[j - 1]) * pow(rlowInv2, 3);
											FCorrShells =
												FCorrShells + FCorrTemp * factorF * rhoShellsT[j - 1] * RShells[j - 1];
											PNCorrTemp =  // 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
												-3. * (rdashInv2 - rlowInv2) +
												2. * (ksi2 - RShells2[j - 1]) *
													(  // 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
														-1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2))) +
												pow((ksi2 - RShells2[j - 1]), 2) *
													(  // sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6))
														-(pow(rdashInv2, 3) - pow(rlowInv2, 3)));
											PNCorrShells = PNCorrShells +
														   PNCorrTemp * factorP * rhoShellsT[j - 1] * RShells[j - 1];
											PTCorrTemp = 6. / 5. * sigma6 * (pow(rdashInv2, 5) - pow(rlowInv2, 5)) -
														 1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2));
											PTCorrShells = PTCorrShells + 4. * ksi2 * PTCorrTemp * factorP *
																			  rhoShellsT[j - 1] * RShells[j - 1];
										}
										UCorrShells =
											UCorrShells + UCorrTemp * factorU * rhoShellsT[j - 1] * RShells[j - 1];
									}
								}
								for (unsigned long j = interS[molID]; j < upperS[molID];
									 j++) {  // Loop over partly contributing Shells
									if (rhoShellsT[j - 1] != 0.0) {
										rlow = rc;
										rlowInv = 1 / rlow;
										rlowInv2 = rlowInv * rlowInv;
										rdashInv = 1 / (RShells[j - 1] + ksi[molID]);
										UCorrTemp = sigma6 * 0.2 * (pow(rdashInv, 10) - pow(rlowInv, 10)) -
													0.5 * (pow(rdashInv, 4) - pow(rlowInv, 4));
										if (rlow < rcmax) {
											rdash = min(rcmax, (RShells[j - 1] + ksi[molID]));
											rdashInv = 1 / rdash;
											rdashInv2 = rdashInv * rdashInv;
											FCorrTemp =
												sigma6 * ((6. / 5. * rdash * rdash + ksi2 - RShells2[j - 1]) *
															  pow(rdashInv2, 6) -
														  (6. / 5. * rlow * rlow + ksi2 - RShells2[j - 1]) *
															  pow(rlowInv2, 6)) -
												(1.5 * rdash * rdash + ksi2 - RShells2[j - 1]) * pow(rdashInv2, 3) +
												(1.5 * rlow * rlow + ksi2 - RShells2[j - 1]) * pow(rlowInv2, 3);
											FCorrShells =
												FCorrShells + FCorrTemp * factorF * rhoShellsT[j - 1] * RShells[j - 1];
											PNCorrTemp =  // 1.5*sigma6 * (pow(rdashInv2,4) - pow(rlowInv2,4))
												-3. * (rdashInv2 - rlowInv2) +
												2. * (ksi2 - RShells2[j - 1]) *
													(  // 6/5*sigma6 * (pow(rdashInv2,5) - pow(rlowInv2,5))
														-1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2))) +
												pow((ksi2 - RShells2[j - 1]), 2) *
													(  // sigma6 * ( pow(rdashInv2,6) - pow(rlowInv2,6))
														-(pow(rdashInv2, 3) - pow(rlowInv2, 3)));
											PNCorrShells = PNCorrShells +
														   PNCorrTemp * factorP * rhoShellsT[j - 1] * RShells[j - 1];
											PTCorrTemp = 6. / 5. * sigma6 * (pow(rdashInv2, 5) - pow(rlowInv2, 5)) -
														 1.5 * (pow(rdashInv2, 2) - pow(rlowInv2, 2));
											PTCorrShells = PTCorrShells + 4. * ksi2 * PTCorrTemp * factorP *
																			  rhoShellsT[j - 1] * RShells[j - 1];
										}
										UCorrShells =
											UCorrShells + UCorrTemp * factorU * rhoShellsT[j - 1] * RShells[j - 1];
									}
								}
							} else if ((tau1 == 0.0) || (tau2 == 0.0)) {  // Center-Site
								std::cout << "weird center site" << std::endl;
							} else if ((tau1 != 0.0) && (tau2 != 0.0)) {  // Site-Site
								std::cout << "weird site site" << std::endl;
							} else {
								std::cout << "weird since else and not cc" << std::endl;
							}
						}
					}
				}
				PTCorrShells -= PNCorrShells;
				UShells_Mean[PartShells[molID]] += UCorrShells;
				FShells_Mean[PartShells[molID]] += FCorrShells;
				PNShells_Mean[PartShells[molID]] += PNCorrShells;
				PTShells_Mean[PartShells[molID]] += PTCorrShells;
			}
		}

		// Distribution of Shell Corrections to every node
		_domainDecomposition->collCommInit(4 * NShells);
		for (unsigned i = 0; i < NShells; i++) {
			_domainDecomposition->collCommAppendDouble(UShells_Mean[i]);
			_domainDecomposition->collCommAppendDouble(FShells_Mean[i]);
			_domainDecomposition->collCommAppendDouble(PNShells_Mean[i]);
			_domainDecomposition->collCommAppendDouble(PTShells_Mean[i]);
		}
		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i = 0; i < NShells; i++) {
			UShells_Mean_global[i] = _domainDecomposition->collCommGetDouble();
			FShells_Mean_global[i] = _domainDecomposition->collCommGetDouble();
			PNShells_Mean_global[i] = _domainDecomposition->collCommGetDouble();
			PTShells_Mean_global[i] = _domainDecomposition->collCommGetDouble();
		}
		_domainDecomposition->collCommFinalize();

		for (unsigned int i = 0; i < NShells; i++) {
			if (rhoShellsTemp_global[i] != 0.0) {
				UShells_Mean_global[i] /= (rhoShellsTemp_global[i] * VShells[i]);
				FShells_Mean_global[i] /= (rhoShellsTemp_global[i] * VShells[i]);
				PNShells_Mean_global[i] /= (rhoShellsTemp_global[i] * VShells[i]);
				PTShells_Mean_global[i] /= (rhoShellsTemp_global[i] * VShells[i]);
			}
		}
		// Only Root writes to files
		if (rank == 0) {
			if (simstep % writeFreq == 0) {  // reduced from 100 to 100000, should be enough, no?

				// calculating and writing thermoData:

				// int inside_from = NShells / 10;
				// int inside_to = NShells / 6;
				// int outside_from = NShells - NShells / 20 - 1;
				// int outside_to = NShells - 1;

				// double pInside = 0;
				// double rhoInside = 0;
				// for (unsigned int i = inside_from; i < inside_to; i++) {
				// 	pInside += _T * rhoShellsAvg_global[i];
				// 	pInside += (VirShells_T_global[i] + VirShells_N_global[i]) / 2;
				// 	rhoInside += rhoShellsAvg_global[i];
				// }
				// pInside /= (inside_to - inside_from);
				// rhoInside /= (inside_to - inside_from);

				// double pOutside = 0;
				// double rhoOutside = 0;
				// for (unsigned int i = outside_from; i < outside_to; i++) {
				// 	pOutside += _T * rhoShellsAvg_global[i] + VirShells_N_global[i];
				// 	pOutside += _T * rhoShellsAvg_global[i] + VirShells_T_global[i];
				// 	rhoOutside += rhoShellsAvg_global[i];
				// }
				// pOutside /= (outside_to - outside_from) * 2;
				// rhoOutside /= outside_to - outside_from;

				// double rhoDiff_Avg = 0;
				// double pDiff_Avg = 0;
				// if (rhoInside > rhoOutside) {              // bubble
				// 	rhoDiff_Avg = rhoOutside - rhoInside;  // rho_vap - rho_liq
				// 	pDiff_Avg = pInside - pOutside;        // p_liq - p_vap
				// } else {
				// 	rhoDiff_Avg = rhoInside - rhoOutside;  // rho_vap - rho_liq
				// 	pDiff_Avg = pOutside - pInside;        // p_liq - p_vap
				// }
				// double pDiff2_Avg = pow(pInside - pOutside, 2);

				// double dpN_Avg[NShells];
				// memset(dpN_Avg, 0, NShells * sizeof(double));
				// double drho_Avg[NShells];
				// memset(drho_Avg, 0, NShells * sizeof(double));
				// double R_e3 = 0.;
				// double R_e = 0;

				// // double integral_term_Avg[NShells] = { 0 };
				// // double gamma_Avg_iterative[NShells] = { 0 };
				// double gamma_integral_Avg[NShells];
				// memset(gamma_integral_Avg, 0, NShells * sizeof(double));
				// double gamma_Avg[NShells];
				// memset(gamma_Avg, 0, NShells * sizeof(double));

				// // std::cout << "--------------------------------------------------------- "<< std::endl;
				// // std::cout << "--------------Simstep = "<< simstep << "----------------------"<< std::endl;
				// gamma_integral_Avg[0] = 0.;
				// gamma_Avg[0] = 0.;
				// // calculation of gamma / R_e:
				// for (unsigned i = 1; i < NShells; i++) {  // not considering index 0!
				// 	dpN_Avg[i] = (_T * rhoShellsAvg_global[i] + VirShells_N_global[i]) -
				// 				 (_T * rhoShellsAvg_global[i - 1] + VirShells_N_global[i - 1]);
				// 	drho_Avg[i] = rhoShellsAvg_global[i] - rhoShellsAvg_global[i - 1];

				// 	R_e3 += RShells3[i] * drho_Avg[i];
				// 	// //------- TODO: remove these two, once everything checks out:
				// 	// integral_term_Avg[i] = integral_term_Avg[i-1] + RShells3[i] * dpN_Avg[i];
				// 	// gamma_Avg_iterative[i] =  pow(-(integral_term_Avg[i] * pDiff2_Avg)/8., 1/3.);
				// 	// // ----------------- //
				// 	gamma_integral_Avg[i] = gamma_integral_Avg[i - 1] + RShells3[i] * dpN_Avg[i];
				// 	gamma_Avg[i] = -std::cbrt((pDiff2_Avg * gamma_integral_Avg[i]) / 8.);
				// 	// std::cout << "i = "<< i << ""<< std::endl;
				// 	// std::cout << "dpN_Avg[ "<<i<<"] = " << dpN_Avg[i] << ""<< std::endl;
				// 	// std::cout << "drho_Avg [ "<<i<<"] = " << drho_Avg[i] << ""<< std::endl;
				// 	// std::cout << "R_e3 = "<< R_e3 << ""<< std::endl;
				// 	// std::cout << "gamma_integral_Avg = "<< gamma_integral_Avg << std::endl;
				// }

				// // std::cout << "--------------------------------------------------------- "<< std::endl;
				// // std::cout << "rhoDiff = "<< rhoDiff_Avg << ""<< std::endl;
				// // std::cout << "pDiff   = "<< pDiff_Avg << "; pOutside = "<<pOutside <<"; pInside="<< pInside<<
				// // std::endl;
				// R_e3 /= rhoDiff_Avg;
				// // std::cout << "R_e3    = "<< R_e3 << ""<< std::endl;
				// R_e = std::cbrt(R_e3);
				// // std::cout << "R_e     = "<< R_e << ""<< std::endl;
				// // std::cout << "gamma   = "<< gamma_Avg << ""<< std::endl;
				// // gamma_Avg[i] =  -std::cbrt((pDiff2_Avg*gamma_integral_Avg[i])/8.);

				// double R_gamma = 2 * gamma_Avg[NShells - 2] / pDiff_Avg;
				// // std::cout << "R_gamma = "<< R_gamma << ""<< std::endl;

				// ofstream outfilestreamThermData(filenameThermData, ios::app);
				// outfilestreamThermData << std::setw(24) << simstep << "";
				// // outfilestreamThermData << std::setw(24) << gamma_Avg_iterative[NShells-1] << "";
				// outfilestreamThermData << std::setw(24) << gamma_Avg[NShells - 3] << "";
				// outfilestreamThermData << std::setw(24) << gamma_Avg[NShells - 2] << "";
				// outfilestreamThermData << std::setw(24) << gamma_Avg[NShells - 1] << "";
				// outfilestreamThermData << std::setw(24) << R_gamma << "";
				// outfilestreamThermData << std::setw(24) << R_e << "";
				// outfilestreamThermData << std::setw(24) << R_e - R_gamma << "";
				// outfilestreamThermData << std::setw(24) << rhoInside << "";
				// outfilestreamThermData << std::setw(24) << rhoOutside << "";
				// outfilestreamThermData << std::setw(24) << pInside << "";
				// outfilestreamThermData << std::setw(24) << pOutside << "";
				// outfilestreamThermData << std::setw(24) << inside_from << "";
				// outfilestreamThermData << std::setw(24) << inside_to << "";
				// outfilestreamThermData << std::setw(24) << outside_from << "";
				// outfilestreamThermData << std::setw(24) << outside_to << "";
				// outfilestreamThermData << std::setw(24) << nSamples << "";
				// outfilestreamThermData << std::endl;

				// outfilestreamThermData.close();

				// output for GlobalCorrs:
				filenameGlobalCorrs = _outputPrefix + "_globalCorrections_" + std::to_string(simstep) + ".csv";
				ofstream outfilestreamGlobalCorrs(filenameGlobalCorrs, ios::out);
				outfilestreamGlobalCorrs << std::setw(24) << "Radius";
				outfilestreamGlobalCorrs << std::setw(24) << "RhoShells_Avg";
				outfilestreamGlobalCorrs << std::setw(24) << "UShells_Mean";
				outfilestreamGlobalCorrs << std::setw(24) << "FShells_Mean";
				outfilestreamGlobalCorrs << std::setw(24) << "PNShells";  // Spherical Teil der Korrektur, im aktuellen Zeitschritt ???
				outfilestreamGlobalCorrs << std::setw(24) << "PTShells";  // Spherical Teil der Korrektur, im aktuellen Zeitschritt ???
				outfilestreamGlobalCorrs << std::setw(24) << "PNShells_Avg";  // Spherical Teil der Korrektur, gemittelt ???
				outfilestreamGlobalCorrs << std::setw(24) << "PTShells_Avg";  // Spherical Teil der Korrektur, gemittelt ???
				outfilestreamGlobalCorrs << std::setw(24) << "P_N_Avg";
				outfilestreamGlobalCorrs << std::setw(24) << "P_T_Avg";
				outfilestreamGlobalCorrs << std::setw(24) << "T_Avg";
				outfilestreamGlobalCorrs << std::setw(24) << "VirialKorrLJ";  // Homogener Teil der Korrektur
				outfilestreamGlobalCorrs << std::setw(24) << "rhoShellsT";
				outfilestreamGlobalCorrs << std::setw(24) << "rhoShells_global";
				// outfilestreamGlobalCorrs << std::setw(24) << "gamma";

				outfilestreamGlobalCorrs << std::endl;
				for (unsigned int i = 0; i < NShells; i++) {
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10) << RShells[i]
											 << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << rhoShellsAvg_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << UShells_Mean_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << FShells_Mean_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << PNShells_Mean_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << PTShells_Mean_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << PNShellsCorr_Avg_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << PTShellsCorr_Avg_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << _T * rhoShellsAvg_global[i] + VirShells_N_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << _T * rhoShellsAvg_global[i] + VirShells_T_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << TShellsAvg_global[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << VirialKorrLJ << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << rhoShellsT[i] << "";
					outfilestreamGlobalCorrs << std::setw(24)
											 << std::setprecision(std::numeric_limits<double>::digits10)
											 << rhoShells_global[i] << "";
					// outfilestreamGlobalCorrs << std::setw(24) <<
					// std::setprecision(std::numeric_limits<double>::digits10) << gamma_Avg_iterative[i]<< "" ;
					outfilestreamGlobalCorrs << std::endl;
				}
				outfilestreamGlobalCorrs.close();
			}
		}
		nSamples++;
	}

	// Lesen der Korrekturen in jedem Zeitschritt
	
	// Variablen für globale Korrektur in Domain
	double UCorrSum = 0.0;
	double VirCorrSum = 0.0;

	for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid();
		 ++tempMol) {
		unsigned long molID = tempMol->getID();

		// Get virial from particle data which was calculated in potforce.h
		// The "normal" cartisian virial is not relevant here
		VirShells_N[PartShells[molID]] += tempMol->ViN(); // * 0.5; // *.5 to keep consistency with legacy (potforce.h didnt calculate exact virn)
		VirShells_T[PartShells[molID]] += tempMol->ViT();//*.25;// *.25 to keep consistency with legacy (potforce.h didnt calculate exact virt)
		// Apply LRC for each shell, first spherical specific part, then homogeneous part
		if (!disableLRC) {
			VirShells_N[PartShells[molID]] -= 0.5 * PNShells_Mean_global[PartShells[molID]];
			VirShells_T[PartShells[molID]] -= 0.25 * PTShells_Mean_global[PartShells[molID]];
			VirShells_N[PartShells[molID]] -= VirialKorrLJ;
			VirShells_T[PartShells[molID]] -= VirialKorrLJ;
		}
		FcorrX[molID] = FShells_Mean_global[PartShells[molID]] * FcorrX[molID] / ksi[molID];
		FcorrY[molID] = FShells_Mean_global[PartShells[molID]] * FcorrY[molID] / ksi[molID];
		FcorrZ[molID] = FShells_Mean_global[PartShells[molID]] * FcorrZ[molID] / ksi[molID];

		/* // Reset virial NOW DONE IN simulation.cpp
		double zero3[] = {.0,.0,.0};
		tempMol->setVi(zero3);
		tempMol->setViN(0.0);
		tempMol->setViT(0.0); 
		*/

		PNShellsCorr_Avg[PartShells[molID]] -= 0.5 * PNShells_Mean_global[PartShells[molID]];
		PTShellsCorr_Avg[PartShells[molID]] -= 0.25 * PTShells_Mean_global[PartShells[molID]];

		// Accumulate global corrections
		UCorrSum += UShells_Mean_global[PartShells[molID]];
		VirCorrSum -= 0.5 * PNShells_Mean_global[PartShells[molID]];
		VirCorrSum -= 0.5 * PTShells_Mean_global[PartShells[molID]];
		VirCorrSum -= 3. * VirialKorrLJ;
	}

	//	std::cout << "Local rank " << rank << " UCorrSum: " << UCorrSum << " VirCorrSum: " << VirCorrSum <<  std::endl;

	// Each particle was considered twice ?
	UCorrSum = 0.5 * UCorrSum;

	// // Distribution of the Force, Energy and Virial to every Node
	_domainDecomposition->collCommInit(4 * NShells);
	for (unsigned i = 0; i < NShells; i++) {
		_domainDecomposition->collCommAppendDouble(PNShellsCorr_Avg[i]);
		_domainDecomposition->collCommAppendDouble(PTShellsCorr_Avg[i]);
		_domainDecomposition->collCommAppendDouble(VirShells_N[i]);
		_domainDecomposition->collCommAppendDouble(VirShells_T[i]);
	}
	_domainDecomposition->collCommAllreduceSum();
	for (unsigned i = 0; i < NShells; i++) {
		PNShellsCorr_Avg_global[i] = _domainDecomposition->collCommGetDouble();
		PTShellsCorr_Avg_global[i] = _domainDecomposition->collCommGetDouble();
		VirShells_N_global[i] = _domainDecomposition->collCommGetDouble();
		VirShells_T_global[i] = _domainDecomposition->collCommGetDouble();
	}
	_domainDecomposition->collCommFinalize();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendDouble(UCorrSum);
	_domainDecomposition->collCommAllreduceSum();
	double UCorrSum_global = _domainDecomposition->collCommGetDouble();
	_domainDecomposition->collCommFinalize();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendDouble(VirCorrSum);
	_domainDecomposition->collCommAllreduceSum();
	double VirCorrSum_global = _domainDecomposition->collCommGetDouble();
	_domainDecomposition->collCommFinalize();

	for (unsigned int i = 0; i < NShells; i++) {
		PNShellsCorr_Avg_global[i] /= (3 * VShells[i] * (simstep + 1));
		PTShellsCorr_Avg_global[i] /= (3 * VShells[i] * (simstep + 1));
		VirShells_N_global[i] /= (VShells[i] * (simstep + 1));
		VirShells_T_global[i] /= (VShells[i] * (simstep + 1));
	}

	// Add homogeneous part to upot correction
	// For virial for domain correction this was done above (TODO move here)
	UCorrSum_global += globalNumMols * UpotKorrLJ;

	if (!disableLRC) {
		// Apply corrections to each particle
		for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid();
			 ++tempMol) {
			unsigned long molID = tempMol->getID();
			double Fa[3] = {0.0, 0.0, 0.0};
			Fa[0] = FcorrX[molID];
			Fa[1] = FcorrY[molID];
			Fa[2] = FcorrZ[molID];
			tempMol->Fadd(Fa);

			// tempMol->VirNadd(-0.5 * PNShells_Mean_global[PartShells[molID]]);
			// tempMol->VirTadd(-0.25 * PTShells_Mean_global[PartShells[molID]]);
			// tempMol->VirNadd(-VirialKorrLJ);
			// tempMol->VirTadd(-VirialKorrLJ);
			const double virNCorr = -0.5 * PNShells_Mean_global[PartShells[molID]] - VirialKorrLJ;
			const double virTCorr = -0.25 * PTShells_Mean_global[PartShells[molID]] - VirialKorrLJ;
			tempMol->setViNConstCorr(virNCorr);
			tempMol->setViTConstCorr(virTCorr);
			tempMol->Uadd(UShells_Mean_global[PartShells[molID]]);
		}
		// Apply corrections to global values in domain class
		_domain->setUpotCorr(UCorrSum_global);
		_domain->setVirialCorr(VirCorrSum_global);
	} else {
		_domain->setUpotCorr(0.);
		_domain->setVirialCorr(0.);
	}
}

// void Spherical::centerCenter(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj) {
// 	double a = 1;
// }

// void Spherical::centerSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
// 	// double sig2=sig*sig;
// 	// double sig3=sig2*sig;
// 	// double t = eLong[numLJSum2[ci]+si] + eLong[numLJSum2[cj]+sj]; // one of them is equal to zero.
// 	global_log->error() << "[LongRangeCorrection] Center-Site correction not yet supported. Program exit ..." <<
// std::endl;
//     Simulation::exit(-1);
// }

// void Spherical::siteSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
// 	double a = 1;
// 	// double sig2=sig*sig;
// 	// double sig3=sig2*sig;
// 	// double sig4=sig2*sig2;
// 	// double t1 = eLong[numLJSum2[ci]+si];
// 	// double t2 = eLong[numLJSum2[cj]+sj];
// 	// double tP = t1 + t2; // tau+
// 	// double tM = t1 - t2; // tau-

// 	// global_log->error() << "[LongRangeCorrection] Site-Site correction not yet supported. Program exit ..." <<
// std::endl;
//     // Simulation::exit(-1);
// }

void Spherical::writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstepI) { double a = 1; }

double Spherical::RhoP(double r, double rhov, double rhol, double D0, double R0) {
	return 0.5 * (rhol + rhov) - 0.5 * (rhol - rhov) * tanh(2 * (r - R0) / D0);
}

// double Spherical::SICCu(int n, double r) {

//   return -( pow(r,(2*n + 2)) ) /( (n+1) );

// }

double Spherical::SICSu(int n, double r, double tau) {
	return (pow((r + tau), (2 * n + 3)) - pow((r - tau), (2 * n + 3))) / (4. * tau * (n + 1) * (2 * n + 3));
}

double Spherical::CS(int n, double r, double tau) {
	return (pow((r + tau), (2 * n + 2)) - pow((r - tau), (2 * n + 2))) / (4. * r * tau * (n + 1));
}

double Spherical::SISSu(int n, double r, double tau1, double tau2) {
	double tauPlus = tau1 + tau2;
	double tauMinus = tau1 - tau2;

	return (pow((r + tauPlus), (2 * n + 4)) - pow((r + tauMinus), (2 * n + 4)) - pow((r - tauMinus), (2 * n + 4)) +
			pow((r - tauPlus), (2 * n + 4))) /
		   (8. * tau1 * tau2 * (n + 1) * (2 * n + 3) * (2 * n + 4));
}

double Spherical::SS(int n, double r, double tau1, double tau2) {
	double tauPlus = tau1 + tau2;
	double tauMinus = tau1 - tau2;

	return (pow((r + tauPlus), (2 * n + 3)) - pow((r + tauMinus), (2 * n + 3)) - pow((r - tauMinus), (2 * n + 3)) +
			pow((r - tauPlus), (2 * n + 3))) /
		   (8. * tau1 * tau2 * r * (n + 1) * (2 * n + 3));
}

double Spherical::SSLN(double r, double tau1, double tau2) {
	double tauPlus = tau1 + tau2;
	double tauMinus = tau1 - tau2;

	return ((pow((r + tauPlus), (-1)) - pow((r + tauMinus), (-1)) - pow((r - tauMinus), (-1)) +
			 pow((r - tauPlus), (-1))) *
				r -
			log((r * r - tauPlus * tauPlus) / (r * r - tauMinus * tauMinus))) /
		   (48. * tau1 * tau2);
}

// Functions for homogeneous corrections

double Spherical::TICCu(int n, double rcutoff, double sigma2) {
	return -(pow(rcutoff, (2 * n + 3))) / (pow(sigma2, n) * (2 * n + 3));
}

double Spherical::TICSu(int n, double rcutoff, double sigma2, double tau) {
	return -(pow((rcutoff + tau), (2 * n + 3)) - pow((rcutoff - tau), (2 * n + 3))) * rcutoff /
			   (pow(4. * sigma2, n) * tau * (n + 1) * (2 * n + 3)) +
		   (pow((rcutoff + tau), (2 * n + 4)) - pow((rcutoff - tau), (2 * n + 4))) /
			   (4. * pow(sigma2, n) * tau * (n + 1) * (2 * n + 3) * (2 * n + 4));
}

double Spherical::TISSu(int n, double rcutoff, double sigma2, double tau1, double tau2) {
	double tauPlus = tau1 + tau2;
	double tauMinus = tau1 - tau2;
	return -(pow((rcutoff + tauPlus), (2 * n + 4)) - pow((rcutoff + tauMinus), (2 * n + 4)) -
			 pow((rcutoff - tauMinus), (2 * n + 4)) + pow((rcutoff - tauPlus), (2 * n + 4))) *
			   rcutoff / (8. * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3) * (2 * n + 4)) +
		   (pow((rcutoff + tauPlus), (2 * n + 5)) - pow((rcutoff + tauMinus), (2 * n + 5)) -
			pow((rcutoff - tauMinus), (2 * n + 5)) + pow((rcutoff - tauPlus), (2 * n + 5))) /
			   (8. * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3) * (2 * n + 4) * (2 * n + 5));
}

double Spherical::TICCp(int n, double rcutoff, double sigma2) { return 2 * n * TICCu(n, rcutoff, sigma2); }

double Spherical::TICSp(int n, double rcutoff, double sigma2, double tau) {
	return -(pow((rcutoff + tau), (2 * n + 2)) - pow((rcutoff - tau), (2 * n + 2))) * pow(rcutoff, (2)) /
			   (4. * pow(sigma2, n) * tau * (n + 1)) -
		   3. * TICSu(n, rcutoff, sigma2, tau);
}

double Spherical::TISSp(int n, double rcutoff, double sigma2, double tau1, double tau2) {
	double tauPlus = tau1 + tau2;
	double tauMinus = tau1 - tau2;
	return -(pow((rcutoff + tauPlus), (2 * n + 3)) - pow((rcutoff + tauMinus), (2 * n + 3)) -
			 pow((rcutoff - tauMinus), (2 * n + 3)) + pow((rcutoff - tauPlus), (2 * n + 3))) *
			   pow(rcutoff, 2) / (8. * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3)) -
		   3. * TISSu(n, rcutoff, sigma2, tau1, tau2);
}
