/*
 * CubicGridGeneratorInternal.cpp
 *
 *  Created on: Jun 9, 2017
 *      Author: seckler
 */

/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1

#include "CubicGridGeneratorInternal.h"

#include "Domain.h"
#include "IOHelpers.h"
#include "WrapOpenMP.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/Random.h"

#include <cmath>
#include <algorithm>

CubicGridGeneratorInternal::CubicGridGeneratorInternal() :
		_numMolecules(0), _binaryMixture(false)//, _RNG(0)
{
}

void CubicGridGeneratorInternal::readXML(XMLfileUnits& xmlconfig) {
	Log::global_log->info() << "------------------------------------------------------------------------" << std::endl;
	Log::global_log->info() << "CubicGridGeneratorInternal (CGG)" << std::endl;
	xmlconfig.getNodeValue("numMolecules", _numMolecules);
	if(_numMolecules)Log::global_log->info() << "numMolecules: " << _numMolecules << std::endl;
	double density = -1.;
	xmlconfig.getNodeValue("density", density);
	Log::global_log->info() << "density: " << density << std::endl;
	xmlconfig.getNodeValue("binaryMixture", _binaryMixture);
	Log::global_log->info() << "binaryMixture: " << _binaryMixture << std::endl;
	// setting both or none is not allowed!
	if((_numMolecules == 0 && density == -1.) || (_numMolecules != 0 && density != -1.) ){
		Log::global_log->error() << "Error in CubicGridGeneratorInternal: You have to set either density or numMolecules!" << std::endl;
		Simulation::exit(2341);
	}

	if(density != -1.){
		// density has been set
		if(density <= 0){
			Log::global_log->error()
					<< "Error in CubicGridGeneratorInternal: Density has to be positive and non-zero!"
					<< std::endl;
			Simulation::exit(2342);
		}
		double vol = 1.0;
		for (int d = 0; d < 3; ++d)
			vol *= global_simulation->getDomain()->getGlobalLength(d);
		_numMolecules = density * vol;
		Log::global_log->info() << "numMolecules: " << _numMolecules << std::endl;
	}
}

unsigned long CubicGridGeneratorInternal::readPhaseSpace(ParticleContainer *particleContainer, Domain *domain, DomainDecompBase *domainDecomp) {
	global_simulation->timers()->start("CUBIC_GRID_GENERATOR_INPUT");
	Log::global_log->info() << "Reading phase space file (CubicGridGenerator)." << std::endl;

	if(_numMolecules == 0){
		Log::global_log->error() << "Error in CubicGridGeneratorInternal: numMolecules is not set!"
				<< std::endl << "Please make sure to run readXML()!" << std::endl;
		Simulation::exit(2341);
	}

	// create a body centered cubic layout, by creating by placing the molecules on the
	// vertices of a regular grid, then shifting that grid by spacing/2 in all dimensions.

	std::array<double, 3> simBoxLength{};
	for (int d = 0; d < 3; ++d) {
		simBoxLength[d] = global_simulation->getDomain()->getGlobalLength(d);
	}

	std::array<unsigned long, 3> numMoleculesPerDim = determineMolsPerDimension(_numMolecules, simBoxLength);


	if (_binaryMixture) {
		global_simulation->getEnsemble()->getComponents()->at(1).updateMassInertia();
	}

	unsigned long int id = particleContainer->initCubicGrid(
		numMoleculesPerDim, simBoxLength, static_cast<size_t>(domainDecomp->getRank()) * mardyn_get_max_threads());

	Log::global_log->info() << "Finished reading molecules: 100%" << std::endl;

	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(id); //number of local molecules
	domainDecomp->collCommScanSum();
	unsigned long idOffset = domainDecomp->collCommGetUnsLong() - id;
	domainDecomp->collCommFinalize();
	// fix ID's to be unique:
	Log::global_log->info() << "CGG: ids" << std::endl;
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{

		for (auto mol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
			mol->setid(mol->getID() + idOffset);
		}
	}
	Log::global_log->info() << "CGG: ids done" << std::endl;
	//std::cout << domainDecomp->getRank()<<": #num local molecules:" << id << std::endl;
	//std::cout << domainDecomp->getRank()<<": offset:" << idOffset << std::endl;

	Log::global_log->info() << "CGG: remove momentum" << std::endl;
	IOHelpers::removeMomentum(particleContainer, *(global_simulation->getEnsemble()->getComponents()), domainDecomp);
	Log::global_log->info() << "CGG: momentum done" << std::endl;
	domain->evaluateRho(particleContainer->getNumberOfParticles(), domainDecomp);
	Log::global_log->info() << "Calculated Rho=" << domain->getglobalRho() << std::endl;
	global_simulation->timers()->stop("CUBIC_GRID_GENERATOR_INPUT");
	global_simulation->timers()->setOutputString("CUBIC_GRID_GENERATOR_INPUT", "Initial IO took:                 ");
	Log::global_log->info() << "Initial IO took:                 "
			<< global_simulation->timers()->getTime("CUBIC_GRID_GENERATOR_INPUT") << " sec" << std::endl;
	Log::global_log->info() << "------------------------------------------------------------------------" << std::endl;
	return id + idOffset;
}

std::array<unsigned long, 3> CubicGridGeneratorInternal::determineMolsPerDimension(
		unsigned long targetTotalNumMols,
		std::array<double, 3> boxLength) const {

	unsigned long long numMoleculesHalf = targetTotalNumMols / 2; // two grids

	double vol = 1.0;

	for (int d = 0; d < 3; ++d) {
		vol *= boxLength[d];
	}

	std::array<unsigned long, 3> ret;
	for (int d = 0; d < 3; ++d) {
		double L = boxLength[d];
		double frac = L * L * L / vol;
		double oneThird = 1.0 / 3.0;
		unsigned long answer = round(pow(numMoleculesHalf * frac, oneThird));

		mardyn_assert(answer >= 1);

		if (answer < 1) {
			Log::global_log->error() << "computed num Molecules along dimension " << d << ": " << answer << std::endl;
			Log::global_log->error() << "Should be larger than 1. Exiting." << std::endl;
			mardyn_exit(1);
		}

		ret[d] = answer;
	}

	// quality of approximation:
	unsigned long larger = std::max(ret[0] * ret[1] * ret[2] * 2, targetTotalNumMols);
	unsigned long smaller = std::min(ret[0] * ret[1] * ret[2] * 2, targetTotalNumMols);
	unsigned long diff = larger - smaller;

	// now iterate over {-1, 0, +1} for the three found numbers, to find best overall match
	std::array<long, 3> ret2;
	std::array<long, 3> temp;
	for (int d = 0; d < 3; ++d)
		temp[d] = static_cast<long>(ret[d]);

	for (long z = -1; z <= +1; ++z) {
		for (long y = -1; y <= +1; ++y) {
			for (long x = -1; x <= +1; ++x) {

				ret2[0] = temp[0] + x;
				ret2[1] = temp[1] + y;
				ret2[2] = temp[2] + z;

				unsigned long larger2 = std::max(static_cast<unsigned long>(ret2[0] * ret2[1] * ret2[2] * 2), targetTotalNumMols);
				unsigned long smaller2 = std::min(static_cast<unsigned long>(ret2[0] * ret2[1] * ret2[2] * 2), targetTotalNumMols);
				if (larger2 - smaller2 < diff) {
					diff = larger2 - smaller2;
					for (int d = 0; d < 3; ++d)
						ret[d] = static_cast<unsigned long>(ret2[d]);
				}
			}
		}
	}

	return ret;
}
