/*
 * CubicGridGeneratorInternal.cpp
 *
 *  Created on: Jun 9, 2017
 *      Author: seckler
 */

/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1

#include "CubicGridGeneratorInternal.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

#include "utils/Random.h"
#include "WrapOpenMP.h"

#include <cmath>
#include <algorithm>

CubicGridGeneratorInternal::CubicGridGeneratorInternal() :
		_numMolecules(0), _binaryMixture(false)//, _RNG(0)
{
}

void CubicGridGeneratorInternal::readXML(XMLfileUnits& xmlconfig) {
	global_log->info() << "------------------------------------------------------------------------" << std::endl;
	global_log->info() << "CubicGridGeneratorInternal (CGG)" << std::endl;
	xmlconfig.getNodeValue("numMolecules", _numMolecules);
	if(_numMolecules)global_log->info() << "numMolecules: " << _numMolecules << std::endl;
	double density = -1.;
	xmlconfig.getNodeValue("density", density);
	global_log->info() << "density: " << density << std::endl;
	xmlconfig.getNodeValue("binaryMixture", _binaryMixture);
	global_log->info() << "binaryMixture: " << _binaryMixture << std::endl;
	// setting both or none is not allowed!
	if((_numMolecules == 0 && density == -1.) || (_numMolecules != 0 && density != -1.) ){
		global_log->error() << "Error in CubicGridGeneratorInternal: You have to set either density or numMolecules!" << std::endl;
		Simulation::exit(2341);
	}

	if(density != -1.){
		// density has been set
		if(density <= 0){
			global_log->error()
					<< "Error in CubicGridGeneratorInternal: Density has to be positive and non-zero!"
					<< std::endl;
			Simulation::exit(2342);
		}
		double vol = 1.0;
		for (int d = 0; d < 3; ++d)
			vol *= global_simulation->getDomain()->getGlobalLength(d);
		_numMolecules = density * vol;
		global_log->info() << "numMolecules: " << _numMolecules << std::endl;
	}
}

unsigned long CubicGridGeneratorInternal::readPhaseSpace(ParticleContainer *particleContainer, Domain *domain, DomainDecompBase *domainDecomp) {
	global_simulation->timers()->start("CUBIC_GRID_GENERATOR_INPUT");
	Log::global_log->info() << "Reading phase space file (CubicGridGenerator)." << std::endl;

	if(_numMolecules == 0){
		global_log->error() << "Error in CubicGridGeneratorInternal: numMolecules is not set!"
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

	unsigned long int id = 0;

    id = particleContainer->initCubicGrid(numMoleculesPerDim, simBoxLength, domainDecomp->getRank()*mardyn_get_max_threads());

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
	removeMomentum(particleContainer, *(global_simulation->getEnsemble()->getComponents()), domainDecomp);
	Log::global_log->info() << "CGG: momentum done" << std::endl;
	domain->evaluateRho(particleContainer->getNumberOfParticles(), domainDecomp);
	Log::global_log->info() << "Calculated Rho=" << domain->getglobalRho() << std::endl;
	global_simulation->timers()->stop("CUBIC_GRID_GENERATOR_INPUT");
	global_simulation->timers()->setOutputString("CUBIC_GRID_GENERATOR_INPUT", "Initial IO took:                 ");
	Log::global_log->info() << "Initial IO took:                 "
			<< global_simulation->timers()->getTime("CUBIC_GRID_GENERATOR_INPUT") << " sec" << std::endl;
	global_log->info() << "------------------------------------------------------------------------" << std::endl;
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
			global_log->error() << "computed num Molecules along dimension " << d << ": " << answer << std::endl;
			global_log->error() << "Should be larger than 1. Exiting." << std::endl;
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

//bool CubicGridGeneratorInternal::addMolecule(double x, double y, double z, unsigned long id,
//		ParticleContainer* particleContainer) {
//	std::vector<double> velocity = getRandomVelocity(global_simulation->getEnsemble()->T());
//
//	//double orientation[4] = {1, 0, 0, 0}; // default: in the xy plane
//	// rotate by 30° along the vector (1/1/0), i.e. the angle bisector of x and y axis
//	// o = cos 30° + (1 1 0) * sin 15°
//	double orientation[4];
//	getOrientation(15, 10, orientation);
//
//	int componentType = 0;
//	if (_binaryMixture) {
//		componentType = randdouble(0, 1.999999);
//	}
//
//	double I[3] = { 0., 0., 0. };
//	I[0] = global_simulation->getEnsemble()->getComponents()->at(0).I11();
//	I[1] = global_simulation->getEnsemble()->getComponents()->at(0).I22();
//	I[2] = global_simulation->getEnsemble()->getComponents()->at(0).I33();
//	/*****  Copied from animake - initialize anular velocity *****/
//	double w[3];
//	for (int d = 0; d < 3; d++) {
//		w[d] = (I[d] == 0) ?
//				0.0 : ((randdouble(0, 1) > 0.5) ? 1 : -1) * sqrt(2.0 * randdouble(0, 1) * global_simulation->getEnsemble()->T() / I[d]);
//		double fs_2_mardyn = 0.030619994;
//		w[d] = w[d] * fs_2_mardyn;
//	}
//	/************************** End Copy **************************/
//
//	Molecule m(id, &(global_simulation->getEnsemble()->getComponents()->at(componentType)), x, y, z, // position
//			velocity[0], -velocity[1], velocity[2], // velocity
//			orientation[0], orientation[1], orientation[2], orientation[3], w[0], w[1], w[2]);
//	return particleContainer->addParticle(m);
//}

void CubicGridGeneratorInternal::removeMomentum(ParticleContainer* particleContainer,
		const std::vector<Component>& components, DomainDecompBase* domainDecomp) {
	double mass_sum = 0.;
	//double momentum_sum[3] = { 0., 0., 0. };
	double momentum_sum0 = 0.;
	double momentum_sum1 = 0.;
	double momentum_sum2 = 0.;


	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:mass_sum,momentum_sum0,momentum_sum1,momentum_sum2)
	#endif
	{
		for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid();
			 ++molecule) {
			double mass = components[molecule->componentid()].m();
			mass_sum += mass;
			momentum_sum0 += mass * molecule->v(0);
			momentum_sum1 += mass * molecule->v(1);
			momentum_sum2 += mass * molecule->v(2);
		}
	}

	domainDecomp->collCommInit(4);
	domainDecomp->collCommAppendDouble(mass_sum);
	domainDecomp->collCommAppendDouble(momentum_sum0);
	domainDecomp->collCommAppendDouble(momentum_sum1);
	domainDecomp->collCommAppendDouble(momentum_sum2);
	domainDecomp->collCommAllreduceSum();
	mass_sum = domainDecomp->collCommGetDouble();
	momentum_sum0 = domainDecomp->collCommGetDouble();
	momentum_sum1 = domainDecomp->collCommGetDouble();
	momentum_sum2 = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();

	Log::global_log->info() << "momentumsum: " << momentum_sum0 << " " << momentum_sum1<< " " << momentum_sum2 << std::endl;
	Log::global_log->info() << "mass_sum: " << mass_sum << std::endl;
	double v_sub0 = momentum_sum0 / mass_sum;
	double v_sub1 = momentum_sum1 / mass_sum;
	double v_sub2 = momentum_sum2 / mass_sum;
	{
		auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		if(iter.isValid()) {
			Log::global_log->info() << "v_sub: " << v_sub0 << " " << v_sub1 << " " << v_sub2 << std::endl;
			Log::global_log->info() << "m1 v: " << std::setprecision(10) << iter->v(0) << " " << iter->v(1) << " "
									<< iter->v(2) << std::endl;
		}
	}
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{

		for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid(); ++molecule) {
			molecule->vsub(v_sub0, v_sub1, v_sub2);
		}
	}
	{
		auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		if(iter.isValid()) {
			Log::global_log->info() << "m1 v: " << iter->v(0) << " " << iter->v(1) << " " << iter->v(2)
									<< std::setprecision(5) << std::endl;
		}
	}
#ifndef NDEBUG
	//test
	momentum_sum0 = 0.;
	momentum_sum1 = 0.;
	momentum_sum2 = 0.;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:momentum_sum0,momentum_sum1,momentum_sum2)
	#endif
	{

		for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid(); ++molecule) {
			double mass = components[molecule->componentid()].m();
			momentum_sum0 += mass * molecule->v(0);
			momentum_sum1 += mass * molecule->v(1);
			momentum_sum2 += mass * molecule->v(2);
		}
	}
	Log::global_log->info() << "momentumsum: " << momentum_sum0 << " " << momentum_sum1<< " " << momentum_sum2 << std::endl;
	// Leave commented out - as there are many molecules, mass_sum is large, leading to small corrections to  the molecule velocities
	// and since molecule velocities are stored in only single precision, this is likely not
	// going to be easy to fix, across all ranges of magnitudes, which can appear here (1 to 10^9 molecules per process? to 10^13 molecules for total simulation?).

	// should we try to MPI-this? range will be even laaarger.

//	mardyn_assert(fabs(momentum_sum0)<1e-7);
//	mardyn_assert(fabs(momentum_sum1)<1e-7);
//	mardyn_assert(fabs(momentum_sum2)<1e-7);
#endif
	//printf("momentum_sum[0] from removeMomentum is %lf\n", momentum_sum[0]);
	//printf("momentum_sum[1] from removeMomentum is %lf\n", momentum_sum[1]);
	//printf("momentum_sum[2] from removeMomentum is %lf\n", momentum_sum[2]);
}


//void CubicGridGeneratorInternal::getOrientation(int base, int delta, double orientation[4]) {
//	double offset = randdouble(-delta / 2., delta / 2.) / 180. * M_PI;
//	double rad = base / 180. * M_PI;
//	double angle = rad + offset;
//
//	double cosinePart = cos(angle);
//	double sinePart = sin(angle);
//
//	double length = sqrt(cosinePart * cosinePart + 2 * (sinePart * sinePart));
//	orientation[0] = cosinePart / length;
//	orientation[1] = sinePart / length;
//	orientation[2] = sinePart / length;
//	orientation[3] = 0;
//}
//
//std::vector<double> CubicGridGeneratorInternal::getRandomVelocity(double temperature) {
//	std::vector<double> v_;
//	v_.resize(3);
//
//	// Velocity
//	for (int dim = 0; dim < 3; dim++) {
//		v_[dim] = randdouble(-0.5, 0.5);
//	}
//	double dotprod_v = 0;
//	for (unsigned int i = 0; i < v_.size(); i++) {
//		dotprod_v += v_[i] * v_[i];
//	}
//	// Velocity Correction
//	double vCorr = sqrt(3.0 * temperature / dotprod_v);
//	for (unsigned int i = 0; i < v_.size(); i++) {
//		v_[i] *= vCorr;
//	}
//
//	return v_;
//}
