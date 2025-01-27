#include "IOHelpers.h"

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

void IOHelpers::removeMomentum(ParticleContainer* particleContainer, const std::vector<Component>& components,
							   DomainDecompBase* domainDecomp) {
	Log::global_log->info() << "Removing momentum from Input Data." << std::endl;

	double mass_sum = 0.;
	double momentum_sum[3] {0., 0., 0.};
#if defined(_OPENMP)
#pragma omp parallel reduction(+ : mass_sum, momentum_sum[:])
#endif
	{
		for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid();
			 ++molecule) {
			double mass = components[molecule->componentid()].m();
			mass_sum += mass;
			momentum_sum[0] += mass * molecule->v(0);
			momentum_sum[1] += mass * molecule->v(1);
			momentum_sum[2] += mass * molecule->v(2);
		}
	}

	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), mass_sum, momentum_sum[0], momentum_sum[1], momentum_sum[2]);
	collComm.communicate();
	std::tie(mass_sum, momentum_sum[0], momentum_sum[1], momentum_sum[2]) = collComm.get();

	Log::global_log->info() << "momentumsum prior to removal: " << momentum_sum[0] << " " << momentum_sum[1] << " "
							<< momentum_sum[2] << std::endl;
	Log::global_log->info() << "mass_sum: " << mass_sum << std::endl;
	double v_sub0 = momentum_sum[0] / mass_sum;
	double v_sub1 = momentum_sum[1] / mass_sum;
	double v_sub2 = momentum_sum[2] / mass_sum;
#if defined(_OPENMP)
#pragma omp parallel
#endif
	{
		for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid();
			 ++molecule) {
			molecule->vsub(v_sub0, v_sub1, v_sub2);
		}
	}
	// test
	momentum_sum[0] = 0.;
	momentum_sum[1] = 0.;
	momentum_sum[2] = 0.;

#if defined(_OPENMP)
#pragma omp parallel reduction(+ : momentum_sum[:])
#endif
	{
		for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid();
			 ++molecule) {
			double mass = components[molecule->componentid()].m();
			momentum_sum[0] += mass * molecule->v(0);
			momentum_sum[1] += mass * molecule->v(1);
			momentum_sum[2] += mass * molecule->v(2);
		}
	}
	Log::global_log->info() << "momentumsum after removal: " << momentum_sum[0] << " " << momentum_sum[1] << " "
							<< momentum_sum[2] << std::endl;
}

void IOHelpers::initializeVelocityAccordingToTemperature(ParticleContainer* particleContainer,
														 DomainDecompBase* /*domainDecomp*/, double temperature) {
	Random rng;
	for (auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter) {
		auto velocity = getRandomVelocityAccordingToTemperature(temperature, rng);
		for (unsigned short d = 0; d < 3; ++d) {
			iter->setv(d, velocity[d]);
		}
	}
}

std::array<double, 3> IOHelpers::getRandomVelocityAccordingToTemperature(double temperature, Random& RNG) {
	std::array<double, 3> ret{};
	double dotprod_v = 0;
	// Doing this multiple times ensures that we have equally distributed directions for the velocities!
	do {
		// Velocity
		for (int dim = 0; dim < 3; dim++) {
			ret[dim] = RNG.uniformRandInRange(-0.5f, 0.5f);
		}
		dotprod_v = ret[0] * ret[0] + ret[1] * ret[1] + ret[2] * ret[2];
	} while (dotprod_v < 0.0625);

	// Velocity Correction
	double vCorr = sqrt(3. * temperature / dotprod_v);
	for (unsigned int i = 0; i < ret.size(); i++) {
		ret[i] *= vCorr;
	}

	return ret;
}

unsigned long IOHelpers::makeParticleIdsUniqueAndGetTotalNumParticles(ParticleContainer* particleContainer,
																	  DomainDecompBase* domainDecomp) {
	unsigned long localNumParticles = particleContainer->getNumberOfParticles();

	auto collCommScan = makeCollCommObjScanAdd(domainDecomp->getCommunicator(), localNumParticles);
	collCommScan.communicate();
	auto [idOffset] = collCommScan.get();
	idOffset -= localNumParticles;

	// fix ID's to be unique:
#if defined(_OPENMP)
#pragma omp parallel
#endif
	{
		for (auto mol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
			mol->setid(mol->getID() + idOffset);
		}
	}

	auto collComm = makeCollCommObjAllreduceAdd<2>(domainDecomp->getCommunicator(), localNumParticles);
	collComm.communicate();
	auto [globalNumParticles] = collComm.get();

	return globalNumParticles;
}
