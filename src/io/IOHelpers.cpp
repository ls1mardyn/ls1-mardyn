#include "IOHelpers.h"

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/generator/EqualVelocityAssigner.h"

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

	domainDecomp->collCommInit(4);
	domainDecomp->collCommAppendDouble(mass_sum);
	domainDecomp->collCommAppendDouble(momentum_sum[0]);
	domainDecomp->collCommAppendDouble(momentum_sum[1]);
	domainDecomp->collCommAppendDouble(momentum_sum[2]);
	domainDecomp->collCommAllreduceSum();
	mass_sum = domainDecomp->collCommGetDouble();
	momentum_sum[0] = domainDecomp->collCommGetDouble();
	momentum_sum[1] = domainDecomp->collCommGetDouble();
	momentum_sum[2] = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();

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
	EqualVelocityAssigner eqVeloAssigner(temperature);
	for (auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter) {
		eqVeloAssigner.assignVelocity(&(*iter));
	}
}

unsigned long IOHelpers::makeParticleIdsUniqueAndGetTotalNumParticles(ParticleContainer* particleContainer,
																	  DomainDecompBase* domainDecomp) {
	unsigned long localNumParticles = particleContainer->getNumberOfParticles();

	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(localNumParticles);
	domainDecomp->collCommScanSum();
	unsigned long idOffset = domainDecomp->collCommGetUnsLong() - localNumParticles;
	domainDecomp->collCommFinalize();
	// fix ID's to be unique:
#if defined(_OPENMP)
#pragma omp parallel
#endif
	{
		for (auto mol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
			mol->setid(mol->getID() + idOffset);
		}
	}

	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(localNumParticles);
	domainDecomp->collCommAllreduceSum();
	unsigned long globalNumParticles = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
	return globalNumParticles;
}
