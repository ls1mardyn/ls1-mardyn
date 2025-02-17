#include "IOHelpers.h"

#include "Simulation.h"
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

	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), mass_sum, momentum_sum[0], momentum_sum[1], momentum_sum[2]);
	collComm.communicate();
	std::tie(mass_sum, momentum_sum[0], momentum_sum[1], momentum_sum[2]) = collComm.get();

	Log::global_log->info() << "momentumsum prior to removal: " << momentum_sum[0] << " " << momentum_sum[1] << " "
							<< momentum_sum[2] << std::endl;
	Log::global_log->info() << "mass_sum: " << mass_sum << std::endl;
	const double v_sub0 = momentum_sum[0] / mass_sum;
	const double v_sub1 = momentum_sum[1] / mass_sum;
	const double v_sub2 = momentum_sum[2] / mass_sum;
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

	// After subtraction, the temperature was effectively changed
	// Therefore, the velocities have to be scaled again to match target temperature
	const double temperature_target = global_simulation->getEnsemble()->T();
	double ekin_sum = 0.0;
	unsigned long numDOFs = 0;

#if defined(_OPENMP)
#pragma omp parallel reduction(+ : ekin_sum, numDOFs)
#endif
	{
		for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid();
			 ++molecule) {
			ekin_sum += molecule->U_kin();
			numDOFs += (3+molecule->component()->getRotationalDegreesOfFreedom());
		}
	}

	domainDecomp->collCommInit(2);
	domainDecomp->collCommAppendDouble(ekin_sum);
	domainDecomp->collCommAppendUnsLong(numDOFs);
	domainDecomp->collCommAllreduceSum();
	ekin_sum = domainDecomp->collCommGetDouble();
	numDOFs = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();

	const double temperature_current = 2*ekin_sum/numDOFs;
	const double scaleFactor = std::sqrt(temperature_target/temperature_current);
	Log::global_log->info() << "current temperature: " << temperature_current
							<< "scale velocities by: " << scaleFactor << std::endl;
#if defined(_OPENMP)
#pragma omp parallel
#endif
	{
		for (auto mol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid();
			 ++mol) {
			for (int d = 0; d < 3; d++) {
				// Scale velocities to match target temperature
				mol->setv(d,mol->v(d)*scaleFactor);
			}
		}
	}
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
