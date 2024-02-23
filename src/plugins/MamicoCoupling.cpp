#include "MamicoCoupling.h"
#include "Domain.h"

void MamicoCoupling::readXML(XMLfileUnits& xmlconfig)
{
	return;
}

void MamicoCoupling::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
#ifdef MAMICO_COUPLING
	//code to print to log that plugin is initialised
	Log::global_log->info() << "MaMiCo coupling plugin initialized" << std::endl;
#endif
}

void MamicoCoupling::beforeEventNewTimestep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{

}

void MamicoCoupling::beforeForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{
#ifdef MAMICO_COUPLING
	if(_couplingEnabled)
	{
		// This object should be set by MaMiCo after the plugins are created in the simulation readxml file
		// Even though this method is called before the object is set, at this point the coupling switch is always off
		_macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();
		_macroscopicCellService->distributeMass(simstep);
		_macroscopicCellService->applyTemperatureToMolecules(simstep);
#ifndef MARDYN_AUTOPAS
		particleContainer->deleteOuterParticles();
#endif
	}
#endif
}

void MamicoCoupling::afterForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{
#ifdef MAMICO_COUPLING
	if(_couplingEnabled)
	{
		_macroscopicCellService->distributeMomentum(simstep);
		_macroscopicCellService->applyBoundaryForce(simstep);
	}
#endif
}

void MamicoCoupling::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep)
{

}

void MamicoCoupling::finish(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{

}
