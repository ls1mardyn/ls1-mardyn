#ifdef MAMICO_COUPLING

#include "MamicoCoupling.h"
#include "Domain.h"

void MamicoCoupling::readXML(XMLfileUnits &xmlconfig) { return; }

void MamicoCoupling::init(ParticleContainer *particleContainer,
                          DomainDecompBase *domainDecomp, Domain *domain) {

  Log::global_log->info() << "MaMiCo coupling plugin initialized" << std::endl;

}

void MamicoCoupling::beforeEventNewTimestep(
    ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
    unsigned long simstep) {}

void MamicoCoupling::beforeForces(ParticleContainer *particleContainer,
                                  DomainDecompBase *domainDecomp,
                                  unsigned long simstep) {
  if (_couplingEnabled) {
    // This object should be set by MaMiCo after the plugins are created in the
    // simulation readxml file Even though this method is called before the
    // object is set, at this point the coupling switch is always off
    _couplingCellService->processInnerCouplingCellAfterMDTimestep();
    _couplingCellService->distributeMass(simstep);
    _couplingCellService->applyTemperatureToMolecules(simstep);
#ifndef MARDYN_AUTOPAS
    particleContainer->deleteOuterParticles();
#endif
  }
}

void MamicoCoupling::afterForces(ParticleContainer *particleContainer,
                                 DomainDecompBase *domainDecomp,
                                 unsigned long simstep) {
  if (_couplingEnabled) {
    _couplingCellService->distributeMomentum(simstep);
    _couplingCellService->applyBoundaryForce(simstep);
  }
}

void MamicoCoupling::endStep(ParticleContainer *particleContainer,
                             DomainDecompBase *domainDecomp, Domain *domain,
                             unsigned long simstep) {}

void MamicoCoupling::finish(ParticleContainer *particleContainer,
                            DomainDecompBase *domainDecomp, Domain *domain) {}

#endif
