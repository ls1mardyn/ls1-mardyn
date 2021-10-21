#include "MamicoCoupling.h"
#include "Domain.h"
#include "coupling/interface/impl/ls1/LS1MamicoCouplingSwitch.h"

void MamicoCoupling::readXML(XMLfileUnits& xmlconfig)
{	
    return;
}

void MamicoCoupling::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
    if(_macroscopicCellService==NULL) 
    {
        _macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<ParticleCell,3>*)
                coupling::interface::MamicoInterfaceProvider<ParticleCell,3>::getInstance().getMacroscopicCellService();
	}
    //since using mamico thermostat, switch off ls1 thermostat
    domain->thermostatOff();
	//code to print to log that plugin is initialised
}

void MamicoCoupling::beforeEventNewTimestep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{
	
}

void MamicoCoupling::beforeForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)        
{
	if(coupling::interface::LS1MamicoCouplingSwitch::getInstance().getCouplingState())
    {
        if(_macroscopicCellService==NULL) 
        {
            _macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<ParticleCell,3>*)
                    coupling::interface::MamicoInterfaceProvider<ParticleCell,3>::getInstance().getMacroscopicCellService();
        }
        _macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();
        _macroscopicCellService->distributeMass(simstep);
        _macroscopicCellService->applyTemperatureToMolecules(simstep);
    }
}

void MamicoCoupling::afterForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(coupling::interface::LS1MamicoCouplingSwitch::getInstance().getCouplingState())
    {
        if(_macroscopicCellService==NULL) 
        {
            _macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<ParticleCell,3>*)
                    coupling::interface::MamicoInterfaceProvider<ParticleCell,3>::getInstance().getMacroscopicCellService();
        }
        _macroscopicCellService->distributeMomentum(simstep);
        _macroscopicCellService->applyBoundaryForce(simstep);
    }
}   

void MamicoCoupling::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep)
{
	
}

void MamicoCoupling::finish(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
	
}