#include "MamicoCoupling.h"
#include "Domain.h"

#ifdef MAMICO_COUPLING
#include <coupling/interface/impl/ls1/LS1MamicoCouplingSwitch.h>
//#include <LS1MamicoCouplingSwitch.h>
#endif

void MamicoCoupling::readXML(XMLfileUnits& xmlconfig)
{	
    return;
}

void MamicoCoupling::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
    #ifdef MAMICO_COUPLING
    if(_macroscopicCellService==NULL) 
    {
        _macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>*)
                coupling::interface::MamicoInterfaceProvider<ls1::LS1RegionWrapper,3>::getInstance().getMacroscopicCellService();
	}
    //since using mamico thermostat, switch off ls1 thermostat 
    domain->thermostatOff();
	//code to print to log that plugin is initialised
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
	if(coupling::interface::LS1MamicoCouplingSwitch::getInstance().getCouplingState())
    {
        if(_macroscopicCellService==NULL) 
        {
            _macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>*)
                    coupling::interface::MamicoInterfaceProvider<ls1::LS1RegionWrapper,3>::getInstance().getMacroscopicCellService();
        }
        _macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();
        _macroscopicCellService->distributeMass(simstep);
        _macroscopicCellService->applyTemperatureToMolecules(simstep);
    }
    #endif
}

void MamicoCoupling::afterForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{
    #ifdef MAMICO_COUPLING
    if(coupling::interface::LS1MamicoCouplingSwitch::getInstance().getCouplingState())
    {
        if(_macroscopicCellService==NULL) 
        {
            _macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>*)
                    coupling::interface::MamicoInterfaceProvider<ls1::LS1RegionWrapper,3>::getInstance().getMacroscopicCellService();
        }
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