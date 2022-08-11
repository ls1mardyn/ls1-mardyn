#include "MamicoCoupling.h"
#include "Domain.h"

#ifdef MAMICO_COUPLING
#include <coupling/interface/MamicoInterfaceProvider.h>
#include <coupling/interface/impl/ls1/LS1MamicoCouplingSwitch.h>
#endif

void MamicoCoupling::readXML(XMLfileUnits& xmlconfig)
{	
	return;
}

/** @brief Ensures _macroscopicCellService is initialized, and switches the ls1 thermostat off.
 * 
 * MaMiCo uses its own thermostat, which will conflict with the ls1 thermostat. */
void MamicoCoupling::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
	#ifdef MAMICO_COUPLING
	//since using mamico thermostat, switch off ls1 thermostat 
	domain->thermostatOff();
	//code to print to log that plugin is initialised
	global_log->info() << "MaMiCo coupling plugin initialized" << std::endl;
	#endif
}

void MamicoCoupling::beforeEventNewTimestep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{

}

/** @brief Takes coupling steps such as particle insertion, to make sure they are accounted for before forces are calculated.
 * 
 * Following steps are taken, if coupling is switched on:
 * - Iterate over cells to average values like momentum and mass, to pass to macroscopic solvers
 * - Distribute incoming mass from macroscopic solver by inserting perticles (if enabled)
 * - Run the MaMiCo thermostat cell by cell
 * */
#pragma GCC push_options
#pragma GCC optimize("O0")
void MamicoCoupling::beforeForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)        
{
	#ifdef MAMICO_COUPLING
	if(coupling::interface::LS1MamicoCouplingSwitch::getInstance().getCouplingState()) //only perform coupling steps if coupling is enabled
	{
		// This object should be set by MaMiCo after the plugins are created in the simulation readxml file
		// Even though this method is called before the object is set, at this point the coupling switch is always off
		_macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();
		_macroscopicCellService->distributeMass(simstep);
		_macroscopicCellService->applyTemperatureToMolecules(simstep);
    }
	#endif
}

/** @brief Performs adjustments after force calculation 
 * 
 * Following steps are taken, if coupling is switched on:
 * - Distribute incoming momentum among affected cells
 * - Apply boundary force to molecules near microscopic domain boundary
 * */
void MamicoCoupling::afterForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{
	#ifdef MAMICO_COUPLING
	if(coupling::interface::LS1MamicoCouplingSwitch::getInstance().getCouplingState())
	{
		_macroscopicCellService->distributeMomentum(simstep);
		_macroscopicCellService->applyBoundaryForce(simstep);
	}
	#endif
}
#pragma GCC pop_options

void MamicoCoupling::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep)
{

}

void MamicoCoupling::finish(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
	
}