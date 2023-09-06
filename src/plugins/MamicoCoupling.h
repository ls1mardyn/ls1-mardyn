#pragma once

#include "PluginBase.h"
#ifdef MAMICO_COUPLING
#include <coupling/services/MacroscopicCellService.h>
#include <coupling/interface/impl/ls1/LS1RegionWrapper.h>
#endif

/** @brief Allows execution of MaMiCo code to enable coupling with MaMiCo.
 * 
 * When MaMiCo is coupled with any simulation software, it needs to control the simulation from both the outside and inside.
 * From the outside, MaMiCo needs to be able to start and stop the simulation, and keep unique simulations distinct.
 * From the inside, MaMiCo code needs to be executed at specific points within the same simulation step. This plugin enables the "inside" behaviour.
 * 
 * MaMiCo coupling requires the MAMICO_COUPLING flag to be set, and the MAMICO_SRD_DIR flag to point to MaMiCo source files. With these enabled, ls1 compiles as a library.
 * To make sure the program compiles and works with tests, the relevant MaMiCo portions for the code are put in #ifdef regions.
 * 
 * */

class MamicoCoupling: public PluginBase {

public:
	MamicoCoupling() {}
	virtual ~MamicoCoupling() {}

	/** @brief Prints to log that mamico coupling is inited
	 * 
	 * */
	void init(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) override;


	void readXML(XMLfileUnits& xmlconfig) override;

	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	/** @brief Takes coupling steps such as particle insertion, to make sure they are accounted for before forces are calculated.
	 * 
	 * Following steps are taken, if coupling is switched on:
	 * - Iterate over cells to average values like momentum and mass, to pass to macroscopic solvers
	 * - Distribute incoming mass from macroscopic solver by inserting perticles (if enabled)
	 * - Run the MaMiCo thermostat cell by cell
	 * 
	 * The distributeMass method calls the updateParticleContainerAndDecomposition() function at the end, so we end up with a container with halo particles present.
	 * Hence a manual halo clearance is done to restore the state of the container.
	 * */
	void beforeForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	/** @brief Performs adjustments after force calculation 
	 * 
	 * Following steps are taken, if coupling is switched on:
	 * - Distribute incoming momentum among affected cells
	 * - Apply boundary force to molecules near microscopic domain boundary
	 * */	
	void afterForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	void endStep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep) override;

	void finish(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) override;

	std::string getPluginName() override {
		return std::string("MamicoCoupling");
	}

	static PluginBase* createInstance() { return new MamicoCoupling(); }

#ifdef MAMICO_COUPLING
/** @brief Sets the macroscopicCellService object that controls the inner coupling logic.
 * 
 * MaMiCo extracts the MamicoCoupling plugin object from the simulation object after initialization and uses this function to set the macroscopicCellCervice.
 * The code for this object can be found in https://github.com/HSU-HPC/MaMiCo/blob/master/coupling/services/MacroscopicCellService.cpph
 * */
	void setMamicoMacroscopicCellService(coupling::services::MacroscopicCellService<3>* macroscopicCellService){
		_macroscopicCellService = static_cast<coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>*>
			(macroscopicCellService);
	}
#endif

/** @brief Enables coupling logic, allowing coupling steps to run while simulation is running.
 * 
 * Set from within MaMiCo, check https://github.com/HSU-HPC/MaMiCo/blob/master/coupling/interface/MDSimulationFactory.h, class LS1MDSimulation
 * */
	void switchOnCoupling(){ _couplingEnabled = true; }

/** @brief Disables coupling logic, not allowing coupling steps to run. Typically done when equilibrating the simulation initially.
 * 
 * Set from within MaMiCo, check https://github.com/HSU-HPC/MaMiCo/blob/master/coupling/interface/MDSimulationFactory.h, class LS1MDSimulation
 * */
	void switchOffCoupling(){ _couplingEnabled = false; }

/** @brief Getter method for the coupling state.
 * 
 * Not used currently, but may be used in future
 * */
	bool getCouplingState() { return _couplingEnabled;}


private:
#ifdef MAMICO_COUPLING
	coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>* _macroscopicCellService;
#endif
	bool _couplingEnabled = false;
};
