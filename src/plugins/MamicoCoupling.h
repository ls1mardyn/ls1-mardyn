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

	void init(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) override;

	void readXML(XMLfileUnits& xmlconfig) override;

	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	void beforeForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

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
/** @brief sets the macroscopicCellService object that controls the inner coupling logic.
 * 
 * MaMiCo extracts the MamicoCoupling plugin object from the simulation object after initialization and uses this function to set the macroscopicCellCervice.
 * The code for this object can be found in https://github.com/HSU-HPC/MaMiCo/blob/master/coupling/services/MacroscopicCellService.cpph
 * */
	void setMamicoMacroscopicCellService(coupling::services::MacroscopicCellService<3>* macroscopicCellService){
		_macroscopicCellService = static_cast<coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>*>
			(macroscopicCellService);
	}

	void switchOnCoupling(){ _couplingEnabled = true; }
	void switchOffCoupling(){ _couplingEnabled = false; }
	bool getCouplingState() { return _couplingEnabled;}
#endif

private:
#ifdef MAMICO_COUPLING
	coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>* _macroscopicCellService;
	bool _couplingEnabled = false;
#endif
};
