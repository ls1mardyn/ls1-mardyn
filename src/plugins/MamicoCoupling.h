#pragma once

#include "PluginBase.h"
#ifdef MAMICO_COUPLING
#include <coupling/services/MacroscopicCellService.h>
#include <coupling/interface/impl/ls1/LS1RegionWrapper.h>
#endif

/** @brief Allows execution of MaMiCo code to enable coupling with MaMiCo.
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

private:
#ifdef MAMICO_COUPLING
	coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>* _macroscopicCellService = nullptr;
#endif
};
