#ifndef SRC_PLUGINS_MAMICOCOUPLING_H_
#define SRC_PLUGINS_MAMICOCOUPLING_H_

#include "PluginBase.h"
#ifdef MAMICO_COUPLING
#include <climits> //for INT_MAX in indexconversion
#include <coupling/services/MacroscopicCellService.h>
#include <coupling/interface/MamicoInterfaceProvider.h>
#include <coupling/interface/impl/ls1/LS1RegionWrapper.h>
//#include <MacroscopicCellService.h>
//#include <MamicoInterfaceProvider.h>
#endif

//#include "particleContainer/LinkedCells.h"


class MamicoCoupling: public PluginBase {

public:
	MamicoCoupling() {}
	virtual ~MamicoCoupling() {}

	void init(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

    void readXML(XMLfileUnits& xmlconfig);

	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	);

    void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep);

    void finish(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp, Domain* domain);

    std::string getPluginName() {
    	return std::string("MamicoCoupling");
    }

	static PluginBase* createInstance() { return new MamicoCoupling(); }

private:
#ifdef MAMICO_COUPLING
    coupling::services::MacroscopicCellServiceImpl<ls1::LS1RegionWrapper,3>* _macroscopicCellService = nullptr;
#endif
};

#endif /* SRC_PLUGINS_MAMICOCOUPLING_H_ */