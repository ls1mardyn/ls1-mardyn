#ifndef SRC_PLUGINS_MAMICOCOUPLING_H_
#define SRC_PLUGINS_MAMICOCOUPLING_H_

#include "PluginBase.h"

#ifdef MAMICO_COUPLING
#include "coupling/services/MacroscopicCellService.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#endif

#include "particleContainer/LinkedCells.h"


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
    coupling::services::MacroscopicCellServiceImpl<ParticleCell,3>* _macroscopicCellService;
#endif
};

#endif /* SRC_PLUGINS_MAMICOCOUPLING_H_ */