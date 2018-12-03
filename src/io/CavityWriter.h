#ifndef CAVITYWRITER_H_
#define CAVITYWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

class CavityEnsemble;

class CavityWriter : public PluginBase {
public:
    CavityWriter(){}
	CavityWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~CavityWriter(){};
	//! @todo comment
	
	void readXML(XMLfileUnits& xmlconfig);

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain);

	/** @brief Method will be called first thing in a new timestep. */
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) final;

	/** @brief Method afterForces will be called after forcefields have been applied
     *  no sitewise Forces can be applied here
     */
	void afterForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) final;

	//! @todo comment
	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    );
	//! @todo comment
	void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain);
	
	std::string getPluginName() {
		return std::string("CavityWriter");
	}
	static PluginBase* createInstance() { return new CavityWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	bool _incremental;
	std::map<unsigned, CavityEnsemble> _mcav;
};

#endif /* CAVITYWRITER_H_ */
