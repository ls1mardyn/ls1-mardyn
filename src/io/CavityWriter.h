#ifndef CAVITYWRITER_H_
#define CAVITYWRITER_H_

#include <string>

#include "ensemble/CavityEnsemble.h"
#include "ensemble/GrandCanonical.h"
#include "utils/PluginBase.h"

class CavityWriter : public PluginBase {
public:
    CavityWriter(){}
	CavityWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~CavityWriter();
	//! @todo comment
	
	void readXML(XMLfileUnits& xmlconfig);

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain);
	//! @todo comment
	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep, std::list<ChemicalPotential> *lmu, std::map<unsigned, CavityEnsemble> *mcav
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
};

#endif /* CAVITYWRITER_H_ */
