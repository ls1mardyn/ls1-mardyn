#ifndef CAVITYWRITER_H_
#define CAVITYWRITER_H_

#include <string>

#include "ensemble/CavityEnsemble.h"
#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"

class CavityWriter : public OutputBase {
public:
    CavityWriter(){}
	CavityWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~CavityWriter();
	//! @todo comment
	
	void readXML(XMLfileUnits& xmlconfig);

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	//! @todo comment
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu, std::map<unsigned, CavityEnsemble>* mcav
	);
	//! @todo comment
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	static std::string getPluginName() {
		return std::string("CavityWriter");
	}
	static OutputBase* createInstance() { return new CavityWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	bool _incremental;
};

#endif /* CAVITYWRITER_H_ */
