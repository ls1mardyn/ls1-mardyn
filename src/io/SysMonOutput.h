#ifndef SYSMONOUTPUT_H_
#define SYSMONOUTPUT_H_

#include "plugins/PluginBase.h"

#include <fstream>

class SysMonOutput : public PluginBase {
public:
	SysMonOutput();
	~SysMonOutput(){}

	void readXML(XMLfileUnits& xmlconfig);

	//! @todo comment
	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain);
	//! @todo comment

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                 unsigned long simstep);

	//! @todo comment
	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() {
		return std::string("SysMonOutput");
	}
	static PluginBase* createInstance() { return new SysMonOutput(); }

private:
	//! prefix for the names of all output files
	std::ofstream _resultStream;
	unsigned long _writeFrequency;
};

#endif /* SYSMONOUTPUT_H_ */
