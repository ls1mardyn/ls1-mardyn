#ifndef MMSPDBINWRITER_H_
#define MMSPDBINWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

class MmspdBinWriter : public PluginBase{
  public:
    MmspdBinWriter(){};
	//! @brief: writes a mmspd file used by MegaMol
	//!
	//! Depending on write frequency (for example: every timestep, or every 10th, 100th, 1000th ...) number of frames
	//! can be controlled. The *.mmspd-file can be visualized by visualization software like MegaMol.
	//! (for detail information visit: https://svn.vis.uni-stuttgart.de/trac/megamol/)
	//!
	//! @param filename Name of the *.mmspd-file (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
    MmspdBinWriter(unsigned long writeFrequency, std::string outputPrefix);
	~MmspdBinWriter();

	void readXML(XMLfileUnits& xmlconfig);

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);
	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    );
	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() {
		return std::string("MmspdBinWriter");
	}
	static PluginBase* createInstance() { return new MmspdBinWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
};

#endif /* MMSPDBINWRITER_H_ */
