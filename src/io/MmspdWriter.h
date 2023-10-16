#ifndef MMSPDWRITER_H_
#define MMSPDWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

class MmspdWriter : public PluginBase{
  public:
    MmspdWriter(){};
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
    MmspdWriter(unsigned long writeFrequency, std::string outputPrefix);
	~MmspdWriter();

	void readXML(XMLfileUnits& xmlconfig);

    void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);
	//! @todo comment
    void endStep(ParticleContainer *particleContainer,
                 DomainDecompBase *domainDecomp, Domain *domain,
                 unsigned long simstep
    );
	//! @todo comment
    void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);
	std::string getPluginName() {
		return std::string("MmspdWriter");
	}
	static PluginBase* createInstance() { return new MmspdWriter(); }
  private:
      std::string _outputPrefix;
	  std::string _filename;
      unsigned long _writeFrequency;
      bool _appendTimestamp;
};

#endif /* MMSPDWRITER_H_ */
