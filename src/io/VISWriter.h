#ifndef VISWRITER_H_
#define VISWRITER_H_

#include "plugins/PluginBase.h"
#include "Domain.h"
#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase;
class Domain;

class VISWriter : public PluginBase {
public:
    VISWriter(){}
	//! @brief Writes out a file (using *.vis-format) containing coordinates + orientation (using quaternions)
	//! of each molecule for several timesteps.
	//!
	//! Depending on write frequency (for example: every timestep, or every 10th, 100th, 1000th ...) number of frames
	//! can be controlled. The *.vis-file can be visualized by visualization software like:
	//!   - MolCloud (visit: http://www.visus.uni-stuttgart.de/index.php?id=995)
	//!   - MegaMol  (visit: https://svn.vis.uni-stuttgart.de/trac/megamol/)
	//!
	//! @param filename Name of the *.vis-file (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
	VISWriter(unsigned long writeFrequency, std::string outputPrefix);
	~VISWriter();

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
		return std::string("VISWriter");
	}
	static PluginBase* createInstance() { return new VISWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	bool _wroteVIS;
};

#endif /* VISWRITER_H_ */
