#ifndef XYZWRITER_H_
#define XYZWRITER_H_

#include <string>

#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"

//! @brief Writes out a file (using *.xyz-format) containing coordinates of each molecule.
//!
//! Depending on write frequency (for example: every timestep, or every 10th, 100th, 1000th ...) number of frames
//! can be controlled. The *.xyz-file can be visualized by visualization software like vmd.
//! (for detail information visit: http://www.ks.uiuc.edu/Research/vmd/)
//!
//! @param filename Name of the *.xyz-file (including path)
//! @param particleContainer The molecules that have to be written to the file are stored here
//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
//!                     Methods to achieve this are available in domainDecomp
//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
class XyzWriter : public OutputBase {
public:
    XyzWriter(){}
	XyzWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~XyzWriter();
	//! @todo comment
	
	void readXML(XMLfileUnits& xmlconfig);

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	//! @todo comment
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu
	);
	//! @todo comment
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("XyzWriter");
	}
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	bool _incremental;
};

#endif /* XYZWRITER_H_ */
