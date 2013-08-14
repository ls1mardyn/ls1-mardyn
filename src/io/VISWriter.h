#ifndef VISWRITER_H_
#define VISWRITER_H_

#include "io/OutputBase.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
class ChemicalPotential;

class VISWriter : public OutputBase {
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
	VISWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~VISWriter();

	void readXML(XMLfileUnits& xmlconfig);

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu
	);
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("VISWriter");
	}
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _incremental;
	bool _appendTimestamp;
	bool _wroteVIS;
};

#endif /* VISWRITER_H_ */
