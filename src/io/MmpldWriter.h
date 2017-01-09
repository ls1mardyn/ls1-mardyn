#ifndef MMPLDWRITER_H_
#define MMPLDWRITER_H_

#include <string>

#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"

class MmpldWriter : public OutputBase{
  public:
    MmpldWriter(){};
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
    MmpldWriter(unsigned long writeFrequency, std::string outputPrefix);
	~MmpldWriter();

	void readXML(XMLfileUnits& xmlconfig);

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu,
			std::map<unsigned, CavityEnsemble>* mcav
	);
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("MmpldWriter");
	}
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	std::string _timestampString;
	uint32_t _numSeekEntries;
	uint32_t _frameCount;
#ifdef ENABLE_MPI
	uint64_t *_seekTable;
#endif
};

#endif /* MMPLDWRITER_H_ */
