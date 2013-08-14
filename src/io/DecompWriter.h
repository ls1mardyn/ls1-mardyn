#ifndef DECOMPWRITER_H_
#define DECOMPWRITER_H_

#include <string>

#include "io/OutputBase.h"


//! @brief writes out information about decomposition of the simulation domain.
//!
//! @param filename Name of the written file (including path)
//! @param particleContainer The molecules that are contained in the simulation domain
//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
//!                     Methods to achieve this are available in domainDecomp
//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
class DecompWriter : public OutputBase {
public:
    DecompWriter(){}
	DecompWriter(unsigned long writeFrequency, std::string mode, std::string outputPrefix, bool incremental);
	~DecompWriter();

	void readXML(XMLfileUnits& xmlconfig);

	//! @todo comment
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
		return std::string("DecompWriter");
	}
private:
	unsigned long _writeFrequency;
	std::string _mode;
	bool _appendTimestamp;
	bool _incremental;
	std::string _outputPrefix;
};

#endif /*DECOMPWRITER_H_*/
