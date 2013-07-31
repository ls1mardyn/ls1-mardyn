#ifndef DECOMPWRITER_H_
#define DECOMPWRITER_H_

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <fstream>
#include <sstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
//class Molecule;
//class DecompWriter; 

//! @brief writes out information about decomposition of the simulation domain.
//!
//! @param filename Name of the written file (including path)
//! @param particleContainer The molecules that are contained in the simulation domain
//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
//!                     Methods to achieve this are available in domainDecomp
//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
class DecompWriter : public OutputBase {
public:
	DecompWriter(unsigned long writeFrequency, std::string mode, std::string filename, unsigned long numberOfTimesteps, bool incremental);
	~DecompWriter();
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
	unsigned long _numberOfTimesteps;
	unsigned long _writeFrequency;
	std::string _mode;
	bool _filenameisdate;
	bool _incremental;
	std::string _filename;
};

#endif /*DECOMPWRITER_H_*/
