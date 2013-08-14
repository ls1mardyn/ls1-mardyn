#ifndef POVWRITER_H_
#define POVWRITER_H_

#include <string>

#include "io/OutputBase.h"


class PovWriter : public OutputBase {
public:
    PovWriter(){}
	//! @brief writes a POVray file of the current state of the simluation
	//!
	//! The file can be used to visualize with POVray software (for detail information visit: www.povray.org)
	//!
	//! @param filename Name of the POV file (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
	PovWriter(unsigned long writeFrequency, std::string filename, bool incremental);
	~PovWriter();

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
		return std::string("PovWriter");
	}
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool  _incremental;
	bool  _appendTimestamp;
};

#endif /* POVWRITER_H_ */
