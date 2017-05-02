

#ifndef BINARYCHECKPOINTWRITER_H_
#define BINARYCHECKPOINTWRITER_H_

#include <string>

#include "io/OutputBase.h"

class BinaryCheckpointWriter: public OutputBase {
public:

	BinaryCheckpointWriter();
	//! @brief writes a checkpoint file that can be used to continue the simulation
	//!
	//! The format of the checkpointfile written by this method is the same as the format
	//! of the input file.
	//! @param filename Name of the checkpointfile (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
	BinaryCheckpointWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~BinaryCheckpointWriter();

	void readXML(XMLfileUnits& xmlconfig);

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void doOutput(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, list<ChemicalPotential>* /*lmu*/, map<unsigned, CavityEnsemble>* /*mcav*/
	);
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

	std::string getPluginName() {
		return std::string("CheckpointWriter");
	}
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool	_incremental;
	bool	_appendTimestamp;
};

#endif /* BINARYCHECKPOINTWRITER_H_ */
