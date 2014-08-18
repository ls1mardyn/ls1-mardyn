/** \file MPICheckpointWriter.h
  * \brief temporary CheckpointWriter alternative using MPI-IO
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#ifndef MPICHECKPOINTWRITER_H_
#define MPICHECKPOINTWRITER_H_

#include <string>

#include "io/OutputBase.h"

class MPICheckpointWriter : public OutputBase {
public:
	
    MPICheckpointWriter(){}
	//! @brief writes a checkpoint file that can be used to continue the simulation using MPIIO
	//!
	//! The format of the checkpointfile written by this method is
	//! Byte offset  0- 5:	string	"MarDyn"
	//! Byte offset  6-55:	string	version "YYYYMMDD\0"
	//! Byte offset 56-63:	unsigned int	gap_to_data=displacement-64
	//! Byte offset 64-(63+gap_to_data)  :	Metadata like tuple structure (ICRVQD) or offsets and cardinality to molecule sets within cubes
	//! Byte offset (64+gap_to_data)- :	data tuples
	//! 
	//! @param filename Name of the checkpointfile (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
	MPICheckpointWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~MPICheckpointWriter();
	
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
		return std::string("MPICheckpointWriter");
	}
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool	_incremental;
	bool	_appendTimestamp;
};

#endif /*MPICHECKPOINTWRITER_H_*/

