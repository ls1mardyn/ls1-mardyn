/*
 * MPI_IOWriter.h
 *
 *  Created on: Aug 18, 2014
 *      Author: andal
 */

#ifndef MPI_IOCHECKPOINTWRITER_H_
#define MPI_IOCHECKPOINTWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

class MPI_IOCheckpointWriter  : public PluginBase{
public:

	MPI_IOCheckpointWriter(){}

	//! @brief writes a checkpoint file that can be used to continue the simulation
	//!
	//! The format of the checkpointfile written by this method is the same as the format
	//! of the input file.
	//! @param filename Name of the checkpointfile (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
	MPI_IOCheckpointWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~MPI_IOCheckpointWriter();

	void readXML(XMLfileUnits& xmlconfig);

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);
	void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep, list<ChemicalPotential> * /*lmu*/, map<unsigned, CavityEnsemble> * /*mcav*/
    );
	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);

	void handle_error(int i);

	std::string getPluginName() {
		return std::string("MPI_IOCheckpointWriter");
	}
	static PluginBase* createInstance() { return new MPI_IOCheckpointWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool	_incremental;
	bool	_appendTimestamp;

};

#endif /* MPI_IOCHECKPOINTWRITER_H_ */
