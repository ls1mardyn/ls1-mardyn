/** \file MPICheckpointWriter.h
  * \brief temporary CheckpointWriter alternative using MPI-IO
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#ifndef MPICHECKPOINTWRITER_H_
#define MPICHECKPOINTWRITER_H_

#include <string>

#include "plugins/PluginBase.h"
#ifdef ENABLE_MPI
#include "utils/MPI_Info_object.h"
#endif

class MPICheckpointWriter : public PluginBase {
public:

    MPICheckpointWriter(){}
	//! @brief writes a checkpoint file that can be used to continue the simulation using MPIIO
	//!
	//! The format of the checkpointfile written by this method is
	//! Byte offset  0 -  5,  6:	string	"MarDyn"
	//! Byte offset  6 - 19, 14:	string	version: "20150122trunk\0"
	//! Byte offset 20 - 51, 32:		reserved for version string extension
	//! Byte offset 52 - 55,  4:	int	store 0x0a0b0c0d=168496141 to check endianess
	//! Byte offset 56 - 63,  8:	unsigned int	gap_to_data=data_displ-64=18+numBB*64
	//! Byte offset 64 - 70,  7:	string	tuple structure "ICRVQD\0"
	//! Byte offset 71 - 73,  3:	string	"BB\0"
	//! Byte offset 74 - 81,  8:	unsigned long	number of bounding boxes
	//! Byte offset 82 - (81+numBB*(6*8+2*8)), numBB*64:	numBB*(6*double+2*unsigned long)	bounding boxes
	//! Byte offset (64+gap_to_data) - :	data tuples
	//!
	//! @param writeFrequency	Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
	//! @param outputPrefix	path and prefix for file name used
	//! @param incremental	add simulation step to file name
	//! @param datarep	data representation (e.g. "external32")
	MPICheckpointWriter(unsigned long writeFrequency
	                   , std::string outputPrefix, bool incremental=true
	                   , std::string datarep=std::string(""));

	/** @brief Read in XML configuration for MPICheckpointWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
		<outputplugin name="MPICheckpointWriter">
			<writefrequency>INTEGER</writefrequency> <!-- Frequency in which the output is written; Default: 1 -->
			<outputprefix>STRING</outputprefix> <!-- Prefix of the output file; Default: "mardyn" -->
			<incremental>BOOL</incremental> <!-- Checkpoint files will get individual numbers; Default: false -->
			<appendTimestamp>BOOL</appendTimestamp> <!-- Append timestamp to checkpoint files; Default: false -->
			<datarep>STRING</datarep> <!-- MPI I/O output representation to use, valid values are "native", "internal", "external32"; Default: "" -->
			<mpi_info><!-- see MPI_Info_object class documentation --></mpi_info> <!-- MPI infos to be used for MPI file I/O writing the checkpoint files -->
		</outputplugin>
	   \endcode
	 */
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
		return std::string("MPICheckpointWriter");
	}
	static PluginBase* createInstance() { return new MPICheckpointWriter(); }
private:
	static const char _magicVersion[56];
	static const int _endiannesstest;

	std::string	_outputPrefix;
	unsigned long	_writeFrequency;
	bool	_incremental;
	bool	_appendTimestamp;
	std::string	_datarep;
	bool	_measureTime;
#ifdef ENABLE_MPI
	MPI_Info_object	_mpiinfo;
#endif
	unsigned long	_particlesbuffersize = 0;
};

#endif /*MPICHECKPOINTWRITER_H_*/

