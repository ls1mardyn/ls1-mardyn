#ifndef CHECKPOINTWRITER_H_
#define CHECKPOINTWRITER_H_

#include <string>

#include "io/OutputBase.h"


class CheckpointWriter : public OutputBase {
public:
	
    CheckpointWriter(){}
	//! @brief writes a checkpoint file that can be used to continue the simulation
	//!
	//! The format of the checkpointfile written by this method is the same as the format
	//! of the input file.
	//! @param filename Name of the checkpointfile (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
	//! @param useBinaryFormat  Use binary format to write checkpoint files
	CheckpointWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental, bool useBinaryFormat);
	~CheckpointWriter();
	

	/** @brief Read in XML configuration for CheckpointWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="CheckpointWriter">
	     <type>ASCII|binary</type>
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <incremental>INTEGER</incremental>
	     <appendTimestamp>INTEGER</appendTimestamp>
	   </outputplugin>
	   \endcode
	 */
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
		return std::string("CheckpointWriter");
	}
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
    bool    _useBinaryFormat;
	bool	_incremental;
	bool	_appendTimestamp;
};

#endif /*CHECKPOINTWRITER_H_*/
