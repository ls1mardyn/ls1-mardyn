#ifndef MMPLDWRITER_H_
#define MMPLDWRITER_H_

#include <string>
#include <vector>
#include <array>
#include <cstdint>

#include "plugins/PluginBase.h"
#include "molecules/MoleculeForwardDeclaration.h"

#ifdef ENABLE_MPI
#include "utils/MPI_Info_object.h"
#endif

enum InitSphereData : std::uint8_t
{
	ISD_USE_DEFAULT = 1,
	ISD_READ_FROM_FILE = 2,
	ISD_READ_FROM_XML = 3,
};

enum MMPLD_Vertex_type : std::uint8_t {
	MMPLD_VERTEX_NONE       = 0,
	MMPLD_VERTEX_FLOAT_XYZ  = 1,
	MMPLD_VERTEX_FLOAT_XYZR = 2,
	MMPLD_VERTEX_SHORT_XYZ  = 3,
};

enum MMPLD_Color_type : std::uint8_t {
	MMPLD_COLOR_NONE       = 0,
	MMPLD_COLOR_UINT8_RGB  = 1,
	MMPLD_COLOR_UINT8_RGBA = 2,
	MMPLD_COLOR_FLOAT_I    = 3,
	MMPLD_COLOR_FLOAT_RGB  = 4,
	MMPLD_COLOR_FLOAT_RGBA = 5,
};

/** @brief Output plugin to generate a MegaMol™ Particle List Data file (*.mmpld).
 */
class MmpldWriter : public PluginBase
{
protected:
	MmpldWriter();
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
	MmpldWriter(std::uint64_t startTimestep, std::uint64_t writeFrequency, std::uint64_t stopTimestep, std::uint64_t numFramesPerFile,
			std::string outputPrefix);
	virtual ~MmpldWriter() {}

	virtual void SetNumSphereTypes() {}
	virtual void CalcNumSpheresPerType(ParticleContainer* particleContainer, std::uint64_t* numSpheresPerType) {}
	virtual bool GetSpherePos(float *spherePos, Molecule* mol, std::uint8_t& nSphereTypeIndex) { return false; }

	void InitSphereData();

public:
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
		return std::string("MmpldWriter");
	}
	static PluginBase* createInstance() { return new MmpldWriter(); }

	void SetInitSphereDataParameters(const std::uint8_t &bInitSphereData, const std::string &strSphereDataFilename) {
		_bInitSphereData = bInitSphereData; _strSphereDataFilename = strSphereDataFilename;
	}

protected:
	std::string getOutputFilename();
	void MultiFileApproachReset(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void PrepareWriteControl();
	long get_data_frame_header_size();
	long get_seekTable_size();
	void writeSeekTableEntry(int id, std::uint64_t offset);
	long get_data_list_header_size();
	long get_particle_data_size();
	long get_data_list_size(std::uint64_t particle_count);
	void write_frame_header(std::uint32_t num_data_lists);
	void write_particle_list_header(std::uint64_t particle_count, int sphereId);
	void write_frame(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp);

protected:
	/** First time step to be recorded */
	std::uint64_t _startTimestep;
	/** time steps between two records*/
	std::uint64_t _writeFrequency;
	/** Max time step up to which shall be recorded */
	std::uint64_t _stopTimestep;
	long _writeBufferSize;
	std::string _outputPrefix;
	std::string _timestampString;
	std::uint32_t _frameCount;
	std::uint8_t  _numComponents;
	std::uint8_t  _numSitesTotal;
	std::uint8_t  _numSphereTypes;
	std::vector<std::uint64_t> _seekTable;
	std::vector<std::uint8_t> _numSitesPerComp;
	std::vector<std::uint8_t> _nCompSitesOffset;
	std::string _strSphereDataFilename;
	std::uint8_t _bInitSphereData;
	bool _bWriteControlPrepared;

	long _fileCount;
	std::uint32_t _numFramesPerFile;

	std::uint16_t _mmpldversion;
	std::uint32_t _numSeekEntries;
	MMPLD_Vertex_type _vertex_type;
	MMPLD_Color_type _color_type;
	std::vector<float> _global_radius;
	std::vector< std::array<std::uint8_t, 4> > _global_rgba;
	std::vector< std::array<float, 2> > _global_intensity_range;


#ifdef ENABLE_MPI
	MPI_File _mpifh;
	MPI_Info_object _mpiinfo;
#endif
};

class MmpldWriterSimpleSphere : public MmpldWriter
{
public:
	MmpldWriterSimpleSphere() {}
	MmpldWriterSimpleSphere(std::uint64_t startTimestep, std::uint64_t writeFrequency, std::uint64_t stopTimestep, std::uint64_t numFramesPerFile,
			std::string outputPrefix)
			: MmpldWriter(startTimestep, writeFrequency, stopTimestep, numFramesPerFile, outputPrefix)
	{
//		MmpldWriter::_numSphereTypes = &(MmpldWriter::_numComponents);
	}
	virtual ~MmpldWriterSimpleSphere() {}

	virtual void SetNumSphereTypes() {_numSphereTypes = _numComponents;}
	virtual void CalcNumSpheresPerType(ParticleContainer* particleContainer, std::uint64_t* numSpheresPerType);
	virtual bool GetSpherePos(float *spherePos, Molecule* mol, std::uint8_t& nSphereTypeIndex);
};

class MmpldWriterMultiSphere : public MmpldWriter
{
public:
	MmpldWriterMultiSphere() {}
	MmpldWriterMultiSphere(std::uint64_t startTimestep, std::uint64_t writeFrequency, std::uint64_t stopTimestep, std::uint64_t numFramesPerFile,
			std::string outputPrefix)
			: MmpldWriter(startTimestep, writeFrequency, stopTimestep, numFramesPerFile, outputPrefix)
	{
//		MmpldWriter::_numSphereTypes = &(MmpldWriter::_numSitesTotal);
	}
	virtual ~MmpldWriterMultiSphere() {}

	virtual void SetNumSphereTypes() {_numSphereTypes = _numSitesTotal;}
	virtual void CalcNumSpheresPerType(ParticleContainer* particleContainer, std::uint64_t* numSpheresPerType);
	virtual bool GetSpherePos(float *spherePos, Molecule* mol, std::uint8_t& nSphereTypeIndex);
};

#endif /* MMPLDWRITER_H_ */
