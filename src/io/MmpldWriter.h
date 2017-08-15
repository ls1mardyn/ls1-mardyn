#ifndef MMPLDWRITER_H_
#define MMPLDWRITER_H_

#include <string>
#include <vector>
#include <array>

#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"

#ifdef ENABLE_MPI
#include "utils/MPI_Info_object.h"
#endif

enum InitSphereData : uint8_t
{
	ISD_USE_DEFAULT = 1,
	ISD_READ_FROM_FILE = 2,
	ISD_READ_FROM_XML = 3,
};

enum MMPLD_Vertex_type : uint8_t {
	MMPLD_VERTEX_NONE       = 0,
	MMPLD_VERTEX_FLOAT_XYZ  = 1,
	MMPLD_VERTEX_FLOAT_XYZR = 2,
	MMPLD_VERTEX_SHORT_XYZ  = 3,
};

enum MMPLD_Color_type : uint8_t {
	MMPLD_COLOR_NONE       = 0,
	MMPLD_COLOR_UINT8_RGB  = 1,
	MMPLD_COLOR_UINT8_RGBA = 2,
	MMPLD_COLOR_FLOAT_I    = 3,
	MMPLD_COLOR_FLOAT_RGB  = 4,
	MMPLD_COLOR_FLOAT_RGBA = 5,
};

class Simulation;

/** @brief Output plugin to generate a MegaMol™ Particle List Data file (*.mmpld).
 */
class MmpldWriter : public OutputBase
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
	MmpldWriter(uint64_t startTimestep, uint64_t writeFrequency, uint64_t stopTimestep, uint64_t numFramesPerFile,
			std::string outputPrefix);
	virtual ~MmpldWriter() {}

	virtual void SetNumSphereTypes() {};
	virtual void CalcNumSpheresPerType(ParticleContainer* particleContainer, uint64_t* numSpheresPerType) {};
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex) {};

	void InitSphereData();

public:
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
	
	static std::string getPluginName() {
		return std::string("MmpldWriter");
	}
	static OutputBase* createInstance() { return new MmpldWriter(); }

	void SetInitSphereDataParameters(const uint8_t &bInitSphereData, const std::string &strSphereDataFilename) {
		_bInitSphereData = bInitSphereData; _strSphereDataFilename = strSphereDataFilename;
	}

protected:
	void MultiFileApproachReset(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void PrepareWriteControl();
	long get_data_frame_header_size();
	long get_data_list_header_size();
	long get_data_list_size(uint64_t particle_count);

protected:
	/** First time step to be recorded */
	uint64_t _startTimestep;
	/** time steps between two records*/
	uint64_t _writeFrequency;
	/** Max time step up to which shall be recorded */
	uint64_t _stopTimestep;
	uint64_t _numFramesPerFile;
	std::string _outputPrefix;
	std::string _timestampString;
	uint32_t _numSeekEntries;
	uint32_t _frameCount;
	uint8_t  _numComponents;
	uint8_t  _numSitesTotal;
	uint8_t  _numSphereTypes;
	std::vector<uint64_t> _seekTable;
	std::vector<uint8_t> _numSitesPerComp;
	std::vector<uint8_t> _nCompSitesOffset;
	std::vector<float> _vfSphereRadius;
	std::vector< std::array<uint32_t, 4> > _vaSphereColors;
	std::string _strSphereDataFilename;
	uint8_t _bInitSphereData;
	bool _bWriteControlPrepared;

	long _fileCount;
	int _mmpldversion;
	MMPLD_Vertex_type _vertex_type;
	MMPLD_Color_type _color_type;

	std::string getOutputFilename();

#if ENABLE_MPI
	MPI_Info_object _mpiinfo;
#endif
};

class MmpldWriterSimpleSphere : public MmpldWriter
{
public:
	MmpldWriterSimpleSphere() {}
	MmpldWriterSimpleSphere(uint64_t startTimestep, uint64_t writeFrequency, uint64_t stopTimestep, uint64_t numFramesPerFile,
			std::string outputPrefix)
			: MmpldWriter(startTimestep, writeFrequency, stopTimestep, numFramesPerFile, outputPrefix)
	{
//		MmpldWriter::_numSphereTypes = &(MmpldWriter::_numComponents);
	}
	virtual ~MmpldWriterSimpleSphere() {}

	virtual void SetNumSphereTypes() {_numSphereTypes = _numComponents;}
	virtual void CalcNumSpheresPerType(ParticleContainer* particleContainer, uint64_t* numSpheresPerType);
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex);
};

class MmpldWriterMultiSphere : public MmpldWriter
{
public:
	MmpldWriterMultiSphere() {}
	MmpldWriterMultiSphere(uint64_t startTimestep, uint64_t writeFrequency, uint64_t stopTimestep, uint64_t numFramesPerFile,
			std::string outputPrefix)
			: MmpldWriter(startTimestep, writeFrequency, stopTimestep, numFramesPerFile, outputPrefix)
	{
//		MmpldWriter::_numSphereTypes = &(MmpldWriter::_numSitesTotal);
	}
	virtual ~MmpldWriterMultiSphere() {}

	virtual void SetNumSphereTypes() {_numSphereTypes = _numSitesTotal;}
	virtual void CalcNumSpheresPerType(ParticleContainer* particleContainer, uint64_t* numSpheresPerType);
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex);
};

#endif /* MMPLDWRITER_H_ */




