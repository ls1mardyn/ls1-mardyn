#ifndef MMPLDWRITER_H_
#define MMPLDWRITER_H_

#include <string>
#include <vector>
#include <array>

#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"

enum InitSphereData : uint8_t
{
	ISD_USE_DEFAULT = 1,
	ISD_READ_FROM_FILE = 2,
	ISD_READ_FROM_XML = 3,
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
	virtual ~MmpldWriter();
	virtual void SetNumSphereTypes() = 0;
	virtual void CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol) = 0;
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex) = 0;

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
	
	std::string getPluginName() {
		return std::string("MmpldWriter");
	}

	void SetInitSphereDataParameters(const uint8_t &bInitSphereData, const std::string &strSphereDataFilename) {
		_bInitSphereData = bInitSphereData; _strSphereDataFilename = strSphereDataFilename;
	}

protected:
	void MultiFileApproachReset(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void PrepareWriteControl();

protected:
	uint64_t _startTimestep;
	uint64_t _writeFrequency;
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
	bool _bInitFrameDone;

	// split files
	uint8_t _nFileIndex;
	uint8_t _numFiles;
	std::vector<string> _vecFilePrefixes;
	std::vector<uint64_t> _vecFramesPerFile;
	uint64_t _nextRecSimstep;
	std::string _strOutputPrefixCurrent;
	uint32_t _frameCountMax;
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
	virtual void CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol);
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
	virtual void CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol);
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex);
};

#endif /* MMPLDWRITER_H_ */




