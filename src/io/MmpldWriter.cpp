#include "io/MmpldWriter.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "/usr/include/endian.h"
#include "utils/FileUtils.h"

// set mmpld file version. possible values: 100 or 102
#define MMPLD_FILE_VERSION 100

using Log::global_log;
using namespace std;

MmpldWriter::MmpldWriter()
		:	_startTimestep(0), _writeFrequency(1000), _stopTimestep(0), _numFramesPerFile(0), _outputPrefix("unknown"),
		 	_bInitSphereData(ISD_READ_FROM_XML), _bWriteControlPrepared(false), _bInitFrameDone(false), _nFileIndex(0), _numFiles(1), _nextRecSimstep(0),
		 	_strOutputPrefixCurrent("unknown"), _frameCountMax(1)
{
	//	if (outputPrefix == "default") <-- what for??
		if(false)
			_appendTimestamp = true;
		else
			_appendTimestamp = false;

		// ERROR
		if (0 == _writeFrequency)
			mardyn_exit(-1);
}

MmpldWriter::MmpldWriter(uint64_t startTimestep, uint64_t writeFrequency, uint64_t stopTimestep, uint64_t numFramesPerFile,
		std::string outputPrefix)
		:	_startTimestep(startTimestep), _writeFrequency(writeFrequency), _stopTimestep(stopTimestep), _numFramesPerFile(numFramesPerFile), _outputPrefix(outputPrefix),
		 	_bInitSphereData(ISD_READ_FROM_XML), _bWriteControlPrepared(false), _bInitFrameDone(false), _nFileIndex(0), _numFiles(1), _nextRecSimstep(startTimestep),
		 	_strOutputPrefixCurrent("unknown"), _frameCountMax(1)
{
//	if (outputPrefix == "default") <-- what for??
	if(false)
		_appendTimestamp = true;
	else
		_appendTimestamp = false;

	// ERROR
	if (0 == _writeFrequency)
		mardyn_exit(-1);
}

MmpldWriter::~MmpldWriter()
{
}

void MmpldWriter::readXML(XMLfileUnits& xmlconfig)
{
	// write control
	xmlconfig.getNodeValue("writecontrol/start", _startTimestep);
	xmlconfig.getNodeValue("writecontrol/writefrequency", _writeFrequency);
	xmlconfig.getNodeValue("writecontrol/stop", _stopTimestep);
	xmlconfig.getNodeValue("writecontrol/framesperfile", _numFramesPerFile);
	global_log->info() << "Start sampling from simstep: " << _startTimestep << endl;
	global_log->info() << "Write with frequency: " << _writeFrequency << endl;
	global_log->info() << "Stop sampling at simstep: " << _stopTimestep << endl;
	global_log->info() << "Split files every " << _numFramesPerFile << "th frame."<< endl;

	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;
	
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;

	// sphere params: radius, colors
	uint32_t numSites = 0;
	XMLfile::Query query = xmlconfig.query("spheres/site");
	numSites = query.card();
	global_log->info() << "Number of sites: " << numSites << endl;
	if(numSites < 1) {
		global_log->warning() << "No site parameters specified." << endl;
	}
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator outputSiteIter;
	for( outputSiteIter = query.begin(); outputSiteIter; outputSiteIter++ )
	{
		xmlconfig.changecurrentnode( outputSiteIter );
		float radius;
		uint32_t r, g, b, alpha;
		xmlconfig.getNodeValue("radius", radius);
		xmlconfig.getNodeValue("color/r", r);
		xmlconfig.getNodeValue("color/g", g);
		xmlconfig.getNodeValue("color/b", b);
		xmlconfig.getNodeValue("color/alpha", alpha);

		_vfSphereRadius.push_back(radius);
		std::array<uint32_t, 4> arrColors = {r, g, b, alpha};
		_vaSphereColors.push_back(arrColors);
	}
}

//Header Information
void MmpldWriter::initOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
	// only executed once
	this->PrepareWriteControl();

	_frameCount = 0;
	_strOutputPrefixCurrent = _vecFilePrefixes.at(_nFileIndex);
	_frameCountMax = _vecFramesPerFile.at(_nFileIndex);

	// number of components / sites
	_numComponents = (uint8_t)domain->getNumberOfComponents();
	_numSitesPerComp.resize(_numComponents);
	_nCompSitesOffset.resize(_numComponents);
	_numSitesTotal = 0;

	vector<Component>& vComponents = *(_simulation.getEnsemble()->getComponents() );
	vector<Component>::iterator cit;
	cit=vComponents.begin();

	_nCompSitesOffset.at(0) = 0;
	for(uint8_t cid=0; cid<_numComponents; cid++)
	{
		uint8_t numSites = (*cit).numLJcenters();
		_numSitesPerComp.at(cid) = numSites;
		if(cid>0)
			_nCompSitesOffset.at(cid) = _nCompSitesOffset.at(cid-1) + _numSitesPerComp.at(cid-1);
		// total number of sites (all components)
		_numSitesTotal += numSites;
		cit++;
	}

#ifndef NDEBUG
	for(uint8_t cid=0; cid<_numComponents; cid++)
		cout << "_nCompSitesOffset[" << (uint32_t)cid << "] = " << (uint32_t)_nCompSitesOffset.at(cid) << endl;
#endif

	// init radius and color of spheres
	this->InitSphereData();
	this->SetNumSphereTypes();

	stringstream filenamestream;
	filenamestream << _strOutputPrefixCurrent;

	/*
	 * not in use actually
	 *
	if(_appendTimestamp) {
		filenamestream << "-" << _timestampString;
	}
	*/
	filenamestream << ".mmpld";

	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

	int rank = domainDecomp->getRank();
#ifndef NDEBUG
	cout << "rank[" << rank << "]: _nFileIndex = " << (uint32_t)_nFileIndex
			<< ", filename = " << filenamestream.str() << endl;
#endif

#ifdef ENABLE_MPI
//	int rank = domainDecomp->getRank();
	if (rank == 0){
#endif
	ofstream mmpldfstream(filename, ios::binary|ios::out);

	//format marker
	uint8_t magicIdentifier[6] = {0x4D, 0x4D, 0x50, 0x4C, 0x44, 0x00};
	mmpldfstream.write((char*)magicIdentifier, sizeof(magicIdentifier));

	//version number
	uint16_t versionNumber;
	switch (MMPLD_FILE_VERSION)
	{
	case 100:
		versionNumber = htole16(100);
		break;

	case 102:
		versionNumber = htole16(102);
		break;

	default:
		cout << "Error mmpld-writer: file version " << MMPLD_FILE_VERSION << " not supported." << endl;
		return;
		break;
	}

	mmpldfstream.write((char*)&versionNumber, sizeof(versionNumber));

	//calculate the number of frames
	uint32_t numframes;
	uint32_t numframes_le;
	numframes = _vecFramesPerFile.at(_nFileIndex);

	_numSeekEntries = numframes+1;
	numframes_le = htole32(numframes);
	mmpldfstream.write((char*)&numframes_le,sizeof(numframes_le));

	_seekTable.resize(_numSeekEntries);

	//boundary box
	float minbox[3] = {0, 0, 0};
	float maxbox[3];
	for (unsigned short d = 0; d < 3; ++d) maxbox[d] = domain->getGlobalLength(d);
	mmpldfstream.write((char*)&minbox,sizeof(minbox));
	mmpldfstream.write((char*)&maxbox,sizeof(maxbox));

	//clipping box
	float inflateRadius = *(_vfSphereRadius.begin() );
	std::vector<float>::iterator it;
	for(it=_vfSphereRadius.begin(); it!=_vfSphereRadius.end(); ++it)
		if(inflateRadius < (*it) ) inflateRadius = (*it);

	for (unsigned short d = 0; d < 3; ++d){
		maxbox[d] = maxbox[d] + inflateRadius;
		minbox[d] = minbox[d] - inflateRadius;
	}
	mmpldfstream.write((char*)&minbox,sizeof(minbox));
	mmpldfstream.write((char*)&maxbox,sizeof(maxbox));

	//preallocate seektable
	uint64_t seekNum = 0;
	for (uint32_t i = 0; i <= numframes; ++i)
		mmpldfstream.write((char*)&seekNum,sizeof(seekNum));

	mmpldfstream.close();
#ifdef ENABLE_MPI
	}
#endif

	// eventually write initial frame
	if(false == _bInitFrameDone && (_startTimestep == _simulation.getNumInitTimesteps() ) )
		this->doOutput(particleContainer, domainDecomp, domain, _nextRecSimstep, NULL, NULL);
	_bInitFrameDone = true;
}

void MmpldWriter::doOutput( ParticleContainer* particleContainer,
		   DomainDecompBase* domainDecomp, Domain* domain,
		   unsigned long simstep, std::list<ChemicalPotential>* /*lmu*/,
		   map<unsigned, CavityEnsemble>* /*mcav*/)
{
	// control writing with desired write frequency
	if( _nextRecSimstep != simstep)
		return;
	_nextRecSimstep += _writeFrequency;

	stringstream filenamestream, outputstream;
	filenamestream << _strOutputPrefixCurrent;

	/*
	 * not in use actually
	 *
	if(_appendTimestamp) {
		filenamestream << "-" << _timestampString;
	}
	*/
	filenamestream << ".mmpld";

	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	int numprocs = domainDecomp->getNumProcs();
	unsigned long numberParticles = particleContainer->getNumberOfParticles();
	long outputsize = 0;

	uint64_t numSpheresPerType[_numSphereTypes];
	for (uint8_t ti = 0; ti < _numSphereTypes; ++ti){
		numSpheresPerType[ti] = 0;
	}

	//calculate number of spheres per component|siteType
	uint32_t molcid = 0;
	/*
	 * should use this version in future
	 */
	for (ParticleIterator mol = particleContainer->iteratorBegin(); mol != particleContainer->iteratorEnd(); ++mol)
			this->CalcNumSpheresPerType(numSpheresPerType, &(*mol));
	/*
	for (Molecule* mol = particleContainer->begin(); mol != particleContainer->end(); mol = particleContainer->next())
		this->CalcNumSpheresPerType(numSpheresPerType, mol);
*/
	//distribute global component particle count
	uint64_t globalNumCompSpheres[_numSphereTypes];
	if (rank == 0){
		for (uint32_t i = 0; i < _numSphereTypes; ++i){
			globalNumCompSpheres[i] = numSpheresPerType[i];
		}
		MPI_Status status;
		uint64_t numSpheresPerTypeTmp[_numSphereTypes];
		for (int source = rank+1; source < numprocs; ++source){
			int recvcount = sizeof(numSpheresPerTypeTmp);
			int recvtag = 1;
			MPI_Recv(numSpheresPerTypeTmp, recvcount, MPI_BYTE, source, recvtag, MPI_COMM_WORLD, &status);
			for (uint8_t ti = 0; ti < _numSphereTypes; ++ti)
				globalNumCompSpheres[ti] = globalNumCompSpheres[ti] + numSpheresPerTypeTmp[ti];
		}
	}else{
			int dest = 0;
			int sendcount = sizeof(numSpheresPerType);
			int sendtag = 1;
			MPI_Request request;
			MPI_Isend(numSpheresPerType, sendcount, MPI_BYTE, dest, sendtag, MPI_COMM_WORLD, &request);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY|MPI_MODE_APPEND|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);


	//write particle list for each component|site (sphere type)
	for (uint8_t nSphereTypeIndex=0; nSphereTypeIndex<_numSphereTypes; ++nSphereTypeIndex){
		//add space for particle data
		outputsize = (long)numSpheresPerType[nSphereTypeIndex]*12;
		
		//add space for particle list header
		if (rank == 0){
			//add particle list header
			outputsize += 18;
			if (nSphereTypeIndex == 0){

				switch (MMPLD_FILE_VERSION){
					case 100:
						//add space for number of particle lists
						outputsize += 4;
						break;

					case 102:
						//add space for timestamp and number of particle lists
						outputsize += 8;
						break;

					default:
						cout << "Error mmpld-writer: file version " << MMPLD_FILE_VERSION << " not supported." << endl;
						return;
						break;
				}
			}
		}
		
		//send outputsize of current rank to next rank
		for (int dest = rank+1; dest < numprocs; ++dest){
			int sendcount = 1;
			int sendtag = 0;
			MPI_Request request;
			MPI_Isend(&outputsize, sendcount, MPI_LONG, dest, sendtag, MPI_COMM_WORLD, &request);
		}
		//accumulate outputsizes of previous ranks and use it as offset for output file
		MPI_Status status;
		long offset = 0;
		long outputsize_get;
		for (int source = 0; source < rank; ++source){
			int recvcount = 1;
			int recvtag = 0;
			MPI_Recv(&outputsize_get, recvcount, MPI_LONG, source, recvtag, MPI_COMM_WORLD, &status);
			offset += outputsize_get;
		}

		global_log->debug() << "MmpldWriter rank: " << rank << "; step: " << simstep << "; sphereTypeIndex: " << nSphereTypeIndex << "; offset: " << offset << endl;

		MPI_File_seek(fh, offset, MPI_SEEK_END);

		MPI_Barrier(MPI_COMM_WORLD);

		//write particle list header
		if (rank == 0){
			
			//write frame header if we are before the first particle list
			if (nSphereTypeIndex == 0){
				//store file position for seek table
				if (_frameCount < _numSeekEntries){
					MPI_Offset entry;
					MPI_File_get_position(fh, &entry);
					_seekTable.at(_frameCount) = (uint64_t)entry;
				}
//					_frameCount = _frameCount + 1;
				
				float frameHeader_timestamp = simstep;
				
				switch (MMPLD_FILE_VERSION){
					case 100:
						//do not write timestamp to frame header
						break;

					case 102:
						//write timestamp to frame header
						MPI_File_write(fh, &frameHeader_timestamp, 1, MPI_FLOAT, &status);
						break;

					default:
						cout << "Error mmpld-writer: file version " << MMPLD_FILE_VERSION << " not supported." << endl;
						return;
						break;
				}
				
				uint32_t frameHeader_numPLists = htole32(_numSphereTypes);
				MPI_File_write(fh, &frameHeader_numPLists, 1, MPI_UNSIGNED, &status);
				
			}
			uint8_t pListHeader_vortexType;
			uint8_t pListHeader_colorType;
			float pListHeader_globalRadius;
			uint8_t pListHeader_red;
			uint8_t pListHeader_green;
			uint8_t pListHeader_blue;
			uint8_t pListHeader_alpha;
			uint64_t pListHeader_particleCount;

			//set vortex data type to FLOAT_XYZ
			pListHeader_vortexType = 1;
			
			//set color data type to NONE (only global color used)
			pListHeader_colorType = 0;

			//select different colors depending on component|site id
			pListHeader_globalRadius = _vfSphereRadius[nSphereTypeIndex];
			pListHeader_red   = _vaSphereColors[nSphereTypeIndex][0];
			pListHeader_green = _vaSphereColors[nSphereTypeIndex][1];
			pListHeader_blue  = _vaSphereColors[nSphereTypeIndex][2];
			pListHeader_alpha = _vaSphereColors[nSphereTypeIndex][3];

			//store componentParticleCount
			pListHeader_particleCount = htole64(globalNumCompSpheres[nSphereTypeIndex]);

			MPI_File_write(fh, &pListHeader_vortexType, 1, MPI_BYTE, &status);
			MPI_File_write(fh, &pListHeader_colorType, 1, MPI_BYTE, &status);
			MPI_File_write(fh, &pListHeader_globalRadius, 1, MPI_FLOAT, &status);
			MPI_File_write(fh, &pListHeader_red, 1, MPI_BYTE, &status);
			MPI_File_write(fh, &pListHeader_green, 1, MPI_BYTE, &status);
			MPI_File_write(fh, &pListHeader_blue, 1, MPI_BYTE, &status);
			MPI_File_write(fh, &pListHeader_alpha, 1, MPI_BYTE, &status);
			MPI_File_write(fh, &pListHeader_particleCount, 1, MPI_LONG_LONG_INT, &status);
		}  // if (rank == 0){

		float spherePos[3];
		/*
		 * should use this version in future
		 */
		for (ParticleIterator mol = particleContainer->iteratorBegin(); mol != particleContainer->iteratorEnd(); ++mol)
		{
			if(true == GetSpherePos(spherePos, &(*mol), nSphereTypeIndex) )
		/*
		for (Molecule* mol = particleContainer->begin(); mol != particleContainer->end(); mol = particleContainer->next())
		{
			if(true == GetSpherePos(spherePos, mol, nSphereTypeIndex) )*/
				MPI_File_write(fh, spherePos, 3, MPI_FLOAT, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}  // for (uint8_t nSphereTypeIndex=0; nSphereTypeIndex<_numSphereTypes; ++nSphereTypeIndex)

	// data of frame is written
	_frameCount++;
#ifndef NDEBUG
	cout << "rank[" << rank << "]: _frameCount = " << _frameCount << endl;
#endif

	// write seek table entry
	if (rank == 0)
	{
		MPI_Status status;
		uint64_t seektablePos = 0x3C+(sizeof(uint64_t)*(_frameCount-1));
		MPI_File_seek(fh, seektablePos, MPI_SEEK_SET);
		uint64_t seekPosition;
		seekPosition = htole64(_seekTable.at(_frameCount-1) );
		MPI_File_write(fh, &seekPosition, 1, MPI_LONG_LONG_INT, &status);
		uint32_t frameCount = htole32(_frameCount-1);  // last frame will be ignored when simulation is aborted
		MPI_File_seek(fh, 0x08, MPI_SEEK_SET);  // 0x08: frame count position in file header
		MPI_File_write(fh, &frameCount, 1, MPI_UNSIGNED, &status);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_close(&fh);
#endif

	// split data to multiple files
	if(_frameCount == _frameCountMax && (_nFileIndex+1) < _numFiles)
		this->MultiFileApproachReset(particleContainer, domainDecomp, domain);
}

void MmpldWriter::finishOutput(ParticleContainer* /*particleContainer*/, DomainDecompBase* domainDecomp, Domain* /*domain*/) {
	//fill seektable
	stringstream filenamestream;
	filenamestream << _vecFilePrefixes.at(_nFileIndex);

	/*
	 * not in use actually
	 *
	if(_appendTimestamp) {
		filenamestream << "-" << _timestampString;
	}
	*/
	filenamestream << ".mmpld";

	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
		
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_seek(fh, 0, MPI_SEEK_END);
		MPI_Offset endPosition;
		MPI_File_get_position(fh, &endPosition);
		_seekTable.at(_numSeekEntries-1) = (uint64_t)endPosition;
		uint64_t seektablePos = 0x3C+(sizeof(uint64_t)*(_numSeekEntries-1));
		MPI_File_seek(fh, seektablePos, MPI_SEEK_SET);
		uint64_t seekPosition;
		MPI_Status status;
		seekPosition = htole64(_seekTable.at(_numSeekEntries-1) );
		MPI_File_write(fh, &seekPosition, 1, MPI_LONG_LONG_INT, &status);
		uint32_t frameCount = htole32(_frameCount);  // set final number of frames
		MPI_File_seek(fh, 0x08, MPI_SEEK_SET);  // 0x08: frame count position in file header
		MPI_File_write(fh, &frameCount, 1, MPI_UNSIGNED, &status);
		MPI_File_close(&fh);
	}else{
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_close(&fh);
	}

#endif
}

void MmpldWriter::InitSphereData()
{
	if(_bInitSphereData == ISD_READ_FROM_XML)
		return;
	else if(_bInitSphereData == ISD_USE_DEFAULT)
	{
		for(uint8_t i=0; i<6; i++)
			_vfSphereRadius.push_back(0.5);

		//                            R    G    B  alpha
		std::array<uint32_t, 4> red = { 255, 0, 0, 255 };
		_vaSphereColors.push_back(red);
		std::array<uint32_t, 4> lightblue = { 0, 205, 255, 255 };
		_vaSphereColors.push_back(lightblue);
		std::array<uint32_t, 4> blue = { 255, 0, 255, 255 };
		_vaSphereColors.push_back(blue);
		std::array<uint32_t, 4> green = { 0, 155, 0, 255 };
		_vaSphereColors.push_back(green);
		std::array<uint32_t, 4> purple = { 105, 0, 205, 255 };
		_vaSphereColors.push_back(purple);
		std::array<uint32_t, 4> orange = { 255, 125, 0, 255 };
		_vaSphereColors.push_back(orange);

		return;
	}

	std::ifstream filein(_strSphereDataFilename.c_str(), ios::in);
	std::string strLine, strToken;
	std::string strTokens[6];
	std::array<uint32_t, 4> arrColors;

	while (getline (filein, strLine))
	{
		stringstream sstr;
		sstr << strLine;
//		cout << sstr.str() << endl;

		uint8_t ti=0;
		while (sstr >> strToken)
			if(ti<6) strTokens[ti++] = strToken;

		if(ti==6 && strTokens[0][0] != '#')
		{
			_vfSphereRadius.push_back( (float)(atof( strTokens[1].c_str() ) ) );
			for(uint8_t ci=0; ci<4; ci++)
				arrColors[ci] = (uint32_t)(atoi(strTokens[ci+2].c_str() ) );
			_vaSphereColors.push_back(arrColors);
		}
	}

#ifndef NDEBUG
	std::vector<float>::iterator it; int i=0; cout << "radii" << endl;
	for(it=_vfSphereRadius.begin(); it!=_vfSphereRadius.end(); ++it)
		cout << i++ << ": " << (*it) << endl;

	std::vector< std::array<uint32_t, 4> >::iterator cit; i=0; cout << "colors" << endl;
	for(cit=_vaSphereColors.begin(); cit!=_vaSphereColors.end(); ++cit)
	{
		cout << i++ << ":";
		for(int j=0; j<4; j++)
			 cout << setw(4) << (int)(*cit).data()[j];
		cout << endl;
	}
#endif
}

void MmpldWriter::MultiFileApproachReset(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
	this->finishOutput(particleContainer, domainDecomp, domain);
	_nFileIndex++;
	this->initOutput(particleContainer, domainDecomp, domain);
}

void MmpldWriter::PrepareWriteControl()
{
	// this method should only be executed once
	if(true == _bWriteControlPrepared)
		return;
	_bWriteControlPrepared = true;

	_startTimestep = ( _startTimestep > _simulation.getNumInitTimesteps() ) ? _startTimestep : _simulation.getNumInitTimesteps();
	_nextRecSimstep = _startTimestep;
	if( 0 == _stopTimestep )
		_stopTimestep = _simulation.getNumTimesteps();
	else
		_stopTimestep  = ( _stopTimestep  < _simulation.getNumTimesteps() ) ? _stopTimestep  : _simulation.getNumTimesteps();

	uint64_t numTimesteps = _stopTimestep - _startTimestep;
	uint64_t numFramesTotal = numTimesteps/_writeFrequency;
	if(_numFramesPerFile >= numFramesTotal || _numFramesPerFile == 0)
		_numFramesPerFile = numFramesTotal;

#ifndef NDEBUG
	cout << "_startTimestep    = " << _startTimestep << endl;
	cout << "_writeFrequency   = " << _writeFrequency << endl;
	cout << "_stopTimestep     = " << _stopTimestep << endl;
	cout << "_numFramesPerFile = " << _numFramesPerFile << endl;
	cout << "numTimesteps      = " << numTimesteps << endl;
	cout << "numFramesTotal    = " << numFramesTotal << endl;
#endif

	// init frames per file vector
	_numFiles = numFramesTotal/_numFramesPerFile;
	_vecFramesPerFile.insert (_vecFramesPerFile.begin(), _numFiles, _numFramesPerFile);
	_vecFramesPerFile.at(0) += 1;  // +1 frame in first file: start configuration

	// init file prefix vector
	for(uint8_t fi=0; fi<_numFiles; ++fi)
	{
		std::stringstream sstrPrefix;
		sstrPrefix << _outputPrefix << "_" << fill_width('0', 4) << (uint32_t)(fi+1);
		_vecFilePrefixes.push_back(sstrPrefix.str() );
	}

	// handle the case: numFramesTotal not a multiple of _numFramesPerFile
	if(0 != numFramesTotal % _numFramesPerFile)
	{
		_numFiles++;
		_vecFramesPerFile.push_back(numFramesTotal % _numFramesPerFile);
		std::stringstream sstrPrefix;
		sstrPrefix << _outputPrefix << "_" << fill_width('0', 4) << (uint32_t)(_numFiles);
		_vecFilePrefixes.push_back(sstrPrefix.str() );
	}

#ifndef NDEBUG
	for(auto fit : _vecFilePrefixes)
		cout << fit << endl;
	for(auto fit : _vecFramesPerFile)
		cout << fit << endl;
#endif
}

// derived classes
void MmpldWriterSimpleSphere::CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol)
{
	uint8_t cid = mol->componentid();
	numSpheresPerType[cid]++;
}

bool MmpldWriterSimpleSphere::GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex)
{
	uint8_t cid = mol->componentid();
	for (unsigned short d = 0; d < 3; ++d) spherePos[d] = (float)mol->r(d);
	return (cid == nSphereTypeIndex);
}


void MmpldWriterMultiSphere::CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol)
{
	uint8_t cid = mol->componentid();
	uint8_t offset = _nCompSitesOffset.at(cid);
	for (uint8_t si = 0; si < _numSitesPerComp.at(cid); ++si)
		numSpheresPerType[offset+si]++;
}

bool MmpldWriterMultiSphere::GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex)
{
	bool ret = false;
	uint8_t cid = mol->componentid();
	uint8_t numSites =  _numSitesPerComp.at(cid);
	uint8_t offset  = _nCompSitesOffset.at(cid);
	for (uint8_t si = 0; si < numSites; ++si)
	{
		if(offset+si == nSphereTypeIndex)
		{
			const std::array<double,3> arrSite = mol->ljcenter_d_abs(si);
			const double* posSite = arrSite.data();
			for (unsigned short d = 0; d < 3; ++d) spherePos[d] = (float)posSite[d];
			ret = true;
		}
	}
	return ret;
}



















