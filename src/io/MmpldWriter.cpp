#include "io/MmpldWriter.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <endian.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <iostream>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/FileUtils.h"


// default version to use for mmpld format writing. possible values: 100 or 102
#define MMPLD_DEFAULT_VERSION 100
#define MMPLD_HEADER_DATA_SIZE 60

using Log::global_log;
using namespace std;

std::string MmpldWriter::getOutputFilename() {
	std::stringstream filenamestream;
	filenamestream << _outputPrefix << "_" << fill_width('0', 4) << _fileCount << ".mmpld";
	return filenamestream.str();
}

MmpldWriter::MmpldWriter() :
		_startTimestep(0), _writeFrequency(1000), _stopTimestep(0), _numFramesPerFile(0), _outputPrefix("unknown"),
		_bInitSphereData(ISD_READ_FROM_XML), _bWriteControlPrepared(false),
		_fileCount(1), _mmpldversion(MMPLD_DEFAULT_VERSION), _vertex_type(MMPLD_VERTEX_FLOAT_XYZ), _color_type(MMPLD_COLOR_NONE)
{}

MmpldWriter::MmpldWriter(uint64_t startTimestep, uint64_t writeFrequency, uint64_t stopTimestep, uint64_t numFramesPerFile,
		std::string outputPrefix)
		:	_startTimestep(startTimestep), _writeFrequency(writeFrequency), _stopTimestep(stopTimestep), _numFramesPerFile(numFramesPerFile), _outputPrefix(outputPrefix),
			_bInitSphereData(ISD_READ_FROM_XML), _bWriteControlPrepared(false), _fileCount(1), _vertex_type(MMPLD_VERTEX_FLOAT_XYZ), _color_type(MMPLD_COLOR_NONE)
{
	if (0 == _writeFrequency) {
		mardyn_exit(-1);
	}
}

void MmpldWriter::readXML(XMLfileUnits& xmlconfig)
{
	// write control
	xmlconfig.getNodeValue("writecontrol/start", _startTimestep);
	xmlconfig.getNodeValue("writecontrol/writefrequency", _writeFrequency);
	xmlconfig.getNodeValue("writecontrol/stop", _stopTimestep);
	xmlconfig.getNodeValue("writecontrol/framesperfile", _numFramesPerFile);
	global_log->info() << "[MMPLD Writer] Start sampling from simstep: " << _startTimestep << endl;
	global_log->info() << "[MMPLD Writer] Write with frequency: " << _writeFrequency << endl;
	global_log->info() << "[MMPLD Writer] Stop sampling at simstep: " << _stopTimestep << endl;
	global_log->info() << "[MMPLD Writer] Split files every " << _numFramesPerFile << "th frame."<< endl;

	xmlconfig.getNodeValue("mmpldversion", _mmpldversion);
	switch(_mmpldversion) {
		case 100:
		case 102:
			break;
		default:
			global_log->error() << "Unsupported MMPLD version:" << _mmpldversion << endl;
			global_simulation->exit(1);
			break;
	}
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "[MMPLD Writer] Output prefix: " << _outputPrefix << endl;

	// sphere params: radius, colors
	uint32_t numSites = 0;
	XMLfile::Query query = xmlconfig.query("spheres/site");
	numSites = query.card();
	global_log->info() << "[MMPLD Writer] Number of sites: " << numSites << endl;
	if(numSites < 1) {
		global_log->warning() << "[MMPLD Writer] No site parameters specified." << endl;
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

	if(xmlconfig.changecurrentnode("mpi_info")) {
#ifdef ENABLE_MPI
		global_log->info() << "[MMPLD Writer] Setting MPI info object for IO" << endl;
		_mpiinfo.readXML(xmlconfig);
#else
		global_log->info() << "[MMPLD Writer] mpi_info only used in parallel/MPI version" << endl;
#endif
		xmlconfig.changecurrentnode("..");
	}
}

//Header Information
void MmpldWriter::initOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
	// only executed once
	this->PrepareWriteControl();

	_frameCount = 0;

	// number of components / sites
	vector<Component> *components = global_simulation->getEnsemble()->getComponents();
	_numComponents = components->size();
	_numSitesPerComp.resize(_numComponents);
	_nCompSitesOffset.resize(_numComponents);
	_numSitesTotal = 0;

	for(int cid = 0; cid < _numComponents; ++cid) {
		Component &component = components->at(cid);
		/** @todo MMPLD writer takes into account only LJ sites at the moment, here */
		int numSites = component.numLJcenters();
		_numSitesPerComp.at(cid) = numSites;
		_nCompSitesOffset.at(cid) = _numSitesTotal; /* offset is total number of sites so far */
		global_log->debug() << "[MMPLD Writer] Component[" << cid << "] numSites=" << numSites << " offset=" << unsigned(_nCompSitesOffset.at(cid)) << endl;
		_numSitesTotal += numSites;
	}
	global_log->debug() << "[MMPLD Writer] Total number of sites taken into account: " << unsigned(_numSitesTotal) << endl;

	// init radius and color of spheres
	this->InitSphereData();
	this->SetNumSphereTypes();

	string filename = getOutputFilename();


#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
#endif
	ofstream mmpldfstream(filename.c_str(), ios::binary|ios::out);

	//format marker
	uint8_t magicIdentifier[6] = {0x4D, 0x4D, 0x50, 0x4C, 0x44, 0x00};
	mmpldfstream.write((char*)magicIdentifier, sizeof(magicIdentifier));

	uint16_t mmpldversion_little_endian = htole16(_mmpldversion);
	mmpldfstream.write((char*)&mmpldversion_little_endian, sizeof(mmpldversion_little_endian));

	//calculate the number of frames
	uint32_t numframes = _numFramesPerFile;
	uint32_t numframes_le = htole32(numframes);
	mmpldfstream.write((char*)&numframes_le,sizeof(numframes_le));

	_numSeekEntries = numframes + 1; /* need additional seek entry for end of file offset marker */
	_seekTable.resize(_numSeekEntries);

	global_log->debug() << "[MMPLD Writer] Writing bounding box data." << endl;
	float minbox[3] = {0, 0, 0};
	float maxbox[3];
	for (unsigned short d = 0; d < 3; ++d) maxbox[d] = domain->getGlobalLength(d);
	mmpldfstream.write((char*)&minbox,sizeof(minbox));
	mmpldfstream.write((char*)&maxbox,sizeof(maxbox));

	global_log->debug() << "[MMPLD Writer] Writing clipping box data." << endl;
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

	global_log->debug() << "[MMPLD Writer] Preallocating " << _numSeekEntries << " seek table entries for frames" << endl;
	uint64_t seekNum = 0;
	for (uint32_t i = 0; i <= numframes; ++i)
		mmpldfstream.write((char*)&seekNum,sizeof(seekNum));

	mmpldfstream.close();
#ifdef ENABLE_MPI
	}
#endif
}

void MmpldWriter::doOutput( ParticleContainer* particleContainer,
		   DomainDecompBase* domainDecomp, Domain* domain,
		   unsigned long simstep, std::list<ChemicalPotential>* /*lmu*/,
		   map<unsigned, CavityEnsemble>* /*mcav*/)
{
	if((simstep < _startTimestep) || (simstep > _stopTimestep) || (0 != ((simstep - _startTimestep) % _writeFrequency)) ) {
		return;
	}
	if(_frameCount == _numFramesPerFile) {
		MultiFileApproachReset(particleContainer, domainDecomp, domain);  // begin new file
	}

	string filename = getOutputFilename();
	global_log->info() << "[MMPLD Writer] Writing MMPLD frame " << _frameCount << " for simstep " << simstep << " to file " << filename << endl;

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();

	//calculate local number of spheres per component|siteType
	std::vector<uint64_t> numSpheresPerType(_numSphereTypes);
	this->CalcNumSpheresPerType(particleContainer, numSpheresPerType.data());

	//distribute global component particle count
	std::vector<uint64_t> globalNumCompSpheres(_numSphereTypes);
	MPI_Reduce(numSpheresPerType.data(), globalNumCompSpheres.data(), _numSphereTypes, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY|MPI_MODE_APPEND|MPI_MODE_CREATE, _mpiinfo, &fh);

	//write particle list for each component|site (sphere type)
	for (uint8_t nSphereTypeIndex = 0; nSphereTypeIndex < _numSphereTypes; ++nSphereTypeIndex){
		long outputsize = 0;
		if (rank == 0){
			// add space for initial data frame  header
			if (nSphereTypeIndex == 0){
				outputsize += get_data_frame_header_size();
			}
			// add space for particle list header
			outputsize += get_data_list_header_size();
		}
		//add space for particle data
		outputsize += get_data_list_size(numSpheresPerType[nSphereTypeIndex]); /* add space for local size */

		MPI_Status status;

		//accumulate outputsizes of previous ranks and use it as offset for output file
		long offset = 0;
		MPI_Exscan(&outputsize, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		if(rank == 0) {offset = 0;} /* state of offset undefined */
		global_log->debug() << "[MMPLD Writer] rank: " << rank << "; step: " << simstep << "; sphereTypeIndex: " << unsigned(nSphereTypeIndex) << "; offset: " << offset << endl;

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
				
				switch (_mmpldversion){
					case 100:
						//do not write timestamp to frame header
						break;
					case 102:
						//write timestamp to frame header
						MPI_File_write(fh, &frameHeader_timestamp, 1, MPI_FLOAT, &status);
						break;
					default:
						global_log->error() << "[MMPLD Writer] Unsupported MMPLD version: " << _mmpldversion << endl;
						global_simulation->exit(1);
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
		for (ParticleIterator mol = particleContainer->iteratorBegin(); mol != particleContainer->iteratorEnd(); ++mol)
		{
			if(true == GetSpherePos(spherePos, &(*mol), nSphereTypeIndex) )
				MPI_File_write(fh, spherePos, 3, MPI_FLOAT, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}  // for (uint8_t nSphereTypeIndex=0; nSphereTypeIndex<_numSphereTypes; ++nSphereTypeIndex)

	// data of frame is written
	_frameCount++;

	// write seek table entry
	if (rank == 0)
	{
		MPI_Status status;
		uint64_t seektablePos =  MMPLD_HEADER_DATA_SIZE + (sizeof(uint64_t)*(_frameCount-1));
		uint64_t seekPosition;
		seekPosition = htole64(_seekTable.at(_frameCount-1) );
		MPI_File_write_at(fh, seektablePos, &seekPosition, sizeof(seekPosition), MPI_BYTE, &status);
		uint32_t frameCount = htole32(_frameCount-1);  // last frame will be ignored when simulation is aborted
		// 8: frame count position in file header
		MPI_File_write_at(fh, 8, &frameCount, sizeof(frameCount), MPI_BYTE, &status);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_close(&fh);
#endif
}

void MmpldWriter::finishOutput(ParticleContainer* /*particleContainer*/, DomainDecompBase* domainDecomp, Domain* /*domain*/)
{
	string filename = getOutputFilename();

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
		
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_seek(fh, 0, MPI_SEEK_END);
		MPI_Offset endPosition;
		MPI_File_get_position(fh, &endPosition);

		uint64_t seektablePos = MMPLD_HEADER_DATA_SIZE + (_frameCount * sizeof(uint64_t));
		uint64_t seekPosition = htole64(endPosition); /** @todo end of frame offset may not be identical to file end! */
		MPI_Status status;
		MPI_File_write_at(fh, seektablePos, &seekPosition, sizeof(seekPosition), MPI_BYTE, &status);
		uint32_t frameCount = htole32(_frameCount);  // set final number of frames
		// 8: frame count position in file header
		MPI_File_write_at(fh, 8, &frameCount, sizeof(frameCount), MPI_BYTE, &status);
		MPI_File_close(&fh);
	}else{
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
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
	_fileCount++;
	this->initOutput(particleContainer, domainDecomp, domain);
}

void MmpldWriter::PrepareWriteControl()
{
	// this method should only be executed once
	if(true == _bWriteControlPrepared)
		return;
	_bWriteControlPrepared = true;

	_startTimestep = ( _startTimestep > _simulation.getNumInitTimesteps() ) ? _startTimestep : _simulation.getNumInitTimesteps();
	if( 0 == _stopTimestep )
		_stopTimestep = _simulation.getNumTimesteps();

	if(_stopTimestep < _startTimestep) {
		global_log->warning() << "[MMPLD Writer] Invalid time interval. No frames will be recorded!" << endl;
		return;
	}

	uint64_t numTimesteps = _stopTimestep - _startTimestep;
	uint64_t numFramesTotal = numTimesteps/_writeFrequency + 1;
	if(_numFramesPerFile >= numFramesTotal || _numFramesPerFile == 0) {
		_numFramesPerFile = numFramesTotal;
	}
}

long MmpldWriter::get_data_frame_header_size() {
	long data_frame_header_size = 0;
	switch (_mmpldversion){
		case 100: /* number of particle lists (uint32_t) */
			data_frame_header_size = sizeof(float);
			break;
		case 102: /* time stamp (float) | number of particle lists (uint32_t) */
			data_frame_header_size = sizeof(float) + sizeof(uint32_t);
			break;
		default:
			global_log->error() << "[MMPLD Writer] Unsupported MMPLD version: " << _mmpldversion << endl;
			global_simulation->exit(1);
			break;
	}
	return data_frame_header_size;
}

long MmpldWriter::get_data_list_header_size() {
	long data_list_header_size = 0;
	data_list_header_size += sizeof(uint8_t); /* vertex type */
	data_list_header_size += sizeof(uint8_t); /* color type */
	if(_vertex_type == MMPLD_VERTEX_FLOAT_XYZ || _vertex_type == MMPLD_VERTEX_SHORT_XYZ) {
		data_list_header_size += sizeof(float); /* global radius */
	}
	if(_color_type == MMPLD_COLOR_NONE) {
		data_list_header_size += 4 * sizeof(uint8_t); /* global color as UINT8_RGBA */
	} else if(_color_type == MMPLD_COLOR_FLOAT_I) {
		data_list_header_size += 2 * sizeof(float); /* intensityrange min and max */
	}
	data_list_header_size += sizeof(uint64_t); /* particle count */
	return data_list_header_size;
}

long MmpldWriter::get_data_list_size(uint64_t particle_count) {
	long data_list_size = 0;
	long elemsize = 0;
	switch(_vertex_type) {
		case MMPLD_VERTEX_NONE:
			elemsize = 0;
			break;
		case MMPLD_VERTEX_FLOAT_XYZ:
			elemsize = 3 * sizeof(float);
			break;
		case MMPLD_VERTEX_FLOAT_XYZR:
			elemsize = 4 * sizeof(float);
			break;
		case MMPLD_VERTEX_SHORT_XYZ:
			elemsize = 3 * sizeof(uint16_t);
	}
	data_list_size = elemsize * particle_count;
	return data_list_size;
}

// derived classes
void MmpldWriterSimpleSphere::CalcNumSpheresPerType(ParticleContainer* particleContainer, uint64_t* numSpheresPerType)
{
	for (ParticleIterator mol = particleContainer->iteratorBegin(); mol != particleContainer->iteratorEnd(); ++mol) {
		uint8_t cid = mol->componentid();
		numSpheresPerType[cid]++;
	}
}

bool MmpldWriterSimpleSphere::GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex)
{
	uint8_t cid = mol->componentid();
	for (unsigned short d = 0; d < 3; ++d) spherePos[d] = (float)mol->r(d);
	return (cid == nSphereTypeIndex);
}


void MmpldWriterMultiSphere::CalcNumSpheresPerType(ParticleContainer* particleContainer, uint64_t* numSpheresPerType)
{
	for (ParticleIterator mol = particleContainer->iteratorBegin(); mol != particleContainer->iteratorEnd(); ++mol) {
		uint8_t cid = mol->componentid();
		uint8_t offset = _nCompSitesOffset.at(cid);
		for (uint8_t si = 0; si < _numSitesPerComp.at(cid); ++si)
			numSpheresPerType[offset+si]++;
	}
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
