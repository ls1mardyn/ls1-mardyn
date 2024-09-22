#include "io/MmpldWriter.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#ifdef _SX
#include <byteswap.h>
#define htole16(X) bswap2((X))
#define htole32(X) bswap4((X))
#define htole64(X) bswap8((X))
#else
#include <endian.h>
#endif

#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <iostream>
#include <cstdint>

#include "Common.h"
#include "Domain.h"
#include "Simulation.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "utils/FileUtils.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"


// default version to use for mmpld format writing. possible values: 100 or 102
#define MMPLD_DEFAULT_VERSION 100
#define MMPLD_HEADER_DATA_SIZE 60
#define MMPLD_SEEK_TABLE_OFFSET MMPLD_HEADER_DATA_SIZE

std::string MmpldWriter::getOutputFilename() {
	std::stringstream filenamestream;
	filenamestream << _outputPrefix << "_" << fill_width('0', 4) << _fileCount << ".mmpld";
	return filenamestream.str();
}

MmpldWriter::MmpldWriter() :
		_startTimestep(0), _writeFrequency(1000), _stopTimestep(std::numeric_limits<uint64_t>::max()), _writeBufferSize(32768), _outputPrefix("unknown"),
		_bInitSphereData(ISD_READ_FROM_XML), _bWriteControlPrepared(false),
		_fileCount(1), _numFramesPerFile(0), _mmpldversion(MMPLD_DEFAULT_VERSION), _vertex_type(MMPLD_VERTEX_FLOAT_XYZ), _color_type(MMPLD_COLOR_NONE)
{}

MmpldWriter::MmpldWriter(uint64_t startTimestep, uint64_t writeFrequency, uint64_t stopTimestep, uint64_t numFramesPerFile,
		std::string outputPrefix)
		:	_startTimestep(startTimestep), _writeFrequency(writeFrequency), _stopTimestep(stopTimestep), _writeBufferSize(32768),
		_outputPrefix(outputPrefix), _bInitSphereData(ISD_READ_FROM_XML), _bWriteControlPrepared(false),
		_fileCount(1),_numFramesPerFile(numFramesPerFile),  _vertex_type(MMPLD_VERTEX_FLOAT_XYZ),
		_color_type(MMPLD_COLOR_NONE)
{
	if (0 == _writeFrequency) {
		MARDYN_EXIT(-1);
	}
}

void MmpldWriter::readXML(XMLfileUnits& xmlconfig)
{
	// color type
	_color_type = MMPLD_COLOR_NONE;
	uint32_t ctype = 0;
	xmlconfig.getNodeValue("@ctype", ctype);
	_color_type = static_cast<MMPLD_Color_type>(ctype);

	// write control
	xmlconfig.getNodeValue("writecontrol/start", _startTimestep);
	xmlconfig.getNodeValue("writecontrol/writefrequency", _writeFrequency);
	xmlconfig.getNodeValue("writecontrol/stop", _stopTimestep);
	xmlconfig.getNodeValue("writecontrol/framesperfile", _numFramesPerFile);
	xmlconfig.getNodeValue("writecontrol/writeBufferSize", _writeBufferSize);
	Log::global_log->info() << "[MMPLD Writer] Start sampling from simstep: " << _startTimestep << std::endl;
	Log::global_log->info() << "[MMPLD Writer] Write with frequency: " << _writeFrequency << std::endl;
	Log::global_log->info() << "[MMPLD Writer] Stop sampling at simstep: " << _stopTimestep << std::endl;
	Log::global_log->info() << "[MMPLD Writer] Split files every " << _numFramesPerFile << "th frame."<< std::endl;
	Log::global_log->info() << "[MMPLD Writer] Write buffer size: " << _writeBufferSize << " Byte" << std::endl;

	int mmpldversion = MMPLD_DEFAULT_VERSION;
	xmlconfig.getNodeValue("mmpldversion", mmpldversion);
	_mmpldversion = mmpldversion;
	switch(_mmpldversion) {
		case 100:
		case 102:
			break;
		default:
			Log::global_log->error() << "Unsupported MMPLD version:" << _mmpldversion << std::endl;
			MARDYN_EXIT(1);
			break;
	}
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[MMPLD Writer] Output prefix: " << _outputPrefix << std::endl;

	// sphere params: radius, colors
	uint32_t numSites = 0;
	XMLfile::Query query = xmlconfig.query("spheres/site");
	numSites = query.card();
	Log::global_log->info() << "[MMPLD Writer] Number of sites: " << numSites << std::endl;
	if(numSites < 1) {
		Log::global_log->fatal() << "[MMPLD Writer] No site parameters specified." << std::endl;
		MARDYN_EXIT(48973);
	}
	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator outputSiteIter;
	for( outputSiteIter = query.begin(); outputSiteIter; outputSiteIter++ )
	{
		xmlconfig.changecurrentnode( outputSiteIter );
		float radius = 0.;
		int r = 0, g = 0, b = 0, alpha = 0;
		float intensity_min = 0., intensity_max = 0.;
		xmlconfig.getNodeValue("radius", radius);
		xmlconfig.getNodeValue("color/r", r);
		xmlconfig.getNodeValue("color/g", g);
		xmlconfig.getNodeValue("color/b", b);
		xmlconfig.getNodeValue("color/alpha", alpha);
		xmlconfig.getNodeValue("intensity/min", intensity_min);
		xmlconfig.getNodeValue("intensity/max", intensity_max);

		_global_radius.push_back(radius);
		std::array<uint8_t, 4> rgba = { { static_cast<uint8_t>(r), static_cast<uint8_t>(g), static_cast<uint8_t>(b), static_cast<uint8_t>(alpha)} };
		_global_rgba.push_back(rgba);
		std::array<float, 2> intRange = { {intensity_min, intensity_max} };
		_global_intensity_range.push_back(intRange);
	}

	if(xmlconfig.changecurrentnode("mpi_info")) {
#ifdef ENABLE_MPI
		Log::global_log->info() << "[MMPLD Writer] Setting MPI info object for IO" << std::endl;
		_mpiinfo.readXML(xmlconfig);
#else
		Log::global_log->info() << "[MMPLD Writer] mpi_info only used in parallel/MPI version" << std::endl;
#endif
		xmlconfig.changecurrentnode("..");
	}
}

//Header Information
void MmpldWriter::init(ParticleContainer *particleContainer,
						DomainDecompBase *domainDecomp, Domain *domain)
{
	if ( (htole32(1) != 1) || (htole64(1.0) != 1.0) ) {
		Log::global_log->error() << "[MMPLD Writer] The MMPLD Writer currently only supports running on little endian systems." << std::endl;
		MARDYN_EXIT(1);
	}

	// only executed once
	this->PrepareWriteControl();

	_frameCount = 0;

	// number of components / sites
	std::vector<Component> *components = global_simulation->getEnsemble()->getComponents();
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
		Log::global_log->debug() << "[MMPLD Writer] Component[" << cid << "] numSites=" << numSites << " offset=" << unsigned(_nCompSitesOffset.at(cid)) << std::endl;
		_numSitesTotal += numSites;
	}
	Log::global_log->debug() << "[MMPLD Writer] Total number of sites taken into account: " << unsigned(_numSitesTotal) << std::endl;

	// init radius and color of spheres
	this->InitSphereData();
	this->SetNumSphereTypes();

	std::string filename = getOutputFilename();
	uint8_t magicIdentifier[6] = {0x4D, 0x4D, 0x50, 0x4C, 0x44, 0x00}; // format marker
	uint16_t mmpldversion_le = htole16(_mmpldversion);
	uint32_t numframes = _numFramesPerFile; // number of frames
	uint32_t numframes_le = htole32(numframes);
	_numSeekEntries = numframes + 1; /* need additional seek entry for end of file offset marker */
	_seekTable.resize(_numSeekEntries);
	_seekTable.at(0) = MMPLD_HEADER_DATA_SIZE + get_seekTable_size();

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
#endif
	std::ofstream mmpldfstream(filename.c_str(), std::ios::binary|std::ios::out);
	mmpldfstream.write((char*)magicIdentifier, sizeof(magicIdentifier));
	mmpldfstream.write((char*)&mmpldversion_le, sizeof(mmpldversion_le));
	mmpldfstream.write((char*)&numframes_le,sizeof(numframes_le));
	Log::global_log->debug() << "[MMPLD Writer] Writing bounding box data." << std::endl;
	float minbox[3] = {0, 0, 0};
	float maxbox[3];
	for (unsigned short d = 0; d < 3; ++d) {
		maxbox[d] = domain->getGlobalLength(d);
	}
	mmpldfstream.write((char*)minbox, sizeof(minbox));
	mmpldfstream.write((char*)maxbox, sizeof(maxbox));
	Log::global_log->debug() << "[MMPLD Writer] Writing clipping box data." << std::endl;
	float inflateRadius = 0;
	for(auto radius : _global_radius) {
		if(inflateRadius < radius ) {
			inflateRadius = radius;
		}
	}
	for (unsigned short d = 0; d < 3; ++d){
		maxbox[d] = maxbox[d] + inflateRadius;
		minbox[d] = minbox[d] - inflateRadius;
	}
	mmpldfstream.write((char*)minbox, sizeof(minbox));
	mmpldfstream.write((char*)maxbox, sizeof(maxbox));
	Log::global_log->debug() << "[MMPLD Writer] Preallocating " << _numSeekEntries << " seek table entries for frames" << std::endl;
	for (uint32_t i = 0; i < _numSeekEntries; ++i) {
		uint64_t offset_le = htole64(_seekTable.at(i));
		mmpldfstream.write((char*) &offset_le, sizeof(offset_le));
	}
	mmpldfstream.close();
#ifdef ENABLE_MPI
	}
#endif
}


void MmpldWriter::write_frame(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp) {
	std::string filename = getOutputFilename();

	// calculate local number of spheres per component|siteType
	std::vector<uint64_t> numSpheresPerType(_numSphereTypes);
	this->CalcNumSpheresPerType(particleContainer, numSpheresPerType.data());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();

	MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY|MPI_MODE_CREATE, _mpiinfo, &_mpifh);

	//distribute global component particle count and offset counts for distrubted write
	std::vector<uint64_t> globalNumCompSpheres(_numSphereTypes);
	std::vector<uint64_t> exscanNumCompSpheres(_numSphereTypes);
	MPI_Exscan(numSpheresPerType.data(), exscanNumCompSpheres.data(), _numSphereTypes, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
	for(int i = 0; i < _numSphereTypes; ++i) {
		globalNumCompSpheres[i] = exscanNumCompSpheres[i] + numSpheresPerType[i];
	}
	int lastrank = domainDecomp->getNumProcs() - 1;
	MPI_Bcast(globalNumCompSpheres.data(), _numSphereTypes, MPI_UINT64_T, lastrank, MPI_COMM_WORLD);

	if(rank == 0) {
		MPI_File_seek(_mpifh, _seekTable.at(_frameCount), MPI_SEEK_SET);
		write_frame_header(_numSphereTypes);
	}

	/* positions of data lists relative to frame begin */
	std::vector<uint64_t> dataListBeginOffsets(_numSphereTypes);
	dataListBeginOffsets[0] = get_data_frame_header_size();
	for(int i = 1; i < _numSphereTypes; ++i) {
		dataListBeginOffsets[i] = dataListBeginOffsets[i-1] + get_data_list_header_size() + get_data_list_size(globalNumCompSpheres[i-1]);
	}
	/* calculate write positions inside data lists for this process */
	std::vector<uint64_t> dataListWriteOffsets(_numSphereTypes);
	for(int i = 0; i < _numSphereTypes; ++i) {
		dataListWriteOffsets[i] = get_data_list_header_size() + get_data_list_size(exscanNumCompSpheres[i]);
	}

	char* writeBuffer = new char[_writeBufferSize];
	long buffer_pos = 0;
	/* write particle list for each component|site (sphere type)`*/
	for (uint8_t sphereTypeId = 0; sphereTypeId < _numSphereTypes; ++sphereTypeId){
		//write particle list header
		if(rank == 0) {
			long offset = _seekTable.at(_frameCount) + dataListBeginOffsets[sphereTypeId];
			MPI_File_seek(_mpifh, offset, MPI_SEEK_SET);
			write_particle_list_header(globalNumCompSpheres[sphereTypeId], sphereTypeId);
		}
		long offset = _seekTable.at(_frameCount) + dataListBeginOffsets[sphereTypeId] + dataListWriteOffsets[sphereTypeId];
		MPI_File_seek(_mpifh, offset, MPI_SEEK_SET);
		buffer_pos = 0;
		for (auto moleculeIter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); moleculeIter.isValid(); ++moleculeIter) {
			if(true == GetSpherePos(reinterpret_cast<float*>(&writeBuffer[buffer_pos]), &(*moleculeIter), sphereTypeId)) {
				buffer_pos += get_particle_data_size();
				if(buffer_pos > _writeBufferSize - get_particle_data_size()) {
					MPI_Status status;
					MPI_File_write(_mpifh, writeBuffer, buffer_pos, MPI_BYTE, &status);
					buffer_pos = 0;
				}
			}
		}
		MPI_Status status;
		MPI_File_write(_mpifh, writeBuffer, buffer_pos, MPI_BYTE, &status);
	}
	delete[] writeBuffer;
	// data of frame is written
	_frameCount++;
	uint64_t frame_offset = dataListBeginOffsets.back() + get_data_list_header_size() + get_data_list_size(globalNumCompSpheres.back());
	_seekTable.at(_frameCount) = _seekTable.at(_frameCount - 1) + frame_offset;
	// write seek table entry and update frame count
	if (rank == 0) {
		writeSeekTableEntry(_frameCount, _seekTable.at(_frameCount));
		MPI_Status status;
		uint32_t frameCount = htole32(_frameCount);
		// 8: frame count position in file header
		MPI_File_write_at(_mpifh, 8, &frameCount, sizeof(frameCount), MPI_BYTE, &status);
	}
	MPI_File_close(&_mpifh);
#endif
}

void MmpldWriter::endStep(ParticleContainer *particleContainer,
							DomainDecompBase *domainDecomp, Domain *domain,
							unsigned long simstep)
{
	if((simstep < _startTimestep) || (simstep > _stopTimestep) || (0 != ((simstep - _startTimestep) % _writeFrequency)) ) {
		return;
	}
	if(_frameCount == _numFramesPerFile) {
		MultiFileApproachReset(particleContainer, domainDecomp, domain);  // begin new file
	}

	std::string filename = getOutputFilename();
	Log::global_log->debug() << "[MMPLD Writer] Writing MMPLD frame " << _frameCount << " for simstep " << simstep << " to file " << filename << std::endl;
	write_frame(particleContainer, domainDecomp);
}

void MmpldWriter::InitSphereData()
{
	if(_bInitSphereData == ISD_READ_FROM_XML)
		return;
	else if(_bInitSphereData == ISD_USE_DEFAULT)
	{
		for(uint8_t i=0; i<6; i++) {
			_global_radius.push_back(0.5);
		}
		//                            R    G    B  alpha
		std::array<uint8_t, 4> red = { 255, 0, 0, 255 };
		_global_rgba.push_back(red);
		std::array<uint8_t, 4> lightblue = { 0, 205, 255, 255 };
		_global_rgba.push_back(lightblue);
		std::array<uint8_t, 4> blue = { 255, 0, 255, 255 };
		_global_rgba.push_back(blue);
		std::array<uint8_t, 4> green = { 0, 155, 0, 255 };
		_global_rgba.push_back(green);
		std::array<uint8_t, 4> purple = { 105, 0, 205, 255 };
		_global_rgba.push_back(purple);
		std::array<uint8_t, 4> orange = { 255, 125, 0, 255 };
		_global_rgba.push_back(orange);
		return;
	}
}

void MmpldWriter::MultiFileApproachReset(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain)
{
	this->finish(particleContainer, domainDecomp, domain);
	_fileCount++;
	this->init(particleContainer, domainDecomp, domain);
}

void MmpldWriter::PrepareWriteControl()
{
	// this method should only be executed once
	if(true == _bWriteControlPrepared)
		return;
	_bWriteControlPrepared = true;

	_startTimestep = std::max(_startTimestep, _simulation.getNumInitTimesteps());

	_stopTimestep = std::min(_stopTimestep, _simulation.getNumTimesteps());

	Log::global_log->info() << "[MMPLD Writer] Setting start:stop to " << _startTimestep << ":" << _stopTimestep << std::endl;

	if(_stopTimestep < _startTimestep) {
		Log::global_log->warning() << "[MMPLD Writer] Invalid time interval. No frames will be recorded!" << std::endl;
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
			data_frame_header_size = sizeof(uint32_t);
			break;
		case 102: /* time stamp (float) | number of particle lists (uint32_t) */
			data_frame_header_size = sizeof(float) + sizeof(uint32_t);
			break;
		default:
			Log::global_log->error() << "[MMPLD Writer] Unsupported MMPLD version: " << _mmpldversion << std::endl;
			MARDYN_EXIT(1);
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

long MmpldWriter::get_particle_data_size() {
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
	switch(_color_type) {
		case MMPLD_COLOR_NONE:
			break;
		case MMPLD_COLOR_UINT8_RGB:
			elemsize += 3 * sizeof(uint8_t);
			break;
		case MMPLD_COLOR_UINT8_RGBA:
			elemsize += 4 * sizeof(uint8_t);
			break;
		case MMPLD_COLOR_FLOAT_I:
			elemsize += sizeof(float);
			break;
		case MMPLD_COLOR_FLOAT_RGB:
			elemsize += 3 * sizeof(float);
			break;
		case MMPLD_COLOR_FLOAT_RGBA:
			elemsize += 4 * sizeof(float);
			break;
	}
	return elemsize;
}


long MmpldWriter::get_data_list_size(uint64_t particle_count) {
	return particle_count * get_particle_data_size();
}


void MmpldWriter::write_frame_header(uint32_t num_data_lists) {
#ifdef ENABLE_MPI
	MPI_Status status;
	if (_mmpldversion == 102){
		float frameHeader_timestamp = _simulation.getSimulationTime();
		MPI_File_write(_mpifh, &frameHeader_timestamp, 1, MPI_FLOAT, &status);
	}
	uint32_t num_data_lists_le = htole32(num_data_lists);
	MPI_File_write(_mpifh, &num_data_lists_le, sizeof(num_data_lists_le), MPI_BYTE, &status);
#else
	/** @todo need to implement serial version */
#endif
}

long MmpldWriter::get_seekTable_size(){
	return _numSeekEntries * sizeof(uint64_t);
}

void MmpldWriter::writeSeekTableEntry(int id, uint64_t offset) {
#ifdef ENABLE_MPI
	MPI_Offset seekpos = MMPLD_SEEK_TABLE_OFFSET + id*sizeof(uint64_t);
	uint64_t offset_le = htole64(offset);
	MPI_Status status;
	MPI_File_write_at(_mpifh, seekpos, &offset_le, 8, MPI_BYTE, &status);
#endif
}

void MmpldWriter::write_particle_list_header(uint64_t particle_count, int sphereId) {
#ifdef ENABLE_MPI
	MPI_Status status;
	MPI_File_write(_mpifh, &_vertex_type,  1, MPI_BYTE, &status);
	MPI_File_write(_mpifh, &_color_type,   1, MPI_BYTE, &status);
	if(_vertex_type == MMPLD_VERTEX_FLOAT_XYZ || _vertex_type == MMPLD_VERTEX_SHORT_XYZ) {
		MPI_File_write(_mpifh, &_global_radius[sphereId], 4, MPI_BYTE, &status);
	}
	if(_color_type == MMPLD_COLOR_NONE) {
		MPI_File_write(_mpifh, &_global_rgba[sphereId],   4, MPI_BYTE, &status);
	} else if(_color_type == MMPLD_COLOR_FLOAT_I) {
		MPI_File_write(_mpifh, &_global_intensity_range[sphereId],   8, MPI_BYTE, &status);
	}
	uint64_t particle_count_le = htole64(particle_count);
	MPI_File_write(_mpifh, &particle_count_le, sizeof(particle_count_le), MPI_BYTE, &status);
#else
	/** @todo need to implement serial version */
#endif
}


// derived classes
void MmpldWriterSimpleSphere::CalcNumSpheresPerType(ParticleContainer* particleContainer, uint64_t* numSpheresPerType)
{
	for (auto mol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
		uint8_t cid = mol->componentid();
		numSpheresPerType[cid]++;
	}
}

bool MmpldWriterSimpleSphere::GetSpherePos(float *spherePos, Molecule* mol, uint8_t& nSphereTypeIndex)
{
	uint8_t cid = mol->componentid();
	for (unsigned short d = 0; d < 3; ++d) spherePos[d] = (float)mol->r(d);
	if(MMPLD_COLOR_FLOAT_RGB == _color_type)  // color hack
		for (unsigned short d = 0; d < 3; ++d) spherePos[d+3] = (float)mol->v(d);
	return (cid == nSphereTypeIndex);
}


void MmpldWriterMultiSphere::CalcNumSpheresPerType(ParticleContainer* particleContainer, uint64_t* numSpheresPerType)
{
	for (auto mol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
		uint8_t cid = mol->componentid();
		uint8_t offset = _nCompSitesOffset.at(cid);
		for (uint8_t si = 0; si < _numSitesPerComp.at(cid); ++si)
			numSpheresPerType[offset+si]++;
	}
}

bool MmpldWriterMultiSphere::GetSpherePos(float *spherePos, Molecule* mol, uint8_t& nSphereTypeIndex)
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
