#include "io/MmpldWriter.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <vector>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "/usr/include/endian.h"

// set mmpld file version. possible values: 100 or 102
#define MMPLD_FILE_VERSION 100

// set color (RGBA-Bytes) and radius (FLOAT) for each component

// component id 0
#define CID0_RADIUS 0.5
#define CID0_RED 255
#define CID0_GREEN 0
#define CID0_BLUE 0
#define CID0_ALPHA 255

// component id 1
#define CID1_RADIUS 0.553
#define CID1_RED 0
#define CID1_GREEN 0
#define CID1_BLUE 255
#define CID1_ALPHA 255

// component id 2
#define CID2_RADIUS 1.0
#define CID2_RED 0
#define CID2_GREEN 255
#define CID2_BLUE 0
#define CID2_ALPHA 255

// component id 3
#define CID3_RADIUS 1.0
#define CID3_RED 255
#define CID3_GREEN 255
#define CID3_BLUE 0
#define CID3_ALPHA 255

// component id 4
#define CID4_RADIUS 1.0
#define CID4_RED 255
#define CID4_GREEN 0
#define CID4_BLUE 255
#define CID4_ALPHA 255

// component id 5
#define CID5_RADIUS 1.0
#define CID5_RED 0
#define CID5_GREEN 255
#define CID5_BLUE 255
#define CID5_ALPHA 255


using Log::global_log;
using namespace std;

MmpldWriter::MmpldWriter(unsigned long writeFrequency, string outputPrefix) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

MmpldWriter::~MmpldWriter(){}

void MmpldWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;
	
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

//Header Information
void MmpldWriter::initOutput(ParticleContainer* /*particleContainer*/,
                           DomainDecompBase* domainDecomp, Domain* domain) {
	stringstream filenamestream;
	filenamestream << _outputPrefix;

	if(_appendTimestamp) {
		_timestampString = gettimestring();
		filenamestream << "-" << _timestampString;
	}
	filenamestream << ".mmpld";

	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
#endif
	ofstream mmpldfstream(filename, ios::binary|ios::out);

  //format marker
  uint8_t magicIdentifier[6] = {0x4D, 0x4D, 0x50, 0x4C, 0x44, 0x00};
  mmpldfstream.write((char*)magicIdentifier, sizeof(magicIdentifier));
  
  //version number
  uint16_t versionNumber;
  switch (MMPLD_FILE_VERSION){
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
  unsigned long numTimesteps = _simulation.getNumTimesteps();
  uint32_t numframes;
  uint32_t numframes_le;
  if (_writeFrequency == 0){
	  numframes = numTimesteps;
  }else{
	  numframes = (numTimesteps/_writeFrequency)+1;
  }
  _frameCount = 0;
  _numSeekEntries = numframes+1;
  numframes_le = htole32(numframes);
  mmpldfstream.write((char*)&numframes_le,sizeof(numframes_le));
  #ifdef ENABLE_MPI
	_seekTable = new uint64_t[_numSeekEntries];
  #endif
  //boundary box
  float minbox[3] = {0, 0, 0};
  float maxbox[3];
  for (unsigned short d = 0; d < 3; ++d) maxbox[d] = domain->getGlobalLength(d);
  mmpldfstream.write((char*)&minbox,sizeof(minbox));
  mmpldfstream.write((char*)&maxbox,sizeof(maxbox));
    
  //clipping box
  uint32_t numComponents = domain->getNumberOfComponents();
  float inflateRadius = CID0_RADIUS;
  if ((inflateRadius < CID1_RADIUS) && (numComponents > 1)){
	  inflateRadius = CID1_RADIUS;
  }
  if ((inflateRadius < CID2_RADIUS) && (numComponents > 2)){
	  inflateRadius = CID2_RADIUS;
  }
  if ((inflateRadius < CID3_RADIUS) && (numComponents > 3)){
	  inflateRadius = CID3_RADIUS;
  }
  if ((inflateRadius < CID4_RADIUS) && (numComponents > 4)){
	  inflateRadius = CID4_RADIUS;
  }
  if ((inflateRadius < CID5_RADIUS) && (numComponents > 5)){
	  inflateRadius = CID5_RADIUS;
  }
  
  for (unsigned short d = 0; d < 3; ++d){
	  maxbox[d] = maxbox[d] + inflateRadius;
	  minbox[d] = minbox[d] - inflateRadius;
  }
  mmpldfstream.write((char*)&minbox,sizeof(minbox));
  mmpldfstream.write((char*)&maxbox,sizeof(maxbox));
  
  //preallocate seektable
  uint64_t seekNum = 0;
  for (uint32_t i = 0; i <= numframes; ++i){
	  mmpldfstream.write((char*)&seekNum,sizeof(seekNum));
  }
  mmpldfstream.close();
#ifdef ENABLE_MPI
	}
#endif
}

void MmpldWriter::doOutput( ParticleContainer* particleContainer,
		   DomainDecompBase* domainDecomp, Domain* domain,
		   unsigned long simstep, std::list<ChemicalPotential>* /*lmu*/,
		   map<unsigned, CavityEnsemble>* /*mcav*/){
	if (simstep % _writeFrequency == 0) {
		stringstream filenamestream, outputstream;
		filenamestream << _outputPrefix;

		if(_appendTimestamp) {
			filenamestream << "-" << _timestampString;
		}
		filenamestream << ".mmpld";
		
		char filename[filenamestream.str().size()+1];
		strcpy(filename,filenamestream.str().c_str());

#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		int numprocs = domainDecomp->getNumProcs();
		unsigned long numberParticles = particleContainer->getNumberOfParticles();
		long outputsize = 0;
		uint32_t numComponents = domain->getNumberOfComponents();
		uint64_t numCompParticles[numComponents];
		for (uint32_t i = 0; i < numComponents; ++i){
			numCompParticles[i] = 0;
		}
		
		//calculate number of particles per component
		uint32_t molcid = 0;
		for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			if (numComponents > 1){
				molcid = pos->componentid();
			}else{
				molcid = 0;
			}
			if (molcid >= numComponents){
				//cout << "ERROR: Rank: " << rank << " molcid: " << molcid << endl;
				global_log->debug() << "MmpldWriter Error: Molecule ID out of range!" << endl;
				return;
			}
			numCompParticles[molcid] = numCompParticles[molcid] + 1;
		}
		
		//distribute global component particle count
		uint64_t globalNumCompParticles[numComponents];
		if (rank == 0){
			for (uint32_t i = 0; i < numComponents; ++i){
				globalNumCompParticles[i] = numCompParticles[i];
			}
			MPI_Status status;
			uint64_t numCompParticlesTmp[numComponents];
			for (int source = rank+1; source < numprocs; ++source){
				int recvcount = sizeof(numCompParticlesTmp);
				int recvtag = 1;
				MPI_Recv(numCompParticlesTmp, recvcount, MPI_BYTE, source, recvtag, MPI_COMM_WORLD, &status);
				for (uint32_t i = 0; i < numComponents; ++i){
					globalNumCompParticles[i] = globalNumCompParticles[i] + numCompParticlesTmp[i];
				}
			}
		}else{
				int dest = 0;
				int sendcount = sizeof(numCompParticles);
				int sendtag = 1;
				MPI_Request request;
				MPI_Isend(numCompParticles, sendcount, MPI_BYTE, dest, sendtag, MPI_COMM_WORLD, &request);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY|MPI_MODE_APPEND|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
		
		
		//write particle list for each component
		for (uint32_t componentIndex = 0; componentIndex < numComponents; ++componentIndex){
			//add space for particle data
			outputsize = (long)numCompParticles[componentIndex]*12;
			
			//add space for particle list header
			if (rank == 0){
				//add particle list header
				outputsize += 18;
				if (componentIndex == 0){
					
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

			global_log->debug() << "MmpldWriter rank: " << rank << "; step: " << simstep << "; component: " << componentIndex << "; offset: " << offset << endl;

			MPI_File_seek(fh, offset, MPI_SEEK_END);

			MPI_Barrier(MPI_COMM_WORLD);
			
			//write particle list header
			if (rank == 0){
				
				//write frame header if we are before the first particle list
				if (componentIndex == 0){
					//store file position for seek table
					if (_frameCount < _numSeekEntries){
						MPI_Offset entry;
						MPI_File_get_position(fh, &entry);
						_seekTable[_frameCount] = (uint64_t)entry;
					}
					_frameCount = _frameCount + 1;
					
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
					
					uint32_t frameHeader_numPLists = htole32(numComponents);
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
				
				
				//select different colors depending on componend id
				switch (componentIndex) {
					case 0:
						pListHeader_globalRadius = CID0_RADIUS;
						pListHeader_red = CID0_RED;
						pListHeader_green = CID0_GREEN;
						pListHeader_blue = CID0_BLUE;
						pListHeader_alpha = CID0_ALPHA;
						break;
						
					case 1:
						pListHeader_globalRadius = CID1_RADIUS;
						pListHeader_red = CID1_RED;
						pListHeader_green = CID1_GREEN;
						pListHeader_blue = CID1_BLUE;
						pListHeader_alpha = CID1_ALPHA;
						break;
						
					case 2:
						pListHeader_globalRadius = CID2_RADIUS;
						pListHeader_red = CID2_RED;
						pListHeader_green = CID2_GREEN;
						pListHeader_blue = CID2_BLUE;
						pListHeader_alpha = CID2_ALPHA;
						break;
						
					case 3:
						pListHeader_globalRadius = CID3_RADIUS;
						pListHeader_red = CID3_RED;
						pListHeader_green = CID3_GREEN;
						pListHeader_blue = CID3_BLUE;
						pListHeader_alpha = CID3_ALPHA;
						break;
						
					case 4:
						pListHeader_globalRadius = CID4_RADIUS;
						pListHeader_red = CID4_RED;
						pListHeader_green = CID4_GREEN;
						pListHeader_blue = CID4_BLUE;
						pListHeader_alpha = CID4_ALPHA;
						break;
						
					case 5:
						pListHeader_globalRadius = CID5_RADIUS;
						pListHeader_red = CID5_RED;
						pListHeader_green = CID5_GREEN;
						pListHeader_blue = CID5_BLUE;
						pListHeader_alpha = CID5_ALPHA;
						break;
					
					default:
						pListHeader_globalRadius = 1.0;
						pListHeader_red = 128;
						pListHeader_green = 128;
						pListHeader_blue = 128;
						pListHeader_alpha = 255;
						break;
				}
				
				//store componentParticleCount
				pListHeader_particleCount = htole64(globalNumCompParticles[componentIndex]);
				
				MPI_File_write(fh, &pListHeader_vortexType, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_colorType, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_globalRadius, 1, MPI_FLOAT, &status);
				MPI_File_write(fh, &pListHeader_red, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_green, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_blue, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_alpha, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_particleCount, 1, MPI_LONG_LONG_INT, &status);
			}
			
			float molpos[3];
			for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
				if (numComponents > 1){
					molcid = pos->componentid();
				}else{
					molcid = 0;
				}
				if (molcid == componentIndex){
					for (unsigned short d = 0; d < 3; ++d) molpos[d] = (float)pos->r(d);
					MPI_File_write(fh, molpos, 3, MPI_FLOAT, &status);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_File_close(&fh);
#endif
	}
}

void MmpldWriter::finishOutput(ParticleContainer* /*particleContainer*/, DomainDecompBase* domainDecomp, Domain* /*domain*/) {
	//fill seektable
	stringstream filenamestream;
	filenamestream << _outputPrefix;

	if(_appendTimestamp) {
		filenamestream << "-" << _timestampString;
	}
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
		_seekTable[_numSeekEntries-1] = (uint64_t)endPosition;
		MPI_File_seek(fh, 0x3C, MPI_SEEK_SET);
		uint64_t seekPosition;
		MPI_Status status;
		for (uint32_t i = 0; i < _numSeekEntries; ++i){
			seekPosition = htole64(_seekTable[i]);
			MPI_File_write(fh, &seekPosition, 1, MPI_LONG_LONG_INT, &status);
		}
		delete[] _seekTable;
		MPI_File_close(&fh);
	}else{
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_close(&fh);
	}
		
#endif
}
