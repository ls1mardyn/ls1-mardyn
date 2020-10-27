#include "io/MmspdBinWriter.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include <fstream>
#include <sstream>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;

MmspdBinWriter::MmspdBinWriter(unsigned long writeFrequency, string outputPrefix) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

MmspdBinWriter::~MmspdBinWriter(){}

void MmspdBinWriter::readXML(XMLfileUnits& xmlconfig) {
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

void MmspdBinWriter::init(ParticleContainer * /*particleContainer*/,
                          DomainDecompBase *domainDecomp, Domain *domain) {
	stringstream filenamestream;
	filenamestream << _outputPrefix;

	if(_appendTimestamp) {
		filenamestream << "-" << gettimestring();
	}
	filenamestream << ".mmspd";

	std::vector<char> filename (filenamestream.str().size()+1);
	strcpy(filename.data(),filenamestream.str().c_str());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
#endif
	ofstream mmspdfstream(filename.data(), ios::binary|ios::out);

  // format marker
  mmspdfstream << "MMSPDb";
  unsigned char padstart[2] = {0, 255};
  mmspdfstream << padstart[0] << padstart[1];

  unsigned int endi = 2018915346;
  //mmspdfstream << endi;
  mmspdfstream.write((char*)&endi,sizeof(endi));

  unsigned short version[2] = {1, 0};
  mmspdfstream.write((char*)&version,sizeof(version));
  
  unsigned char padend = 128;
  for(int i=0; i<4; ++i) mmspdfstream.write((char*)&padend,sizeof(padend));

  // header line
  bool hasIDs = true;
  mmspdfstream.write((char*)&hasIDs,sizeof(hasIDs));

  double minbox[3] = {0, 0, 0};
  double maxbox[3];
  for (unsigned short d = 0; d < 3; ++d) maxbox[d] = domain->getGlobalLength(d);
  mmspdfstream.write((char*)&minbox,sizeof(minbox));
  mmspdfstream.write((char*)&maxbox,sizeof(maxbox));

  //calculate the number of frames
  unsigned long numTimesteps = _simulation.getNumTimesteps();
  unsigned int numframes = numTimesteps/_writeFrequency;
  mmspdfstream.write((char*)&numframes,sizeof(numframes));

  unsigned int numComponents = domain-> getNumberOfComponents();
  mmspdfstream.write((char*)&numComponents,sizeof(numComponents));

  unsigned long numberMolecules = 0;
  mmspdfstream.write((char*)&numberMolecules,sizeof(numberMolecules));

  //Particle Definition
  char color;
  for(unsigned i = 0; i < domain->getNumberOfComponents() ; ++i){
	  if (i == 0){
		  //Sphere
		  mmspdfstream.write((const char*)"s", 2);

		  //4 fix constants for color
		  unsigned int fixFieldCount = 4;
		  mmspdfstream.write((char*)&fixFieldCount,sizeof(fixFieldCount));
		  //3 variables for position
		  unsigned int varFieldCount = 3;
		  mmspdfstream.write((char*)&varFieldCount,sizeof(varFieldCount));
		  //color
		  mmspdfstream.write((const char*)"cr\0b", 5);
		  color = (char)255;
		  mmspdfstream.write((char*)&color,sizeof(color));
		  mmspdfstream.write((const char*)"cg\0b", 5);
		  color = (char)0;
		  mmspdfstream.write((char*)&color,sizeof(color));
		  mmspdfstream.write((const char*)"cb\0b", 5);
		  color = (char)0;
		  mmspdfstream.write((char*)&color,sizeof(color));
		  //radius
		  mmspdfstream.write((const char*)"r\0f", 4);
		  float radius = 1.518;
		  mmspdfstream.write((char*)&radius,sizeof(radius));
		  //position
		  mmspdfstream.write((const char*)"x\0f\0y\0f\0z\0f", 12);
      }
      else if (i == 1){
		  //Sphere
		  mmspdfstream.write((const char*)"s", 2);

		  //4 fix constants for color
		  unsigned int fixFieldCount = 4;
		  mmspdfstream.write((char*)&fixFieldCount,sizeof(fixFieldCount));
		  //3 variables for position
		  unsigned int varFieldCount = 3;
		  mmspdfstream.write((char*)&varFieldCount,sizeof(varFieldCount));
		  //color
		  mmspdfstream.write((const char*)"cr\0b", 5);
		  color = (char)0;
		  mmspdfstream.write((char*)&color,sizeof(color));
		  mmspdfstream.write((const char*)"cg\0b", 5);
		  color = (char)255;
		  mmspdfstream.write((char*)&color,sizeof(color));
		  mmspdfstream.write((const char*)"cb\0b", 5);
		  color = (char)0;
		  mmspdfstream.write((char*)&color,sizeof(color));
		  //radius
		  mmspdfstream.write((const char*)"r\0f", 4);
		  float radius = 1.553;
		  mmspdfstream.write((char*)&radius,sizeof(radius));
		  //position
		  mmspdfstream.write((const char*)"x\0f\0y\0f\0z\0f", 12);
      }
      else {
    	  mmspdfstream << "**************** Error: Unspecified component!*************\n Possible reason: more than 5 components?\n";
      }
  } // end of particle definitions		
  
  mmspdfstream.close();
#ifdef ENABLE_MPI
	}
#endif
}

void MmspdBinWriter::endStep(ParticleContainer *particleContainer,
                             DomainDecompBase *domainDecomp, Domain *domain,
                             unsigned long simstep){
	if (simstep % _writeFrequency == 0) {
		stringstream filenamestream, outputstream;
		filenamestream << _outputPrefix;

		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}
		filenamestream << ".mmspd";
		
		std::vector<char> filename(filenamestream.str().size()+1);
		strcpy(filename.data(),filenamestream.str().c_str());

#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		int numprocs = domainDecomp->getNumProcs();
		unsigned long numberParticles = particleContainer->getNumberOfParticles();

		long outputsize;
		if (domain->getNumberOfComponents() == 1){
			//in byte (MPI_Byte)
			outputsize = (long)numberParticles*20; //8+4*3
		}else{
			//in byte (MPI_Byte)
			outputsize = (long)numberParticles*24; //8+4+4*3
		}

		if (rank == 0){
			outputsize += 8;
		}

		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename.data(), MPI_MODE_WRONLY|MPI_MODE_APPEND|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		for (int dest = rank+1; dest < numprocs; ++dest){
			int sendcount = 1;
			int sendtag = 0;
			MPI_Request request;
			MPI_Isend(&outputsize, sendcount, MPI_LONG, dest, sendtag, MPI_COMM_WORLD, &request);
		}
		MPI_Status status;
		long offset = 0;
		long outputsize_get;
		for (int source = 0; source < rank; ++source){
			int recvcount = 1;
			int recvtag = 0;
			MPI_Recv(&outputsize_get, recvcount, MPI_LONG, source, recvtag, MPI_COMM_WORLD, &status);
			offset += outputsize_get;
		}

		global_log->debug() << "MmspdBinWriter rank: " << rank << "; step: " << simstep << "; offset: " << offset << endl;

		MPI_File_seek(fh, offset, MPI_SEEK_END);

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0){
			unsigned long globalNumberParticles = domain->getglobalNumMolecules();
			MPI_File_write(fh, &globalNumberParticles, 1, MPI_UNSIGNED_LONG, &status);
		}
		if (domain->getNumberOfComponents() == 1){
			unsigned long molid;
			float molpos[3];
			for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
				molid = pos->getID();
				for (unsigned short d = 0; d < 3; ++d) molpos[d] = (float)pos->r(d);
				MPI_File_write(fh, &molid, 1, MPI_UNSIGNED_LONG, &status);
				MPI_File_write(fh, &molpos, 3, MPI_FLOAT, &status);
			}
		}else{
			unsigned long molid;
			unsigned int molcid;
			float molpos[3];
			for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
				molid = pos->getID();
				molcid = pos->componentid();
				for (unsigned short d = 0; d < 3; ++d) molpos[d] = (float)pos->r(d);
				MPI_File_write(fh, &molid, 1, MPI_UNSIGNED_LONG, &status);
				MPI_File_write(fh, &molcid, 1, MPI_UNSIGNED, &status);
				MPI_File_write(fh, &molpos, 3, MPI_FLOAT, &status);
			}
		}
		MPI_File_close(&fh);
#endif
	}
}

void MmspdBinWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
							Domain * /*domain*/) {}
