#include "io/MmspdWriter.h"

#ifdef ENABLE_MPI
#include <mpi.h>
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

MmspdWriter::MmspdWriter(unsigned long writeFrequency, string outputPrefix) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

MmspdWriter::~MmspdWriter(){}


void MmspdWriter::readXML(XMLfileUnits& xmlconfig) {
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

void MmspdWriter::init(ParticleContainer * /*particleContainer*/,
                       DomainDecompBase *domainDecomp, Domain *domain){
#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
#endif
	stringstream filenamestream;
	filenamestream << _outputPrefix;

	if(_appendTimestamp) {
		filenamestream << "-" << gettimestring();
	}
	filenamestream << ".mmspd";
	_filename = filenamestream.str();
	ofstream mmspdfstream(_filename.c_str(), ios::binary|ios::out);
  
  
  /* writing the header of the mmspd file, i.e. writing the BOM, the format marker (UTF-8),  the header line and defining the particle types */
  // BOM
  short int bom1,bom2,bom3;
  bom1 = 0xef;
  bom2 = 0xbb;
  bom3 = 0xbf;
  
  mmspdfstream.write(reinterpret_cast<const char*>(& bom1), 1);
  mmspdfstream.write(reinterpret_cast<const char*>(& bom2), 1);
  mmspdfstream.write(reinterpret_cast<const char*>(& bom3), 1);
  
  // format marker
  mmspdfstream << "MMSPDu 1.0" << "\n";
  // header line
  unsigned long numTimesteps = _simulation.getNumTimesteps();
  mmspdfstream << "1 " << "0 0 0 " << domain->getGlobalLength(0) <<" "<< domain->getGlobalLength(1)<< " " << domain->getGlobalLength(2) << " "
		       << numTimesteps / _writeFrequency+1    << " " << domain-> getNumberOfComponents() << " " << "0" << "\n";
		       
  
  
  /*mmspdfstream << "1 " << particleContainer->getBoundingBoxMin(0) << " " << particleContainer->getBoundingBoxMin(1) << " " 
		       << particleContainer->getBoundingBoxMin(2) << " " << particleContainer->getBoundingBoxMax(0) << " " 
		       << particleContainer->getBoundingBoxMax(1) << " " << particleContainer->getBoundingBoxMax(2) << " "
		       << _numberOfTimesteps / _writeFrequency+1    << " " << domain-> getNumberOfComponents() << " " << "0" << "\n";*/
		       
  // particle definitions every single line specifies a particular particle type
  for(unsigned i = 0; i < domain->getNumberOfComponents() ; i++){
      if (i == 0){
	mmspdfstream << "s 4 3 cr b 255 cg b 0 cb b 0 r f "; 
      }
      else if (i == 1){
	mmspdfstream << "s 4 3 cr b 0 cg b 102 cb b 0 r f "; 
      }
      else if (i == 2){
	mmspdfstream << "s 4 3 cr b 0 cg b 255 cb b 255 r f "; 
      }
      else if(i == 3){
	mmspdfstream << "s 4 3 cr b 150 cg b 0 cb b 150 r f "; 
      }
      else if (i == 4){
	mmspdfstream << "s 4 3 cr b 100 cg b 100 cb b 100 r f "; 
      }
      else {
	mmspdfstream << "**************** Error: Unspecified component!*************\n Possible reason: more than 5 components?\n"; 
      }
      mmspdfstream<< setprecision(4) << domain->getSigma(i,0)*0.7 << " x f y f z f" << "\n";
  } // end of particle definitions		
  
  mmspdfstream.close();
#ifdef ENABLE_MPI
	}
#endif
} // end init()

void MmspdWriter::endStep(ParticleContainer *particleContainer,
                          DomainDecompBase *domainDecomp, Domain *domain,
                          unsigned long simstep){
	if (simstep % _writeFrequency == 0) {
#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	int tag = 4711;
	if (rank == 0){
#endif
		ofstream mmspdfstream(_filename.c_str(), ios::out|ios::app);
		mmspdfstream << "> " << domain->getglobalNumMolecules() << "\n";
		for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				mmspdfstream << setiosflags(ios::fixed) << setw(8) << pos->getID() << setw(3)
					<< pos->componentid() << setprecision(3) << " ";
				for (unsigned short d = 0; d < 3; d++) mmspdfstream << setw(7) << pos->r(d) << " " ;
				mmspdfstream << "\n";
			}
		}
#ifdef ENABLE_MPI
		for(int fromrank = 1; fromrank < domainDecomp->getNumProcs(); fromrank++) {
			MPI_Status status_probe;
			MPI_Status status_recv;
			MPI_Probe(fromrank, tag, MPI_COMM_WORLD, &status_probe);
			int numchars;
			MPI_Get_count(&status_probe, MPI_CHAR, &numchars);
			char *recvbuff = new char[numchars];
			MPI_Recv(recvbuff, numchars, MPI_CHAR, fromrank, tag, MPI_COMM_WORLD, &status_recv);
			mmspdfstream << string(recvbuff);
			delete[] recvbuff;
		}
#endif
		mmspdfstream.close();
#ifdef ENABLE_MPI
	}
	else {
		stringstream mmspdfstream;
		for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				mmspdfstream << setiosflags(ios::fixed) << setw(8) << pos->getID() << setw(3)
					<< pos->componentid() << setprecision(3) << " ";
				for (unsigned short d = 0; d < 3; d++) mmspdfstream << setw(7) << pos->r(d) << " " ;
				mmspdfstream << "\n";
			}
		}
		
		string sendbuff;
		sendbuff = mmspdfstream.str();
		MPI_Send(sendbuff.c_str(), sendbuff.length() + 1, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
	}
#endif
  }
} // end endStep

void MmspdWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
						 Domain * /*domain*/) {}
