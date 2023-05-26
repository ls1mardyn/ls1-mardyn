#include "io/VISWriter.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "Common.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif


VISWriter::VISWriter(unsigned long writeFrequency, std::string outputPrefix) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;
	_wroteVIS = false;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

VISWriter::~VISWriter(){}

void VISWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "Output prefix: " << _outputPrefix << std::endl;
	
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	Log::global_log->info() << "Append timestamp: " << _appendTimestamp << std::endl;
}

void VISWriter::init(ParticleContainer * /*particleContainer*/,
                     DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/) {
    std::string filename = _outputPrefix + ".vis";
	std::ofstream fileout(filename.c_str(), std::ios::out);
	fileout.close();
}

void VISWriter::endStep(ParticleContainer *particleContainer,
                        DomainDecompBase *domainDecomp, Domain * /*domain*/,
                        unsigned long simstep) {
	if (simstep % _writeFrequency == 0) {
		std::stringstream filenamestream, outputstream;
		filenamestream << _outputPrefix;

		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}
		filenamestream << ".vis";
		
		std::vector<char> filename(filenamestream.str().size()+1);
		strcpy(filename.data(),filenamestream.str().c_str());

#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		int numprocs = domainDecomp->getNumProcs();
		if (rank== 0){
#endif
			if (!_wroteVIS){
				outputstream << "      id t          x          y          z     q0     q1     q2     q3        c\n";
				_wroteVIS = true;
			}
			else
				outputstream << "#" << std::endl;
#ifdef ENABLE_MPI
		}
#endif

		// originally VIS files had a fixed width of 8 (and no t), here I use 12 (with 2 for t)
		//ostrm << "t           x           y           z          q0          q1          q2          q3" << std::endl;
		for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				outputstream << setiosflags(std::ios::fixed) << std::setw(8) << pos->getID() << std::setw(2)
				            << pos->componentid() << std::setprecision(3);
				for (unsigned short d = 0; d < 3; d++) outputstream << std::setw(11) << pos->r(d);
				outputstream << std::setprecision(3) << std::setw(7) << pos->q().qw() << std::setw(7) << pos->q().qx()
				            << std::setw(7) << pos->q().qy()<< std::setw(7) << pos->q().qz()
				            << std::setw(9) << std::right << 0 << "\n";
			}
		}
		long outputsize = outputstream.str().size();

		std::vector<char> output(outputsize+1);
		strcpy(output.data(),outputstream.str().c_str());
#ifdef ENABLE_MPI
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename.data(), MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh);

		for (int dest = rank+1; dest < numprocs; dest++){
			int sendcount = 1;
		    int sendtag = 0;
		    MPI_Request request;
		    MPI_Isend(&outputsize, sendcount, MPI_LONG, dest, sendtag, MPI_COMM_WORLD, &request);
		}
		MPI_Status status;
		long offset = 0;
		long outputsize_get;
		for (int source = 0; source < rank; source++){
			int recvcount = 1;
		    int recvtag = 0;
		    MPI_Recv(&outputsize_get, recvcount, MPI_LONG, source, recvtag, MPI_COMM_WORLD, &status);
		    offset += outputsize_get;
		}

		MPI_File_seek(fh, offset, MPI_SEEK_END);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_File_write(fh, output.data(), outputsize, MPI_CHAR, &status);
		MPI_File_close(&fh);
#else
		std::ofstream fileout(filename.data(), std::ios::out|std::ios::app);
		fileout << output.data();
		fileout.close();
#endif
	}
}

void VISWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
					   Domain * /*domain*/) {}
