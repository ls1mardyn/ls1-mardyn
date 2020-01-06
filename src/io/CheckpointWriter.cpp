#include "io/CheckpointWriter.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#include "utils/MPI_Info_object.h"
#endif

#include <sstream>
#include <string>
#include <cstring>

#include "Common.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"


using Log::global_log;
using namespace std;

void CheckpointWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	if(_writeFrequency == 0) {
		global_log->error() << "Write frequency must be a positive nonzero integer, but is " << _writeFrequency << endl;
		Simulation::exit(-1);
	}
	
	std::string checkpointType = "unknown";
	xmlconfig.getNodeValue("type", checkpointType);
	if("ASCII" == checkpointType) {
		_useBinaryFormat = false;
	}
	else if("binary" == checkpointType) {
		_useBinaryFormat = true;
	}
	else {
		global_log->error() << "Unknown CheckpointWriter type '" << checkpointType << "', expected: ASCII|binary." << endl;
		Simulation::exit(-1);
	}

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	global_log->info() << "Incremental numbers: " << _incremental << endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}else{
		_appendTimestamp = false;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void CheckpointWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                            Domain * /*domain*/) {}

void CheckpointWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                               unsigned long simstep) {
	if( simstep % _writeFrequency == 0 ) {
		stringstream filenamestream;
		filenamestream << _outputPrefix;

		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}

		if(_useBinaryFormat) {
            filenamestream << ".restart";
        }
        else { /* ASCII mode */
            filenamestream << ".restart.dat";
        }

		if(_useBinaryFormat) {
#ifdef ENABLE_MPI
			domainDecomp->assertDisjunctivity(particleContainer);
			// update global number of particles
			domain->updateglobalNumMolecules(particleContainer, domainDecomp);

			int rank = domainDecomp->getRank();
			std::string fileprefix = filenamestream.str();
			std::string filename = fileprefix + ".dat";

			MPI_File mpifh;
			MPI_Info_object mpiinfo;
			MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY|MPI_MODE_CREATE, mpiinfo, &mpifh);

			double regionLowCorner[3], regionHighCorner[3];
			for (unsigned d = 0; d < 3; d++) {
				regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
				regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
			}

			uint64_t numParticles_local = 0;
			uint64_t numParticles_exscan = 0;
			auto begin = particleContainer->regionIterator(regionLowCorner, regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			for(auto it = begin; it.isValid(); ++it)
				numParticles_local++;

			MPI_Exscan(&numParticles_local, &numParticles_exscan, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
//			cout << "[" << rank << "]: numParticles_local=" << numParticles_local << ", numParticles_exscan=" << numParticles_exscan << endl;

			uint16_t particle_data_size = 116;
			uint64_t buffer_size = 32768;
			char* write_buffer = new char[buffer_size];
			uint64_t offset = numParticles_exscan * particle_data_size;
			MPI_File_seek(mpifh, offset, MPI_SEEK_SET);
			uint64_t buffer_pos = 0;
//			auto begin = particleContainer->regionIterator(regionLowCorner, regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			for(auto it = begin; it.isValid(); ++it) {
				uint64_t pid = it->getID();
				uint32_t cid_ub = it->componentid()+1;
				double r[3];
				double v[3];
				double D[3];
				Quaternion Q = it->q();
				double q[4];
				for(auto d=0; d<3; ++d) {
					r[d] = it->r(d);
					v[d] = it->v(d);
					D[d] = it->D(d);
				}
				q[0] = Q.qw();
				q[1] = Q.qx();
				q[2] = Q.qy();
				q[3] = Q.qz();
				memcpy(&write_buffer[buffer_pos], (char*)&pid, sizeof(uint64_t)); buffer_pos += sizeof(uint64_t);
				memcpy(&write_buffer[buffer_pos], (char*)&cid_ub, sizeof(uint32_t)); buffer_pos += sizeof(uint32_t);
				memcpy(&write_buffer[buffer_pos], (char*)&r[0], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&r[1], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&r[2], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&v[0], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&v[1], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&v[2], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&q[0], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&q[1], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&q[2], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&q[3], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&D[0], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&D[1], sizeof(double)); buffer_pos += sizeof(double);
				memcpy(&write_buffer[buffer_pos], (char*)&D[2], sizeof(double)); buffer_pos += sizeof(double);

				if(buffer_pos > buffer_size - particle_data_size) {
					MPI_Status status;
					MPI_File_write(mpifh, write_buffer, buffer_pos, MPI_BYTE, &status);
					buffer_pos = 0;
				}
			}
			MPI_Status status;
			MPI_File_write(mpifh, write_buffer, buffer_pos, MPI_BYTE, &status);

			delete[] write_buffer;
			MPI_File_close(&mpifh);

			// writer Header
			double currentTime = global_simulation->getSimulationTime();
			filename = fileprefix + ".header.xml";
			domain->writeCheckpointHeaderXML(filename, particleContainer, domainDecomp, currentTime);
#else
			string filename = filenamestream.str();
			domain->writeCheckpoint(filename, particleContainer, domainDecomp, _simulation.getSimulationTime(), _useBinaryFormat);
#endif
		}
		else { /* ASCII mode */
			string filename = filenamestream.str();
			domain->writeCheckpoint(filename, particleContainer, domainDecomp, _simulation.getSimulationTime(), _useBinaryFormat);
		}
	}
}

void CheckpointWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
							  Domain * /*domain*/) {
}
