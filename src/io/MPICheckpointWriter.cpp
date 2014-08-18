/** \file MPICheckpointWriter.cpp
  * \brief temporary CheckpointWriter alternative using MPI-IO
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#include "io/MPICheckpointWriter.h"

#include <sstream>
#include <string>

#include "Common.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "parallel/ParticleData.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

using Log::global_log;
using namespace std;


extern Simulation* global_simulation;


MPICheckpointWriter::MPICheckpointWriter(unsigned long writeFrequency, string outputPrefix, bool incremental)
 : _outputPrefix(outputPrefix), _writeFrequency(writeFrequency), _incremental(incremental), _appendTimestamp(false)
{
	if (outputPrefix == "")
	{
		_outputPrefix=global_simulation->getOutputPrefix();
	}
	else if (outputPrefix == "default")
	{
		_appendTimestamp = true;
	}
}

MPICheckpointWriter::~MPICheckpointWriter(){}


void MPICheckpointWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "MPICheckpointWriter\twrite frequency: " << _writeFrequency << endl;
	
	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "MPICheckpointWriter\toutput prefix: " << _outputPrefix << endl;
	
	_incremental = false;
	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	//_incremental = (incremental != 0);
	if(incremental > 0) {
		_incremental = true;
	}
	global_log->info() << "MPICheckpointWriter\tincremental numbers: " << _incremental << endl;
	
	_appendTimestamp = false;
	int appendTimestamp = 0;
	//_appendTimestamp = (appendTimestamp != 0);
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "MPICheckpointWriter\tappend timestamp: " << _appendTimestamp << endl;
}

void MPICheckpointWriter::initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}

void MPICheckpointWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep, list<ChemicalPotential>* lmu ) {
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
		filenamestream << ".MPIrestart.dat";

		string filename = filenamestream.str();
		
		unsigned long nummolecules=particleContainer->getNumberOfParticles();
		int ownrank=0;
#ifdef ENABLE_MPI
		MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
		MPI_File mpifh;
		MPI_CHECK( MPI_File_open(MPI_COMM_WORLD,filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&mpifh) );
		MPI_Offset mpioffset;
		unsigned long gap=128-64;
		MPI_Status mpistat;
		unsigned long startidx;
		if(ownrank==0)
		{
			MPI_CHECK( MPI_File_write(mpifh,"MarDyn20140817",15,MPI_CHAR,&mpistat) );
			mpioffset=64-sizeof(unsigned long);
			MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&gap,1,MPI_UNSIGNED_LONG,&mpistat) );
			mpioffset+=sizeof(unsigned long);
			MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,"ICRVQD",7,MPI_CHAR,&mpistat) );
			startidx=0;
		}
		MPI_CHECK( MPI_Exscan(&nummolecules, &startidx, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD) );
		global_log->debug() << ownrank << ":\t" << nummolecules << ", " << startidx;
		MPI_Datatype mpidtParticle;
		ParticleData::setMPIType(mpidtParticle);
		MPI_Aint mpidtParticlesize=sizeof(ParticleData);	// !=MPI_Type_size
		// MPI_Type_get_extent
		/*
		ParticleData pd[2];
		MPI_Aint addr0,addr1;
		MPI_Get_address(pd,&addr0);
		MPI_Get_address(&pd[1],&addr1);
		mpidtParticlesize=addr1-addr0;
		*/
		mpioffset=64+gap+startidx*mpidtParticlesize;
		MPI_CHECK( MPI_File_set_view(mpifh,mpioffset,mpidtParticle,mpidtParticle,"external32",MPI_INFO_NULL) );
		for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			MPI_CHECK( MPI_File_write(mpifh, pos, 1, mpidtParticle, &mpistat) );
		}
		MPI_CHECK( MPI_File_close(&mpifh) );
#else
		global_log->warning() << "MPICheckpointWriter not implemented for sequential version yet" << endl;
#endif
	}
}

void MPICheckpointWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}

