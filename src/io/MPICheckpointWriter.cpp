/** \file MPICheckpointWriter.cpp
  * \brief temporary CheckpointWriter alternative using MPI-IO
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#include "io/MPICheckpointWriter.h"

#include <sstream>
#include <fstream>
#include <string>
#include <cstring>

#include <sys/time.h>	// gettimeofday()


#include "Common.h"
#include "Domain.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "parallel/DomainDecompBase.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#include "parallel/ParticleData.h"
#endif


extern Simulation* global_simulation;

const char MPICheckpointWriter::_magicVersion[] = "MarDyn20150211trunk";
//    int32_t
const int MPICheckpointWriter::_endiannesstest = 0x0a0b0c0d;

MPICheckpointWriter::MPICheckpointWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental, std::string datarep)
 : _outputPrefix(outputPrefix), _writeFrequency(writeFrequency), _incremental(incremental), _appendTimestamp(false), _datarep(datarep)
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

void MPICheckpointWriter::readXML(XMLfileUnits& xmlconfig)
{
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "[MPICheckpointWriter]\twrite frequency: " << _writeFrequency << std::endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[MPICheckpointWriter]\toutput prefix: " << _outputPrefix << std::endl;

	_incremental = false;
	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	//_incremental = (incremental != 0);
	if(incremental > 0) {
		_incremental = true;
		Log::global_log->info() << "[MPICheckpointWriter]\tusing incremental numbers in file names" << std::endl;
	}

	_appendTimestamp = false;
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	//_appendTimestamp = (appendTimestamp != 0);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
		Log::global_log->info() << "[MPICheckpointWriter]\tappend timestamp to file names" << std::endl;
	}

	_datarep = "";	// -> NULL
	//_datarep = "external32";	// "native", "internal", "external32"
	xmlconfig.getNodeValue("datarep", _datarep);
	if(!_datarep.empty())
		Log::global_log->info() << "[MPICheckpointWriter]\tdata representation: " << _datarep << std::endl;

	_measureTime = false;
	int measureTime = 0;
	xmlconfig.getNodeValue("measureTime", measureTime);
	//_measureTime = (measureTime != 0);
	if(measureTime > 0) {
		_measureTime = true;
		Log::global_log->info() << "[MPICheckpointWriter]\texecution wall time will be measured" << std::endl;
	}

	if(xmlconfig.changecurrentnode("mpi_info")) {
#ifdef ENABLE_MPI
		Log::global_log->info() << "[MPICheckpointWriter] Setting MPI info object for IO" << std::endl;
		_mpiinfo.readXML(xmlconfig);
#else
		Log::global_log->info() << "[MPICheckpointWriter] mpi_info only used in parallel/MPI version" << std::endl;
#endif
		xmlconfig.changecurrentnode("..");
	}

	_particlesbuffersize = 0;
	xmlconfig.getNodeValue("ParticlesBufferSizeMPI", _particlesbuffersize);
	if(_particlesbuffersize)
	{
#ifdef ENABLE_MPI
		Log::global_log->info() << "[MPICheckpointWriter]\tparticles buffer size: " << _particlesbuffersize << std::endl;
#else
		Log::global_log->info() << "[MPICheckpointWriter]\tparticles buffer size (" << _particlesbuffersize << ") only used in parallel/MPI version" << std::endl;
#endif
	}
}

void MPICheckpointWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                               Domain * /*domain*/)
{
		if(_incremental && _appendTimestamp)
		{	// use the same timestamp for all increments: add it to the outputPrefix
			std::stringstream outputPrefixstream;
			outputPrefixstream << _outputPrefix;
			char fmt[] = "%Y%m%dT%H%M%S"; // must have fixed size format for all time values/processes
			char timestring[256];
#ifdef ENABLE_MPI
			int count = 0;
			count = gettimestr(fmt, timestring, sizeof(timestring)/sizeof(timestring[0]));
			MPI_CHECK(MPI_Bcast(timestring, count, MPI_CHAR, 0, MPI_COMM_WORLD));
#else
			gettimestr(fmt, timestring, sizeof(timestring)/sizeof(timestring[0]));
#endif
			outputPrefixstream << "_" << std::string(timestring);
			_outputPrefix = outputPrefixstream.str();
		}
}

void MPICheckpointWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                                  unsigned long simstep) {
#ifdef ENABLE_MPI
	const char *mpidatarep = NULL;
	if (!_datarep.empty()) mpidatarep=_datarep.c_str();
#endif

	if( simstep % _writeFrequency == 0 ) {
		std::stringstream filenamestream;
		filenamestream << _outputPrefix;

		if(_incremental)
		{	/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */

			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		else if(_appendTimestamp)
		{	//  different timestamps if incremental is not used (no file numbers)
			char fmt[] = "%Y%m%dT%H%M%S"; // must have fixed size format for all time values/processes
			char timestring[256];
#ifdef ENABLE_MPI
			int count = 0;
			count = gettimestr(fmt, timestring, sizeof(timestring)/sizeof(timestring[0]));
			MPI_CHECK(MPI_Bcast(timestring, count, MPI_CHAR, 0, MPI_COMM_WORLD));
#else
			gettimestr(fmt, timestring, sizeof(timestring)/sizeof(timestring[0]));
#endif
			filenamestream << "_" << std::string(timestring);
		}
		filenamestream << ".MPIrestart.dat";

		std::string filename = filenamestream.str();
		Log::global_log->info() << "[MPICheckpointWriter]\tfilename: " << filename << std::endl;

		unsigned long numParticles_global = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);
		unsigned long numParticles = particleContainer->getNumberOfParticles();	// local
		unsigned long numbb{1ul};
#ifdef ENABLE_MPI
		Log::global_log->info() << "[MPICheckpointWriter]\tnumber of particles: " << numParticles_global
		                   << "\t(*" << sizeof(ParticleData) << "=" << numParticles_global*sizeof(ParticleData) << " Bytes in memory)"
				   << std::endl;
		//global_log->set_mpi_output_all()
		int num_procs;
		MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &num_procs) );
		unsigned long gap=7+3+sizeof(unsigned long)+num_procs*(6*sizeof(double)+2*sizeof(unsigned long));
		int ownrank;
		MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
		double mpistarttime=0;	// =0 to prevent Jenkins/gcc complaining about uninitialized mpistarttime [-Werror=uninitialized]
		if(_measureTime)
		{
			//if(ownrank==0) Log::global_log->debug() << "MPICheckpointWriter (" << filename << ")\tstart measuring time" << std::endl;
			MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
			mpistarttime=MPI_Wtime();
			// global_simulation->timers()->start("MPI_CHECKPOINT_WRITER_INPUT");	// should use Timer instead
		}
		MPI_File mpifh;
		// arg 2 type cast due to old MPI (<=V2) implementations (should be const char* now)
		// Why does an explicit C cast  (char*)  not work?  -> should be interpreted like a const_cast in the first place
		// (see e.g. http://en.cppreference.com/w/cpp/language/explicit_cast)
		MPI_CHECK( MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY|MPI_MODE_CREATE, _mpiinfo, &mpifh) );
		//MPI_CHECK( MPI_File_preallocate( mpifh, mpifilesize); )	// might ensure that data can be written, but might be slow
		//MPI_Aint mpidtParticleDte; MPI_CHECK( MPI_File_get_type_extent(mpifh,mpidtParticleD,&mpidtParticleDte) );
		// arg 5 type cast due to old MPI (<=V2) implementations (should be const char* now)
		MPI_CHECK( MPI_File_set_view(mpifh, 0, MPI_BYTE, MPI_BYTE, const_cast<char*>(mpidatarep), _mpiinfo) );
		MPI_Offset mpioffset=0;
		MPI_Status mpistat;
		unsigned long startidx;
		if(ownrank==0)
		{	// the first part of header will be written by rank 0 only
			//MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,_magicVersion,strlen(_magicVersion)+1,MPI_CHAR,&mpistat) );
			// arg 2 type cast due to old MPI (<=V2) implementations (should be const void* now)
			MPI_CHECK( MPI_File_write(mpifh, (void*)const_cast<char*>(_magicVersion), strlen(_magicVersion)+1, MPI_CHAR, &mpistat) );
			mpioffset=64-sizeof(unsigned long)-sizeof(int);
			//MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&_endiannesstest,1,MPI_INT,&mpistat) );
			MPI_CHECK( MPI_File_seek(mpifh,mpioffset,MPI_SEEK_SET) );
			// arg 2 type cast due to old MPI (<=V2) implementations (should be const void* now)
			MPI_CHECK( MPI_File_write(mpifh, (void*)const_cast<int*>(&_endiannesstest), 1, MPI_INT, &mpistat) );
			//mpioffset=64-sizeof(unsigned long);
			//MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&gap,1,MPI_UNSIGNED_LONG,&mpistat) );
			//MPI_CHECK( MPI_File_seek(mpifh,mpioffset,MPI_SEEK_SET) );
			MPI_CHECK( MPI_File_write(mpifh, &gap, 1, MPI_UNSIGNED_LONG, &mpistat) );
			//mpioffset+=sizeof(unsigned long);
			//mpioffset=64;
			//MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,"ICRVQD",7,MPI_CHAR,&mpistat) );
			//MPI_CHECK( MPI_File_seek(mpifh,mpioffset,MPI_SEEK_SET) );
			// arg 2 type cast due to old MPI (<=V2) implementations (should be const void* now)
			MPI_CHECK( MPI_File_write(mpifh, (void*)const_cast<char*>("ICRVQD"), 6+1, MPI_CHAR, &mpistat) );
			mpioffset+=7;
			//
			//MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,"BB",3,MPI_CHAR,&mpistat) );
			//MPI_CHECK( MPI_File_seek(mpifh,mpioffset,MPI_SEEK_SET) );
			// arg 2 type cast due to old MPI (<=V2) implementations (should be const void* now)
			MPI_CHECK( MPI_File_write(mpifh, (void*)const_cast<char*>("BB"), 2+1, MPI_CHAR, &mpistat) );
			mpioffset+=3;
			numbb=(unsigned long)(num_procs);
			//MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&numbb,1,MPI_UNSIGNED_LONG,&mpistat) );
			//MPI_CHECK( MPI_File_seek(mpifh,mpioffset,MPI_SEEK_SET) );
			MPI_CHECK( MPI_File_write(mpifh, &numbb, 1, MPI_UNSIGNED_LONG, &mpistat) );
			mpioffset+=sizeof(unsigned long);
			//
			startidx=0;
		}
		MPI_CHECK( MPI_Exscan(&numParticles, &startidx, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD) );
		//
		mpioffset=64+7+3+sizeof(unsigned long)+ownrank*(6*sizeof(double)+2*sizeof(unsigned long));
		double bbmin[3],bbmax[3];
		bbmin[0]=domainDecomp->getBoundingBoxMin(0,domain);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&bbmin[0],1,MPI_DOUBLE,&mpistat) );
		mpioffset+=sizeof(double);
		bbmin[1]=domainDecomp->getBoundingBoxMin(1,domain);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&bbmin[1],1,MPI_DOUBLE,&mpistat) );
		mpioffset+=sizeof(double);
		bbmin[2]=domainDecomp->getBoundingBoxMin(2,domain);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&bbmin[2],1,MPI_DOUBLE,&mpistat) );
		mpioffset+=sizeof(double);
		bbmax[0]=domainDecomp->getBoundingBoxMax(0,domain);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&bbmax[0],1,MPI_DOUBLE,&mpistat) );
		mpioffset+=sizeof(double);
		bbmax[1]=domainDecomp->getBoundingBoxMax(1,domain);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&bbmax[1],1,MPI_DOUBLE,&mpistat) );
		mpioffset+=sizeof(double);
		bbmax[2]=domainDecomp->getBoundingBoxMax(2,domain);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&bbmax[2],1,MPI_DOUBLE,&mpistat) );
		mpioffset+=sizeof(double);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&startidx,1,MPI_UNSIGNED_LONG,&mpistat) );
		mpioffset+=sizeof(unsigned long);
		MPI_CHECK( MPI_File_write_at(mpifh,mpioffset,&numParticles,1,MPI_UNSIGNED_LONG,&mpistat) );
		mpioffset+=sizeof(unsigned long);
		Log::global_log->debug() << "[MPICheckpointWriter](" << ownrank << ")\tBB " << ":\t"
		                    << bbmin[0] << ", " << bbmin[1] << ", " << bbmin[2] << " - "
		                    << bbmax[0] << ", " << bbmax[1] << ", " << bbmax[2]
		                    << "\tstarting index=" << startidx << " numParticles=" << numParticles << std::endl;
		//
		MPI_Datatype mpidtParticleM, mpidtParticleD;
		ParticleData::getMPIType(mpidtParticleM);
		mpidtParticleD=mpidtParticleM;
		MPI_Aint mpidtParticleMsize=sizeof(ParticleData);	// !=MPI_Type_size
		int mpidtParticleDts;
		MPI_CHECK( MPI_Type_size(mpidtParticleM,&mpidtParticleDts) );
		// MPI_Type_get_extent
		/*
		ParticleData pd[2];
		MPI_Aint addr0,addr1;
		MPI_Get_address(pd,&addr0);
		MPI_Get_address(&pd[1],&addr1);
		mpidtParticleMsize=addr1-addr0;
		*/
		mpioffset=64+gap+startidx*mpidtParticleMsize;
		// arg 5 type cast due to old MPI (<=V2) implementations (should be const char* now)
		MPI_CHECK( MPI_File_set_view(mpifh, mpioffset, mpidtParticleM, mpidtParticleM, const_cast<char*>(mpidatarep), _mpiinfo) );
		Log::global_log->debug()
			<< "[MPICheckpointWriter]("
			<< ownrank
			<< ")\twriting molecule data for "
			<< numParticles
			<< " particles of size "
			<< mpidtParticleDts
			<< std::endl;
		//unsigned long writecounter=0;
		if(_particlesbuffersize>0)
		{
			ParticleData* particleStructBuffer=new ParticleData[_particlesbuffersize];
			unsigned long bufidx=0;
			for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
				// Log::global_log->debug()
				// << "MPICheckpointWriter["
				// << ownrank
				// << "]\t"
				// << pos->getID()
				// << "\t"
				// << pos->componentid()
				// << "\t"
				// << pos->r(0)
				// << ","
				// << pos->r(1)
				// << ","
				// << pos->r(2)
				// << std::endl;
				ParticleData::MoleculeToParticleData(particleStructBuffer[bufidx], *pos);
				++bufidx;
				if(bufidx==_particlesbuffersize)
				{
					//global_log->debug() << "MPICheckpointWriter[" << ownrank << "]\twriting" << _particlesbuffersize << " particles" << std::endl
					MPI_CHECK( MPI_File_write(mpifh, particleStructBuffer, _particlesbuffersize, mpidtParticleD, &mpistat) );
					//++writecounter;
					bufidx=0;
				}
			}
			if(bufidx>0)
			{
				//global_log->debug() << "MPICheckpointWriter[" << ownrank << "]\twriting" << bufidx << " particles" << std::endl
				MPI_CHECK( MPI_File_write(mpifh, particleStructBuffer, bufidx, mpidtParticleD, &mpistat) );
				//++writecounter;
			}
			delete[] particleStructBuffer;
		}
		else
		{
			ParticleData particleStruct;
			for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
				// Log::global_log->debug()
				// << "MPICheckpointWriter["
				// << ownrank
				// << "]\t"
				// << pos->getID()
				// << "\t"
				// << pos->componentid()
				// << "\t"
				// << pos->r(0)
				// << ","
				// << pos->r(1)
				// << ","
				// << pos->r(2)
				// << std::endl;
				ParticleData::MoleculeToParticleData(particleStruct, *pos);
				// Log::global_log->debug() << "MPICheckpointWriter[" << ownrank << "]\twriting particle" << std::endl
				MPI_CHECK( MPI_File_write(mpifh, &particleStruct, 1, mpidtParticleD, &mpistat) );
				//++writecounter;
				// saving a struct directly will also save padding zeros...
				//mpioffset+=mpidtParticleMsize;
			}
		}

		MPI_CHECK( MPI_File_close(&mpifh) );
		if(_measureTime)
		{
			MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
			double mpimeasuredtime=MPI_Wtime()-mpistarttime;
			// global_simulation->timers()->stop("MPI_CHECKPOINT_WRITER_INPUT");
			// double mpimeasuredtime=global_simulation->timers()->getTime("MPI_CHECKPOINT_WRITER_INPUT");
			if(ownrank==0) {
				Log::global_log->info() << "[MPICheckpointWriter]\tmeasured time: " << mpimeasuredtime << " sec (par., "
										<< num_procs << " proc.; " << numParticles_global << "*" << mpidtParticleDts
										<< "=" << numParticles_global * mpidtParticleDts << " Bytes)" << std::endl;
			}
		}
#else
		Log::global_log->info() << "[MPICheckpointWriter]\tnumber of particles: " << numParticles_global
		                   << "\t(*" << 2*sizeof(unsigned long)+13*sizeof(double) << "=" << numParticles_global*(2*sizeof(unsigned long)+13*sizeof(double)) << " Bytes in memory)"
				   << std::endl;
		unsigned long gap=7+3+sizeof(unsigned long)+(6*sizeof(double)+2*sizeof(unsigned long));
		unsigned int i;
		unsigned int offset=0;
		if (!_datarep.empty()) {
			Log::global_log->warning() << "[MPICheckpointWriter]\tsetting data representation (" << _datarep
									   << ") is not supported (yet) in sequential version" << std::endl;
		}
		// should use Timer instead
		struct timeval tod_start;
		if(_measureTime) {
			//global_log->debug() << "MPICheckpointWriter (" << filename << ")\tstart measuring time" << std::endl;
			gettimeofday( &tod_start, NULL );
			// global_simulation->timers()->start("MPI_CHECKPOINT_WRITER_INPUT");
		}
		//
		std::ofstream ostrm(filename.c_str(),std::ios::out|std::ios::binary);
		ostrm << _magicVersion;
		offset+=strlen(_magicVersion);
		for(i=0;i<64-offset-sizeof(unsigned long)-sizeof(int);++i) ostrm << '\0';
		ostrm.write((const char*)&_endiannesstest,sizeof(int));
		ostrm.write((char*)&gap,sizeof(unsigned long));
		//offset=64
		//ostrm.seekp(offset);
		ostrm << "ICRVQD" << '\0';
		//offset+=7;
		ostrm << "BB" << '\0';
		ostrm.write((char*)&numbb,sizeof(unsigned long));
		//offset+=3+sizeof(unsigned long);
		double bbmin[3],bbmax[3];
		bbmin[0]=domainDecomp->getBoundingBoxMin(0,domain);
		bbmin[1]=domainDecomp->getBoundingBoxMin(1,domain);
		bbmin[2]=domainDecomp->getBoundingBoxMin(2,domain);
		bbmax[0]=domainDecomp->getBoundingBoxMax(0,domain);
		bbmax[1]=domainDecomp->getBoundingBoxMax(1,domain);
		bbmax[2]=domainDecomp->getBoundingBoxMax(2,domain);
		ostrm.write((char*)bbmin,3*sizeof(double));
		ostrm.write((char*)bbmax,3*sizeof(double));
		//offset+=6*sizeof(double);
		unsigned long startidx=0;
		ostrm.write((char*)&startidx,sizeof(unsigned long));
		ostrm.write((char*)&numParticles,sizeof(unsigned long));
		//offset+=2*sizeof(unsigned long);
		for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
			unsigned long id=pos->getID();
			ostrm.write((char*)&id,sizeof(unsigned long));
			//unsigned int componentid=pos->componentid(); // to be compatible to struct padding
			unsigned long componentid=pos->componentid();
			ostrm.write((char*)&componentid,sizeof(unsigned long));
			double r[3],v[3],q[4],D[3];
			for(i=0;i<3;++i) r[i]=pos->r(i);
			ostrm.write((char*)r,3*sizeof(double));
			for(i=0;i<3;++i) v[i]=pos->v(i);
			ostrm.write((char*)v,3*sizeof(double));
			q[0]=pos->q().qw();
			q[1]=pos->q().qx();
			q[2]=pos->q().qy();
			q[3]=pos->q().qz();
			ostrm.write((char*)q,4*sizeof(double));
			for(i=0;i<3;++i) D[i]=pos->D(i);
			ostrm.write((char*)D,3*sizeof(double));
			//offset+=2*sizeof(unsigned long)+13*sizeof(double);
		}
		ostrm.close();
		if(_measureTime)
		{
			struct timeval tod_end;
			gettimeofday( &tod_end, NULL );
			double measuredtime=(double)(tod_end.tv_sec-tod_start.tv_sec)+(double)(tod_end.tv_usec-tod_start.tv_usec)/1.E6;
			// global_simulation->timers()->stop("MPI_CHECKPOINT_WRITER_INPUT");
			// double measuredtime=global_simulation->timers()->getTime("MPI_CHECKPOINT_WRITER_INPUT");
			Log::global_log->info() << "[MPICheckpointWriter]\tmeasured time: " << measuredtime << " sec (seq.)" << std::endl;
		}
#endif
	}
}

void MPICheckpointWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
								 Domain * /*domain*/) {}

