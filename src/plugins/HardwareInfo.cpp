/*
 * HardwareInfo.cpp
 *
 *  Created on: 16 Dec 2025
 *      Author: amartyads
 */

#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include <sched.h>
#include <unistd.h>

#include <sstream>

#include "HardwareInfo.h"
#include "WrapOpenMP.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"
//#include "utils/compile_info.h"

void HardwareInfo::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	// defaults
	_thread = _rank = 0;
	_totalRanks = _totalThreads = 1;
	gethostname(_processorName, 1023);  // from unistd.h
	// mpi
#ifdef ENABLE_MPI
	auto curCommunicator = domainDecomp->getCommunicator();
	MPI_CHECK(MPI_Comm_size(curCommunicator, &_totalRanks));
	MPI_CHECK(MPI_Comm_rank(curCommunicator, &_rank));
	int name_len;
	MPI_CHECK(MPI_Get_processor_name(_processorName, &name_len));
#endif

	// openmp
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		_thread = mardyn_get_thread_num();
		_totalThreads = mardyn_get_num_threads();
		getcpu(&_openMPCPUID, &_openMPNUMA);  // from sched.h
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	// output
	std::stringstream ss;
	ss << "[HardwareInfo] Thread " << _thread << " out of " << _totalThreads << " threads, NUMA domain " << _openMPNUMA
	   << ", rank " << _rank << " out of " << _totalRanks << " ranks, running on " << _processorName
	   << " with CPU index " << _openMPCPUID << std::endl;
	Log::global_log->set_mpi_output_all();
	Log::global_log->info() << ss.str();
	Log::global_log->set_mpi_output_root(0);
}
