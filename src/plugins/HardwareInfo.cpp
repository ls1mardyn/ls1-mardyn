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

#include <iostream>
#include <sstream>

#include "HardwareInfo.h"
#include "WrapOpenMP.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"

// #include "utils/compile_info.h"

void HardwareInfo::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	// defaults
	_rank = 0;
	_totalRanks = 1;
	_processorName = "";
	populateData(domainDecomp);
	printDataToStdout();
}

void HardwareInfo::populateData(DomainDecompBase* domainDecomp) {
	char cStyleProcName[1024];
	gethostname(cStyleProcName, 1023);	// from unistd.h
	threadData.resize(mardyn_get_max_threads());
	threadData[0] = ThreadwiseInfo();
	// mpi
#ifdef ENABLE_MPI
	auto curCommunicator = domainDecomp->getCommunicator();
	MPI_CHECK(MPI_Comm_size(curCommunicator, &_totalRanks));
	MPI_CHECK(MPI_Comm_rank(curCommunicator, &_rank));
	int name_len;
	MPI_CHECK(MPI_Get_processor_name(cStyleProcName, &name_len));
#endif
	_processorName = std::string(cStyleProcName);
	// openmp
#ifdef _OPENMP
#pragma omp parallel shared(threadData)
#endif
	{
		int thread = mardyn_get_thread_num();
		int totalThreads = mardyn_get_num_threads();
		unsigned int openMPCPUID, openMPNUMA;
		getcpu(&openMPCPUID, &openMPNUMA);	// from sched.h
		threadData[thread] =
			ThreadwiseInfo(_rank, _totalRanks, thread, totalThreads, openMPCPUID, openMPNUMA, _processorName);
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	_dataPopulated = true;
}

void HardwareInfo::printDataToStdout() {
	std::stringstream ss;
	for (auto data : threadData) {
		ss << "[HardwareInfo] Thread " << data.thread << " out of " << data.totalThreads << " threads, NUMA domain "
		   << data.numa << ", rank " << data.rank << " out of " << data.totalRanks << " ranks, running on "
		   << data.processorName << " with CPU index " << data.cpuID << std::endl;
		Log::global_log->set_mpi_output_all();
		Log::global_log->info() << ss.str();
		Log::global_log->set_mpi_output_root(0);
		ss.str(std::string());
	}
}
