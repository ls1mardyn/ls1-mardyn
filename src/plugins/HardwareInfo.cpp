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

#include <fstream>
#include <iostream>
#include <sstream>

#include "HardwareInfo.h"
#include "WrapOpenMP.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"

void HardwareInfo::readXML(XMLfileUnits& xmlconfig) {
	_filename = "";
	xmlconfig.getNodeValue("filename", _filename);
	if (_filename != "") {
		_filename += ".json";
		Log::global_log->info() << "[" << getPluginName() << "] Filename: " << _filename << std::endl;
	}
}

void HardwareInfo::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	// defaults
	_rank = 0;
	_totalRanks = 1;
	_processorName = "";
	populateData(domainDecomp);
	printDataToStdout();
	if (_filename != "")
		writeDataToFile(domainDecomp);
}

void HardwareInfo::populateData(DomainDecompBase* domainDecomp) {
	char cStyleProcName[1024];
	gethostname(cStyleProcName, 1023);	// from unistd.h
	threadData.resize(mardyn_get_max_threads());
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
#if __GLIBC__ > 1 && __GLIBC_MINOR__ >= 29
		getcpu(&openMPCPUID, &openMPNUMA);	// from sched.h
#elif defined(__GLIBC__)
		openMPCPUID = sched_getcpu();  // from sched.h
		openMPNUMA = -1;
#else
		openMPCPUID = openMPNUMA = -1;
#endif
		threadData[thread] = ThreadwiseInfo(thread, totalThreads, openMPCPUID, openMPNUMA);
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	_dataPopulated = true;
}

void HardwareInfo::printDataToStdout() {
	if (!_dataPopulated)
		return;
	std::ostringstream ss;
	for (auto data : threadData) {
		ss << "[" << getPluginName() << "] Thread " << data.thread << " out of " << data.totalThreads << " threads";
		if (data.numa != -1)
			ss << ", NUMA domain " << data.numa;
		ss << ", rank " << _rank << " out of " << _totalRanks << " ranks, running on " << _processorName;
		if (data.cpuID != -1)
			ss << " with CPU index " << data.cpuID << std::endl;
		Log::global_log->set_mpi_output_all();
		Log::global_log->info() << ss.str();
		Log::global_log->set_mpi_output_root(0);
		ss.str(std::string());
	}
}

void HardwareInfo::writeDataToFile(DomainDecompBase* domainDecomp) {
	if (!_dataPopulated)
		return;
	Log::global_log->info() << "[" << getPluginName() << "] Writing to file: " << _filename << std::endl;
	// put all local data into a string ready to write to file
	std::ostringstream outputStringSS;
	// add header data if rank 0
	if (_rank == 0) {
		outputStringSS << "{\n";
		outputStringSS << "\t\"total_ranks\": " << _totalRanks << ",\n";
		outputStringSS << "\t\"total_threads\": " << threadData[0].totalThreads << ",\n";
		outputStringSS << "\t\"rank_data\": {";
		// write header
		std::ofstream outputFile(_filename);
		outputFile << outputStringSS.str();
		outputFile.close();
		outputStringSS.str(std::string());
	}
	// add all thread data
	// since trailing commas are not allowed, put data from rank 0 at end without comma
	std::string outputString = convertFullDataToJson();
#ifdef ENABLE_MPI
	// taken from DomainDecompMPIBase::printDecomp
	MPI_File parFile;
	MPI_File_open(domainDecomp->getCommunicator(), _filename.c_str(),
				  MPI_MODE_WRONLY | MPI_MODE_APPEND | MPI_MODE_CREATE, MPI_INFO_NULL, &parFile);
	unsigned long writeSize = outputString.size();
	unsigned long offset = 0;
	if (_rank == 0) {
		MPI_Offset fileEndOffset;
		MPI_File_seek(parFile, 0, MPI_SEEK_END);
		MPI_File_get_position(parFile, &fileEndOffset);
		writeSize += fileEndOffset;
		MPI_Exscan(&writeSize, &offset, 1, MPI_UINT64_T, MPI_SUM, domainDecomp->getCommunicator());
		offset += fileEndOffset;
	} else {
		MPI_Exscan(&writeSize, &offset, 1, MPI_UINT64_T, MPI_SUM, domainDecomp->getCommunicator());
	}
	MPI_File_write_at(parFile, static_cast<MPI_Offset>(offset), outputString.c_str(),
					  static_cast<int>(outputString.size()), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&parFile);
#endif
	// close remaining braces
	if (_rank == 0) {
		std::ofstream outputFile(_filename, std::ios_base::app);
		outputFile << "\n\t}\n}";  // ending brace if rank 0
		outputFile.close();
	}
}

std::string HardwareInfo::convertFullDataToJson() {
	std::ostringstream rankInfo;
	rankInfo << "\n\t\t\"" << _rank << "\": {\n";
	rankInfo << "\t\t\t\"node_name\": \"" << _processorName << "\",\n";
	rankInfo << "\t\t\t\"thread_data\": {\n";

	for (auto it = threadData.begin(); it != threadData.end(); ++it) {
		if (it != threadData.begin())
			rankInfo << ",\n";
		rankInfo << "\t\t\t\"" << it->thread << "\": {";
		if (it->cpuID != -1)
			rankInfo << "\"cpu_ID\": " << it->cpuID;
		if (it->numa != -1)
			rankInfo << ", \"numa_domain\": " << it->numa;
		rankInfo << "}";
	}

	rankInfo << "\n\t\t\t}";  // close threads
	rankInfo << "\n\t\t}";	  // close rank
	if (_rank != _totalRanks - 1)
		rankInfo << ",";
	return rankInfo.str();
}
