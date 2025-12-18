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

// #include "utils/compile_info.h"

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
		writeDataToFile();
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
		getcpu(&openMPCPUID, &openMPNUMA);	// from sched.h
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
		ss << "[" << getPluginName() << "] Thread " << data.thread << " out of " << data.totalThreads
		   << " threads, NUMA domain " << data.numa << ", rank " << _rank << " out of " << _totalRanks
		   << " ranks, running on " << _processorName << " with CPU index " << data.cpuID << std::endl;
		Log::global_log->set_mpi_output_all();
		Log::global_log->info() << ss.str();
		Log::global_log->set_mpi_output_root(0);
		ss.str(std::string());
	}
}

void HardwareInfo::writeDataToFile() {
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
		outputStringSS << "\t\"ranks\": {";
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
	if (_rank != 0) {
		MPI_File parFile;
		// MPI_File_open
	}
#endif
	// write thread data from rank 0 to file to ensure no comma at end
	if (_rank == 0) {
		std::ofstream outputFile(_filename, std::ios_base::app);
		outputFile << outputString;
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
		rankInfo << "\"cpu_ID\": " << it->cpuID << ", ";
		rankInfo << "\"numa_domain\": " << it->numa;
		rankInfo << "}";
	}

	rankInfo << "\n\t\t\t}";  // close threads
	rankInfo << "\n\t\t}";	  // close rank
	if (_rank != 0)
		rankInfo << ",";
	return rankInfo.str();
}
