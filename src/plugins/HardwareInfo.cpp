/*
 * HardwareInfo.cpp
 *
 *  Created on: 16 Dec 2025
 *      Author: amartyads
 */

#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#ifdef __GLIBC__
#include <sched.h>	// sched_getcpu(), getcpu(int*, int*)
#endif
#ifndef _WIN32		 // only OS not POSIX compliant
#include <unistd.h>	 // gethostname(char*, int)
#endif

#include <fstream>
#include <iostream>
#include <sstream>

#include "HardwareInfo.h"
#include "WrapOpenMP.h"
#include "parallel/DomainDecompBase.h"	// getCommunicator()
#include "utils/Logger.h"

void HardwareInfo::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("filename", _filename);
	if (_filename != "") {
		_filename += ".json";
		Log::global_log->info() << "[" << getPluginName() << "] Filename: " << _filename << std::endl;
	}
}

void HardwareInfo::init(ParticleContainer*, DomainDecompBase* domainDecomp, Domain*) {
	// defaults
	_rank = 0;
	_totalRanks = 1;
	_processorName = "";
	populateData(domainDecomp);
	if (_filename != "")
		writeDataToFile(domainDecomp);
	else
		printDataToStdout();
}

void HardwareInfo::populateData(const DomainDecompBase* domainDecomp) {
	char cStyleProcName[1024] = {"default"};
#ifndef _WIN32
	gethostname(cStyleProcName, 1023);	// from unistd.h
#endif
	_threadData.resize(mardyn_get_max_threads());
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
#pragma omp parallel shared(_threadData)
#endif
	{
		int thread = mardyn_get_thread_num();
		int totalThreads = mardyn_get_num_threads();
		unsigned int openMPCPUID, openMPNUMA;
#if __GLIBC__ >= 2 && __GLIBC_MINOR__ >= 29
		getcpu(&openMPCPUID, &openMPNUMA);	// from sched.h
#elif defined(__GLIBC__)
		openMPCPUID = sched_getcpu();  // from sched.h
		openMPNUMA = -1;
#else
		openMPCPUID = openMPNUMA = -1;
#endif
		_threadData[thread] = ThreadwiseInfo(thread, totalThreads, openMPCPUID, openMPNUMA);
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	_dataPopulated = true;
}

const void HardwareInfo::printDataToStdout() {
	if (!_dataPopulated) {	// sanity check
		std::ostringstream msg;
		msg << "[" << getPluginName() << "] data not populated!" << std::endl;
		MARDYN_EXIT(msg.str());
	}
	std::ostringstream ss;
	for (auto data : _threadData) {
		ss << "[" << getPluginName() << "] Thread " << data.thread << " out of " << data.totalThreads << " threads";
		if (data.numa != -1)
			ss << ", NUMA domain " << data.numa;
		ss << ", rank " << _rank << " out of " << _totalRanks << " ranks, running on " << _processorName;
		if (data.cpuID != -1)
			ss << " with CPU index " << data.cpuID << std::endl;
		Log::global_log->set_mpi_output_all();
		Log::global_log->info() << ss.str();
		Log::global_log->set_mpi_output_root(0);
		ss.str(std::string());	// clear ostringstream
	}
}

const void HardwareInfo::writeDataToFile(const DomainDecompBase* domainDecomp) {
	if (!_dataPopulated) {	// sanity check
		std::ostringstream msg;
		msg << "[" << getPluginName() << "] data not populated!" << std::endl;
		MARDYN_EXIT(msg.str());
	}
	Log::global_log->info() << "[" << getPluginName() << "] Writing to file: " << _filename << std::endl;
	// put all local data into a string ready to write to file
	std::ostringstream outputStringSS;
	// add header data if rank 0
	if (_rank == 0) {
		outputStringSS << "{\n";
		outputStringSS << "\t\"total_ranks\": " << _totalRanks << ",\n";
		outputStringSS << "\t\"rank_data\": {";
		// write header
		std::ofstream outputFile(_filename);
		outputFile << outputStringSS.str();
		outputFile.close();
		outputStringSS.str(std::string());	// clear ostringstream
	}
	// add all thread data
	// since trailing commas are not allowed, put data from rank 0 at end without comma
	std::string outputString = convertFullDataToJson();
#ifdef ENABLE_MPI
	auto curComm = domainDecomp->getCommunicator();
	// taken from DomainDecompMPIBase::printDecomp
	MPI_File parFile;
	MPI_File_open(curComm, _filename.c_str(),
				  MPI_MODE_WRONLY | MPI_MODE_APPEND | MPI_MODE_CREATE, MPI_INFO_NULL, &parFile);
	unsigned long writeSize = outputString.size();
	unsigned long offset = 0;
	if (_rank == 0) {
		MPI_Offset fileEndOffset;
		MPI_File_seek(parFile, 0, MPI_SEEK_END);
		MPI_File_get_position(parFile, &fileEndOffset);
		writeSize += fileEndOffset;
		MPI_Exscan(&writeSize, &offset, 1, MPI_UINT64_T, MPI_SUM, curComm);
		offset += fileEndOffset;
	} else {
		MPI_Exscan(&writeSize, &offset, 1, MPI_UINT64_T, MPI_SUM, curComm);
	}
	MPI_File_write_at(parFile, static_cast<MPI_Offset>(offset), outputString.c_str(),
					  static_cast<int>(outputString.size()), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&parFile);
#endif
	// close remaining braces
	if (_rank == 0) {
		std::ofstream outputFile(_filename, std::ios_base::app);
#ifndef ENABLE_MPI	// write serial data if MPI not enabled
		outputFile << outputString;
#endif
		outputFile << "\n\t}\n}";  // ending brace if rank 0
		outputFile.close();
	}
}

const std::string HardwareInfo::convertFullDataToJson() {
	std::ostringstream rankInfo;
	rankInfo << "\n\t\t\"" << _rank << "\": {\n";
	rankInfo << "\t\t\t\"node_name\": \"" << _processorName << "\",\n";
	rankInfo << "\t\t\t\"total_threads\": \"" << _threadData[0].totalThreads << "\",\n";
	rankInfo << "\t\t\t\"thread_data\": {\n";

	for (auto it = _threadData.begin(); it != _threadData.end(); ++it) {
		if (it != _threadData.begin())
			rankInfo << ",\n";
		rankInfo << "\t\t\t\"" << it->thread << "\": {";
		if (it->cpuID != -1)
			rankInfo << "\"cpu_ID\": " << it->cpuID;
		if (it->numa != -1)
			rankInfo << ", \"numa_domain\": " << it->numa;
		rankInfo << "}";
	}

	rankInfo << "\n\t\t\t}";	   // close threads
	rankInfo << "\n\t\t}";		   // close rank
	if (_rank != _totalRanks - 1)  // no trailing comma for final rank
		rankInfo << ",";
	return rankInfo.str();
}
