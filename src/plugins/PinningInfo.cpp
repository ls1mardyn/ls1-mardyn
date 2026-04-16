/*
 * PinningInfo.cpp
 *
 *  Created on: 16 Dec 2025
 *      Author: amartyads
 */

#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __GLIBC__
#include <sched.h>	// sched_getcpu(), getcpu(int*, int*)
#endif
#ifdef ENABLE_JSON
#include <nlohmann/json.hpp>
using json = nlohmann::json;
#endif

#include <sys/sysinfo.h>  // sysinfo
#include <sys/utsname.h>  // uname

#include <fstream>
#include <iomanip>	// std::fixed, std::setprecision
#include <sstream>

#include "PinningInfo.h"
#include "parallel/DomainDecompBase.h"	// getCommunicator()
#include "utils/Logger.h"
#include "utils/String_utils.h"	 // trim()

void PinningInfo::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("filename", _filename);
	if (_filename != "") {
#ifndef ENABLE_JSON
		std::ostringstream msg;
		msg << "[" << getPluginName() << "] Cannot write to JSON file since ENABLE_JSON not specified in cmake!"
			<< std::endl;
		MARDYN_EXIT(msg.str());
#else
		_filename += ".json";
		Log::global_log->info() << "[" << getPluginName() << "] Filename: " << _filename << std::endl;
#endif
	}
}

void PinningInfo::init(ParticleContainer*, DomainDecompBase* domainDecomp, Domain*) {
	// defaults
	_rank = 0;
	_totalRanks = 1;
	_nodeName = "";
#ifndef __GLIBC__
	Log::global_log->warning()
		<< "[" << getPluginName()
		<< "] sched.h cannot be loaded since glibc was not used! Thread pinning data replaced with -1!" << std::endl;
#elif (__GLIBC__ < 2) || (__GLIBC__ == 2 && __GLIBC_MINOR__ < 29)
	Log::global_log->warning()
		<< "[" << getPluginName()
		<< "] glibc version too low to show NUMA info (>=2.29 required)! NUMA data replaced with -1!" << std::endl;
#endif
	populateData(domainDecomp);
	if (_filename != "")
		writeDataToFile(domainDecomp);
	else
		printDataToStdout();
}

void PinningInfo::populateData(DomainDecompBase* domainDecomp) {
	char cStyleNodeName[1024];
	struct utsname utsnameData;	 // from sys/utsname.h
	if (uname(&utsnameData) < 0) {
		std::ostringstream msg;
		msg << "[" << getPluginName() << "] Error using uname() from sys/utsname.h!" << std::endl;
		MARDYN_EXIT(msg.str());
	}
	strcpy(cStyleNodeName, utsnameData.nodename);

	// rank level data
#ifdef ENABLE_MPI
	auto curComm = domainDecomp->getCommunicator();
	MPI_CHECK(MPI_Comm_size(curComm, &_totalRanks));
	MPI_CHECK(MPI_Comm_rank(curComm, &_rank));
	int name_len;  // value not used, placeholder was needed for get_processor_name
	MPI_CHECK(MPI_Get_processor_name(cStyleNodeName, &name_len));
#endif
	_nodeName = std::string(cStyleNodeName);

	// thread level data
	_maxThreads = 1;
#ifdef _OPENMP
	_maxThreads = omp_get_max_threads();
#endif
	_threadData.resize(_maxThreads);

#ifdef _OPENMP
#pragma omp parallel shared(_threadData)
#endif
	{
		int thread = 0;
#ifdef _OPENMP
		thread = omp_get_thread_num();
#endif
		int openMPCPUID, openMPNUMA;
#if __GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 29)
		unsigned int tempCPU, tempNUMA;
#ifdef _OPENMP
#pragma omp critical(getcpu)  // cannot find documentation explicitly saying this is thread-safe, hence precautionary
#endif
		getcpu(&tempCPU, &tempNUMA);  // from sched.h
		// following casts are only an issue if values > INT_MAX, which is unusual anyway
		openMPCPUID = (int)tempCPU;
		openMPNUMA = (int)tempNUMA;
#elif defined(__GLIBC__)
		openMPCPUID = sched_getcpu();  // from sched.h, documentation says this is thread-safe
		openMPNUMA = -1;
#else
		openMPCPUID = openMPNUMA = -1;
#endif
		_threadData[thread] = ThreadwiseInfo(thread, openMPCPUID, openMPNUMA);
	}

	// RAM information
	_maxRam = "N/A";
	struct sysinfo sysinfoData;	 // from sys/sysinfo.h
	if (sysinfo(&sysinfoData) < 0) {
		std::ostringstream msg;
		msg << "[" << getPluginName() << "] Error using sysinfo() from sys/sysinfo.h!" << std::endl;
		MARDYN_EXIT(msg.str());
	}
	float ramGB = sysinfoData.totalram / 1000.0 / 1000.0 / 1000.0;
	float ramGiB = sysinfoData.totalram / 1024.0 / 1024.0 / 1024.0;
	std::ostringstream ramDataCollect;
	ramDataCollect << std::fixed << std::setprecision(2) << ramGB << " GB (" << ramGiB << " GiB)";
	_maxRam = ramDataCollect.str();

	// CPU information
	_cpuArch = std::string(utsnameData.machine);
	std::ifstream procCPUStream("/proc/cpuinfo");
	std::string fileLine, key, value = "N/A";
	char sep = ':';	 // separator for key-value pairs in /proc/cpuinfo
	if (!procCPUStream) {
		Log::global_log->warning() << "[" << getPluginName() << "] Could not open /proc/cpuinfo!" << std::endl;
	} else {
		if (_cpuArch == "aarch64") {  // ARM, proc/cpuinfo more complex
			std::ostringstream valueBuilder;
			valueBuilder << "ARM ";
			short valsFound = 0;
			while (getline(procCPUStream, fileLine)) {
				key = string_utils::trim(fileLine.substr(0, fileLine.find(sep)));
				if (key.substr(0, 3) == "CPU") {
					key = string_utils::trim(key.substr(4));  // trim out "CPU" from the key
					value = string_utils::trim(fileLine.substr(fileLine.find(sep) + 1));
					valueBuilder << key << ": " << value << " ";
					valsFound++;
				}
				if (valsFound > 4)	// grabs CPU implementer, architecture, variant, part, and revision, and then exits
					break;
			}
			value = string_utils::trim(valueBuilder.str());
		} else if (_cpuArch == "x86_64") {	// Intel/AMD, proc/cpuinfo simpler
			while (getline(procCPUStream, fileLine)) {
				key = string_utils::trim(fileLine.substr(0, fileLine.find(sep)));
				if (key == "model name") {	// only take model name
					value = string_utils::trim(fileLine.substr(fileLine.find(sep) + 1));
					break;
				}
			}
		}  // no other architectures tested, hence no default case
	}
	procCPUStream.close();
	_cpuInfo = value;
	_dataPopulated = true;
}

void PinningInfo::printDataToStdout() {
	if (!_dataPopulated) {	// sanity check
		std::ostringstream msg;
		msg << "[" << getPluginName() << "] Data not populated!" << std::endl;
		MARDYN_EXIT(msg.str());
	}
	std::ostringstream ss;
	for (auto data : _threadData) {
		ss << "[" << getPluginName() << "] Thread " << data.thread << " out of " << _maxThreads << " threads";
		ss << ", NUMA domain " << data.numa;
		ss << ", rank " << _rank << " out of " << _totalRanks << " ranks, running on " << _nodeName;
		ss << " with CPU index " << data.cpuID;
		ss << ", CPU info: " << _cpuInfo << ", architecture: " << _cpuArch;
		ss << ", Max RAM available: " << _maxRam << std::endl;
		Log::global_log->set_mpi_output_all();
		Log::global_log->info() << ss.str();
		Log::global_log->set_mpi_output_root(0);
		ss.str(std::string());	// clear ostringstream
	}
}

void PinningInfo::writeDataToFile(DomainDecompBase* domainDecomp) {
#ifndef ENABLE_JSON	 // should never be reached, but failsafe
	std::ostringstream msg;
	msg << "[" << getPluginName() << "] Cannot write to JSON file since ENABLE_JSON not specified in cmake! Exiting..."
		<< std::endl;
	MARDYN_EXIT(msg.str());
#else
	if (!_dataPopulated) {	// sanity check
		std::ostringstream msg;
		msg << "[" << getPluginName() << "] Data not populated!" << std::endl;
		MARDYN_EXIT(msg.str());
	}
	Log::global_log->info() << "[" << getPluginName() << "] Writing to file: " << _filename << std::endl;

	// create a json object on all ranks to hold the relevant data
	json threadDataJson;

	// write data to json object
	threadDataJson["node_name"] = _nodeName;
	threadDataJson["cpu_info"] = _cpuInfo;
	threadDataJson["cpu_arch"] = _cpuArch;
	threadDataJson["max_avail_ram"] = _maxRam;
	threadDataJson["total_threads"] = _maxThreads;
	for (auto it = _threadData.begin(); it != _threadData.end(); ++it) {
		threadDataJson["thread_data"][std::to_string(it->thread)]["cpu_id"] = it->cpuID;
		threadDataJson["thread_data"][std::to_string(it->thread)]["numa_domain"] = it->numa;
	}

	// collect data and write to file
	if (_rank == 0) {
		// create new json object to collect data, and hold common info
		json collDataJson;
		collDataJson["total_ranks"] = _totalRanks;
		collDataJson["rank_data"]["0"] = threadDataJson;
		// collect data from other ranks
#ifdef ENABLE_MPI
		MPI_Status status;
		std::vector<std::uint8_t> dataToReceive;
		for (int i = 1; i < _totalRanks; i++) {
			MPI_Probe(i, 0, domainDecomp->getCommunicator(), &status);
			int count;
			MPI_Get_count(&status, MPI_UINT8_T, &count);
			dataToReceive.resize(count);
			MPI_Recv(dataToReceive.data(), count, MPI_UINT8_T, i, 0, domainDecomp->getCommunicator(), &status);
			collDataJson["rank_data"][std::to_string(i)] = json::from_bson(dataToReceive);
		}
#endif
		// write to file
		std::ofstream outputFile(_filename);
		outputFile << collDataJson.dump(4);	 // 4 is the tab size
		outputFile.close();
#ifdef ENABLE_MPI
	} else {  // _rank != 0
		// send data to rank 0
		std::vector<std::uint8_t> dataToSend = json::to_bson(threadDataJson);
		MPI_Send(dataToSend.data(), dataToSend.size(), MPI_UINT8_T, 0, 0, domainDecomp->getCommunicator());
#endif
	}
#endif
}
