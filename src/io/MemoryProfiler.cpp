/*
 * MemoryProfiler.cpp
 *
 *  Created on: May 9, 2017
 *      Author: seckler
 */

#include "MemoryProfiler.h"

#include <cstdlib>
#include <fstream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include "sys/sysinfo.h"
#ifndef _SX
#include <sys/sysinfo.h>
#endif
#include "utils/Logger.h"

MemoryProfiler::MemoryProfiler() :
		_list(), _hugePageSize(0) {

	std::ifstream meminfo_file("/proc/meminfo");
	std::string token;
	while(meminfo_file >> token) {
		if (token == "Hugepagesize:") {
			meminfo_file >> _hugePageSize;
			break;
		}
	}
	if (const char *env_p = std::getenv("HUGETLB_DEFAULT_PAGE_SIZE")) {
		std::string hugepsString(env_p);
		_hugePageSize = std::stoi(hugepsString);
		if( hugepsString.find("G") != std::string::npos or hugepsString.find("g") != std::string::npos ){
			_hugePageSize *= 1024 * 1024;
		} else if( hugepsString.find("M") != std::string::npos or hugepsString.find("m") != std::string::npos ){
			_hugePageSize *= 1024;
		}
	}
}

void MemoryProfiler::registerObject(MemoryProfilable** object) {
	_list.push_back(object);
	Log::global_log->debug() << "MemoryProfiler: added object" << std::endl;
}

void MemoryProfiler::doOutput(const std::string& myString) {
	printGeneralInfo(myString);

	// further info
	Log::global_log->debug() << "MemoryProfiler: number of objects: " << _list.size() << std::endl;
	for (auto item : _list) {
		if ((*item) == nullptr) {
			continue;
		}
		Log::global_log->info() << "\t\t" << (*item)->getName() << ": " << (*item)->getTotalSize() / 1.e6 << " MB"
				<< std::endl;
		(*item)->printSubInfo(3);
	}
}

unsigned long long MemoryProfiler::getCachedSize() {
	unsigned long long cached_size = 0;
	std::ifstream meminfo_file("/proc/meminfo");
	std::string token;
	while(meminfo_file >> token) {
		if (token == "Cached:") {
			meminfo_file >> cached_size;
			break;
		}
	}
	return cached_size;
}

int MemoryProfiler::countHugePages() {
	std::ifstream numa_maps_file("/proc/self/numa_maps");

	int hugepagecount = 0;
	std::string line;
	while(std::getline(numa_maps_file, line)) {
		//if (/huge.*dirty=(\d+)/) {
		if (line.find("huge") != std::string::npos) {
			const std::regex re("dirty=(\\d+) ");
			std::smatch match;
			if (std::regex_match(line, match, re)) {
				hugepagecount += std::stoi(match[1].str());
			}
		}
	}
	return hugepagecount;
}

int MemoryProfiler::getOwnMemory() { //Note: this value is in KB!
	int vmrss = 0;
	std::ifstream self_status_file("/proc/self/status");
	std::string token;
	while(self_status_file >> token) {
		if (token == "VmRSS:") {
			self_status_file >> vmrss;
			break;
		}
	}
	return vmrss;
}

void MemoryProfiler::printGeneralInfo(const std::string& myString) {
#ifndef _SX
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	long long totalMem = memInfo.totalram * memInfo.mem_unit / 1024 / 1024;
	long long usedMem = ((memInfo.totalram - memInfo.freeram - memInfo.bufferram) * memInfo.mem_unit / 1024
			- getCachedSize()) / 1024;
	std::stringstream additionalinfo;
	if (myString.length() > 0) {
		additionalinfo << " (" << myString << ")";
	}
	additionalinfo << ":" << std::endl;
	Log::global_log->info() << "Memory consumption" << additionalinfo.str();
	Log::global_log->info() << "\tMemory usage (System total):\t" << usedMem << " MB out of " << totalMem << " MB ("
			<< usedMem * 100. / totalMem << "%)" << std::endl;
	double ownMem = getOwnMemory() / 1024.;
	double hugeMem = countHugePages() * _hugePageSize / 1024.;
	Log::global_log->info() << "\tBy own process:\t\t\t" << ownMem+hugeMem << " MB (" << (ownMem+hugeMem) * 100. / totalMem
			<< "% of total memory)" << std::endl;
	if (hugeMem != 0) {
		Log::global_log->info() << "\t\t\tnormal:\t\t" << ownMem << " MB (" << (ownMem) * 100. / totalMem
				<< "% of total memory)" << std::endl;

		Log::global_log->info() << "\t\t\thugePages:\t" << hugeMem << " MB (" << (hugeMem) * 100. / totalMem
				<< "% of total memory)" << "\tHPS(kB):\t" << _hugePageSize << std::endl;
	}
#else
    Log::global_log->warning() << "MemoryProfiler of ls1 is not available for this platform." << std::endl;
#endif
}
