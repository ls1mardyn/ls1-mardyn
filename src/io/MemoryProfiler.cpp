/*
 * MemoryProfiler.cpp
 *
 *  Created on: May 9, 2017
 *      Author: seckler
 */

#include "MemoryProfiler.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include "sys/types.h"
#include "sys/sysinfo.h"

#include "utils/Logger.h"

MemoryProfiler::MemoryProfiler() :
		_list(), _hugePageSize(0) {

	const size_t MAXLEN = 1024;
	FILE *fp;
	char buf[MAXLEN];
	fp = fopen("/proc/meminfo", "r");
	char *p1;
	while (fgets(buf, MAXLEN, fp)) {
		p1 = strstr(buf, "Hugepagesize:");
		if (p1 != NULL) {
			int colon = ':';
			p1 = strchr(buf, colon) + 1;
			_hugePageSize = atoi(p1);
		}
	}
	if (std::getenv("HUGETLB_DEFAULT_PAGE_SIZE")!=NULL){
		std::string hugepsString = std::getenv("HUGETLB_DEFAULT_PAGE_SIZE");
		_hugePageSize = std::stoi(hugepsString);
		if( hugepsString.find("G") != std::string::npos or hugepsString.find("g") != std::string::npos ){
			_hugePageSize *= 1024 * 1024;
		} else if( hugepsString.find("M") != std::string::npos or hugepsString.find("m") != std::string::npos ){
			_hugePageSize *= 1024;
		}
	}
	fclose(fp);
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
	const size_t MAXLEN = 1024;
	FILE *fp;
	char buf[MAXLEN];
	fp = fopen("/proc/meminfo", "r");
	while (fgets(buf, MAXLEN, fp)) {
		char *p1 = strstr(buf, "Cached:");
		if (p1 != NULL) {
			int colon = ':';
			p1 = strchr(buf, colon) + 1;
			cached_size = strtoull(p1, NULL, 10);
			break;
		}
	}
	fclose(fp);
	return cached_size;
}

int MemoryProfiler::parseLine(char* line) {
	// This assumes that a digit will be found and the line ends in " Kb".
	int i = strlen(line);
	const char* p = line;
	while (*p < '0' || *p > '9')
		p++;
	line[i - 3] = '\0';
	i = atoi(p);
	return i;
}

int MemoryProfiler::countHugePages() {
	FILE* file = fopen("/proc/self/numa_maps", "r");
	if (!file) {
		return 0;
	}
	int hugepagecount = 0;
	char line[1024];
	while (fgets(line, 1024, file) != NULL) {
		//if (/huge.*dirty=(\d+)/) {
		std::string linestring(line);
		if (linestring.find("huge") != std::string::npos) {
			auto pos = linestring.find("dirty");
			if (pos != std::string::npos) {
				linestring = linestring.substr(pos);
				auto eqpos = linestring.find("=");
				auto spacepos = linestring.find(" ");
				linestring = linestring.substr(eqpos + 1, spacepos - eqpos - 1);
				hugepagecount += std::stoi(linestring);
			}
		}
	}
	fclose(file);
	return hugepagecount;
}

int MemoryProfiler::getOwnMemory() { //Note: this value is in KB!
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[1024];

	while (fgets(line, 128, file) != NULL) {
		if (strncmp(line, "VmRSS:", 6) == 0) {
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
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
