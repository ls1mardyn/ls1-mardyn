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

MemoryProfiler::MemoryProfiler() {

}

template<typename T>
void MemoryProfiler::registerObject(T* object) {
	_list.push_back(static_cast<void*>(object));
}

void MemoryProfiler::doOutput(const std::string& string) {
	printGeneralInfo(string);

	// further info
}

unsigned long long MemoryProfiler::getCachedSize() {
	size_t MAXLEN = 1024;
	FILE *fp;
	char buf[MAXLEN];
	fp = fopen("/proc/meminfo", "r");
	while (fgets(buf, MAXLEN, fp)) {
		char *p1 = strstr(buf, "Cached:");
		if (p1 != NULL) {
			int colon = ':';
			char *p1 = strchr(buf, colon) + 1;
			//std::cout << p1 << endl;
			unsigned long long t = strtoull(p1, NULL, 10);
			//std::cout << t << endl;
			return t;
		}
	}
	return 0;
}

int parseLine(char* line) {
	// This assumes that a digit will be found and the line ends in " Kb".
	int i = strlen(line);
	const char* p = line;
	while (*p < '0' || *p > '9')
		p++;
	line[i - 3] = '\0';
	i = atoi(p);
	return i;
}

int getOwnMemory() { //Note: this value is in KB!
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL) {
		if (strncmp(line, "VmRSS:", 6) == 0) {
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

void MemoryProfiler::printGeneralInfo(const std::string& string) {
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	long long totalMem = memInfo.totalram * memInfo.mem_unit / 1024 / 1024;
	long long usedMem = ((memInfo.totalram - memInfo.freeram - memInfo.bufferram) * memInfo.mem_unit / 1024
			- getCachedSize()) / 1024;
	std::stringstream additionalinfo;
	if (string.length() > 0) {
		additionalinfo << " (" << string << ")";
	}
	additionalinfo << ":" << std::endl;
	Log::global_log->info() << "Memory consumption" << additionalinfo.str();
	Log::global_log->info() << "\tMemory usage (System total):\t" << usedMem << " MB out of " << totalMem << " MB ("
			<< usedMem * 100. / totalMem << "%)" << std::endl;
	double ownMem = getOwnMemory() * 1.e-3;
	Log::global_log->info() << "\tBy own process:\t\t\t" << ownMem << " MB (" << ownMem * 100. / totalMem
			<< "% of total memory)" << std::endl;
}
