/*
 * MemoryProfiler.h
 *
 *  Created on: May 9, 2017
 *      Author: seckler
 */

#pragma once

#include <cstddef>
#include <vector>
#include <string>

class MemoryProfilable {
public:
	virtual ~MemoryProfilable() {
	}
	virtual size_t getTotalSize() = 0;
	virtual void printSubInfo(int offset) = 0;
	virtual std::string getName() = 0;
};

class MemoryProfiler {
private:
	std::vector<MemoryProfilable **> _list;
	int _hugePageSize; // Hugepagesize in kB
public:
	MemoryProfiler();

	void registerObject(MemoryProfilable** object);

	void doOutput(const std::string& stepInfo = std::string());

private:
	//returns size of cached memory in kB (0 if error occurs)
	unsigned long long getCachedSize();
	void printGeneralInfo(const std::string& string);
	int parseLine(char* line);
	int countHugePages();
	int getOwnMemory();
};

