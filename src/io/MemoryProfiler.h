/*
 * MemoryProfiler.h
 *
 *  Created on: May 9, 2017
 *      Author: seckler
 */

#pragma once

#include <vector>
#include <string>

class MemoryProfiler {
public:
	MemoryProfiler();

	template<typename T>
	void registerObject(T* object);

	void doOutput(const std::string& stepInfo = std::string());

private:
	std::vector<void *> _list;

	//returns size of cached memory in kB (0 if error occurs)
	unsigned long long getCachedSize();
	void printGeneralInfo(const std::string& string);
};

