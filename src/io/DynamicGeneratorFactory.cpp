/*
 * DynamicGeneratorFactory.cpp
 *
 * @Date: 21.09.2011
 * @Author: eckhardw
 */

#include "DynamicGeneratorFactory.h"
#include "utils/Logger.h"
//#include <sys/types.h>
//#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <dlfcn.h>
#include <cassert>


#ifdef SUPPORT_DL_GENERATOR
#include "generators/MDGenerator.h"
#include "src/Parameters/Parameter.h"
#include "src/Parameters/ParameterCollection.h"
#include "src/IO/WriteOutput.h"
#endif

using namespace Log;
using namespace std;

#ifdef SUPPORT_DL_GENERATOR
std::pair<void*, MDGenerator*> DynamicGeneratorFactory::handlePair = make_pair((void*) NULL, (MDGenerator*) NULL);
#endif

DynamicGeneratorFactory::DynamicGeneratorFactory() {
}

DynamicGeneratorFactory::~DynamicGeneratorFactory() {
}

InputBase* DynamicGeneratorFactory::loadGenerator(std::string generatorName, std::string configFile) {
#ifdef SUPPORT_DL_GENERATOR
	if (handlePair.second != NULL) {
		global_log->error() << "Could not load library " << generatorName << ": " << endl;
		global_log->error() << "A Generator (" << handlePair.second->getName() << ") has already been created!" << endl;
		exit(-1);
	}

	void* handle = dlopen(generatorName.c_str(), RTLD_NOW);
	if (!handle) {
		global_log->error() << "Could not load library " << generatorName << ": " << dlerror() << endl;
		exit(-1);
	}

	Generator* (*factory)();
	factory = (Generator* (*)()) dlsym(handle, "create_generator");
	if (const char* error = dlerror()) {
		global_log->error() << "Could not load factory function for " << generatorName << ":" << error << endl;
		exit(-1);
	}

	MDGenerator* generator = static_cast<MDGenerator*> (factory());
	std::pair<void*,MDGenerator*> pair = make_pair(handle, generator);
	handlePair = pair;
	// somehow, using the global_log in the generator results in a segfault
	// on the first log, somewhere in libstdc++
	// generator->setLogger(global_log);
	generator->load(configFile);

	return generator;
#endif
}


void DynamicGeneratorFactory::unloadGenerator(InputBase* inputGenerator) {
#ifdef SUPPORT_DL_GENERATOR
	assert(reinterpret_cast<Generator*>(inputGenerator) == handlePair.second);

	void* handle = handlePair.first;
	Generator* generator = handlePair.second;

	void (*destructor)(Generator*);
	destructor = (void (*)(Generator*)) dlsym(handle, "destruct_generator");
	if (const char* error = dlerror()) {
		global_log->error() << "Could not load destructor function for " << generator->getName() << ":" << error << endl;
	}

	destructor(generator);
	if (dlclose(handle)) {
		global_log->error() << "Error unloading library: " << dlerror() << endl;
	}
#endif
}
