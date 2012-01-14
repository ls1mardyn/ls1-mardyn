/*
 * DynamicGeneratorFactory.cpp
 *
 * @Date: 21.09.2011
 * @Author: eckhardw
 */

#include "GeneratorFactory.h"
#include "utils/Logger.h"
#include <errno.h>
#include <iostream>
#include <cassert>


#ifdef SUPPORT_GENERATOR
#include "generators/CubicGridGenerator.h"
#include "generators/DropletGenerator.h"
#endif

using namespace Log;
using namespace std;

string GeneratorFactory::generators[] = {"droplet", "cubicgrid"};

InputBase* GeneratorFactory::loadGenerator(std::string generatorName, std::string configFile) {
#ifdef SUPPORT_GENERATOR
	MDGenerator* generator;
	if (generatorName == generators[0]) {
		generator = new DropletGenerator();
	}
	else if (generatorName == generators[1]){
		generator = new CubicGridGenerator();
	}
	else {
		global_log->error() << " Generator " << generatorName << " not known!" << endl;
	}

	// somehow, using the global_log in the generator results in a segfault
	// on the first log, somewhere in libstdc++
	//generator->setLogger(global_log);
	generator->load(configFile);

	return generator;
#else
	global_log->error() << "Generators not supported! (Compile with -DSUPPORT_GENERATORS!)" << endl;
	return NULL;
#endif
}
