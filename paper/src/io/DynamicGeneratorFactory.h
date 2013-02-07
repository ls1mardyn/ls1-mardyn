/*
 * DynamicGeneratorFactory.h
 *
 * @Date: 21.09.2011
 * @Author: eckhardw
 */

#ifndef DYNAMICGENERATORFACTORY_H_
#define DYNAMICGENERATORFACTORY_H_

#include "io/InputBase.h"
#include <vector>
#include <string>

class MDGenerator;

/**
 * This class serves as an interface to the generators in the dynamic libraries.
 * Use it to load and unload dynamic generators.
 */
class DynamicGeneratorFactory {

private:
	// hold the handle to the library.
	// I expect the factory to be called at most once per programme instance,
	// so one pair should be sufficient
#ifdef SUPPORT_DL_GENERATOR
	static std::pair<void*, MDGenerator*> handlePair;
#endif

public:

	DynamicGeneratorFactory();

	virtual ~DynamicGeneratorFactory();

	/**
	 * Construct an input generator from a dynamic library. The library has to
	 * reside in the current working directory.
	 *
	 * @param generatorName the name of the library (just ending on .so)
	 * @param configFile the name (or path) of the configuration file, from which
	 *                   the generator should be initialized.
	 */
	static InputBase* loadGenerator(std::string generatorName, std::string configFile);

	/**
	 * Destruct an input generator.
	 */
	static void unloadGenerator(InputBase* inputGenerator);

};

#endif /* DYNAMICGENERATORFACTORY_H_ */
