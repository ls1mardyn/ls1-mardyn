/*
 * DynamicGeneratorFactory.h
 *
 * @Date: 21.09.2011
 * @Author: eckhardw
 */

#ifndef GENERATORFACTORY_H_
#define GENERATORFACTORY_H_

#include <string>

#include "io/InputBase.h"

/**
 * This class serves as an interface to the generators in the dynamic libraries.
 * Use it to load and unload dynamic generators.
 */
class GeneratorFactory {

private:
	GeneratorFactory() {}

	virtual ~GeneratorFactory() {}

	static std::string generators[];
public:

	/**
	 * Construct an input generator.
	 *
	 * @param generatorName
	 * @param configFile the name (or path) of the configuration file, from which
	 *                   the generator should be initialized.
	 */
	static InputBase* loadGenerator(std::string generatorName, std::string configFile);

};

#endif /* GENERATORFACTORY_H_ */
