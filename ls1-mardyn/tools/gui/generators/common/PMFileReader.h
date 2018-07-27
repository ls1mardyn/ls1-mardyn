/*
 * PMFileReader.h
 *
 * @Date: 12.07.2011
 * @Author: eckhardw
 */

#ifndef PMFILEREADER_H_
#define PMFILEREADER_H_

#include "Generators/Generator.h"
#include "ComponentParameters.h"
#include <string>

/**
 * class to read in PM files containing definition of Mardyn components.
 */
class PMFileReader {

public:

	PMFileReader();

	virtual ~PMFileReader();

	static void readPMFile(const std::string& filename, Generator* generator, ComponentParameters* parameters);

};

#endif /* PMFILEREADER_H_ */
