/*
 * MardynConfigLegacyWriter.h
 *
 * @Date: 13.09.2011
 * @Author: eckhardw
 */

#ifndef MARDYNCONFIGLEGACYWRITER_H_
#define MARDYNCONFIGLEGACYWRITER_H_

#include <string>

class MardynConfiguration;

class MardynConfigLegacyWriter {

public:

	MardynConfigLegacyWriter();

	virtual ~MardynConfigLegacyWriter();

	static void writeConfigFile(const std::string& directory, const std::string& fileName,
			const MardynConfiguration& config);
};

#endif /* MARDYNCONFIGLEGACYWRITER_H_ */
