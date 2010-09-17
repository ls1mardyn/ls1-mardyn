/*
 * FileUtils.cpp
 *
 * @Date: 26.08.2010
 * @Author: Wolfgang Eckhardt
 */
#include "utils/FileUtils.h"
#include "utils/Logger.h"

#include <sys/stat.h>
#include <cstdio>

using namespace Log;

bool fileExists(const char* fileName) {
	struct stat status;
	int retVal = stat(fileName, &status);

	return retVal == 0;
}


void removeFile(const char* fileName) {
	int retVal = remove(fileName);
	if (retVal != 0) {
		global_log->warning() << "Could not remove file " << fileName << std::endl;
	}
}

unsigned int getFileSize(const char* fileName) {
	struct stat status;
	int stat_status = stat(fileName, &status);
	unsigned int retVal = 0;
	if (stat_status == 0) {
		retVal = (unsigned int) status.st_size;
	} else {
		global_log->warning() << "Could not stat file " << fileName << std::endl;
	}
	return retVal;
}
