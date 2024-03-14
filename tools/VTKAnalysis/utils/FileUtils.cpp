/*
 * FileUtils.cpp
 *
 * @Date: 26.08.2010
 * @Author: Wolfgang Eckhardt
 */
#include "FileUtils.h"

#include <sys/stat.h>
#include <errno.h>
#include <cstring>
#include <cstdio>
#include <iostream>

bool fileExists(const char* fileName) {
	struct stat status;
	int retVal = stat(fileName, &status);
	if (retVal != 0) {
		char* error = strerror(errno);
		std::cerr << "File does not exist: " << fileName << std::endl;
		std::cerr << "Error: " << error << std::endl;
	}
	return retVal == 0;
}

std::string getFileExtension(const char* fileName) {
	std::string extension ("");
	const char* pos = strrchr(fileName, '.');
	if(NULL != pos) {
		extension = std::string(pos + 1);
	}
	return extension;
}

bool isAbsolutePath(const std::string& path) {
	return path.find_first_of("/~")==0;
}

void removeFile(const char* fileName) {
	int retVal = remove(fileName);
	if (retVal != 0) {
		char* error = strerror(errno);
		std::cerr << "Could not remove file " << fileName << std::endl;
        std::cerr << "Error: " << error << std::endl;
	}
}

unsigned int getFileSize(const char* fileName) {
	struct stat status;
	int stat_status = stat(fileName, &status);
	unsigned int retVal = 0;
	if (stat_status == 0) {
		retVal = (unsigned int) status.st_size;
	} else {
		std::cerr << "Could not stat file " << fileName << std::endl;
	}
	return retVal;
}

std::ostream& operator<<( std::ostream& o, const fill_width& a )
{
	o.fill( a.fill );
	o.width( a.width );
	return o;
}
