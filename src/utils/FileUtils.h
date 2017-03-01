/*
 * FileUtils.h
 *
 * @Date: 26.08.2010
 * @Author: Wolfgang Eckhardt
 */

#ifndef FILEUTILS_H_
#define FILEUTILS_H_

#include <string>

/**
 * Check if a file exists.
 */
bool fileExists(const char* fileName);

/*
 * Get the file name extension marked by the last '.' in the filename.
 */
std::string getFileExtension(const char* fileName);

/**
 * Delete a file from the system.
 */
void removeFile(const char* fileName);

/**
 * Retrieve the size of a file.
 *
 * @return 0 if an error occurred retrieving the size, so the user has to check
 *         if the file exists.
 */
unsigned int getFileSize(const char* fileName);

/**
 * Adding string of numbers with leading zeros to stream (e.g. simstep in output filename)
 */
struct fill_width
{
	fill_width( char f, uint8_t w )
		: fill(f), width(w) {}
	char fill;
	int width;
};

std::ostream& operator<<( std::ostream& o, const fill_width& a );

#endif /* FILEUTILS_H_ */
