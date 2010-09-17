/*
 * FileUtils.h
 *
 * @Date: 26.08.2010
 * @Author: Wolfgang Eckhardt
 */

#ifndef FILEUTILS_H_
#define FILEUTILS_H_

/**
 * Check if a file exists.
 */
bool fileExists(const char* fileName);

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

#endif /* FILEUTILS_H_ */
