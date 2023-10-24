// Common.cpp

#include "Common.h"

#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>

int gettimestr(const char* fmt, char *buffer, int buffsize) {
	const auto rawtime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	tm unused{};
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
	return strftime(buffer, buffsize, fmt, localtime_r(&rawtime, &unused));
#pragma GCC diagnostic pop
}

std::string gettimestring(const char* fmt) {
	char buffer[80];
	gettimestr(fmt, buffer, 80);
	std::string date(buffer);

	return date;
}

/**
 * Align a number to the right by padding with a filling character.
 *
 * @param number the number value to be aligned
 * @param num_digits the number of digits
 * @param c character used for filling up
 */
std::string aligned_number(int number, int num_digits, char c) {
	std::stringstream numstream;
	numstream.fill(c);
	numstream.width(num_digits);
	numstream << number;
	std::string numstr(numstream.str());
	return numstr;
}


void calculateDistances( float * __restrict__ valuesA[3], float* __restrict__ const valuesB[3], int numValuesA, int numValuesB,
		float ** __restrict__ distances, float *** __restrict__ distanceVectors ) {

/*	for (int i = 0; i < numValuesA; i++) {
		std::cout << "A Points: " << valuesA[0][i] << "," << valuesA[1][i] << "," << valuesA[2][i] << std::endl;
	}
	for (int i = 0; i < numValuesB; i++) {
			std::cout << "B Points: " << valuesB[0][i] << "," << valuesB[1][i] << "," << valuesB[2][i] << std::endl;
	}
*/

	for (int i = 0; i < numValuesA; i++) {
		float xA = valuesA[0][i];
		float yA = valuesA[1][i];
		float zA = valuesA[2][i];
#ifdef __INTEL_COMPILER
		#pragma ivdep
#endif
		for (int j = 0; j < numValuesB; j++) {
			distanceVectors[0][i][j] = xA - valuesB[0][j];
			distanceVectors[1][i][j] = yA - valuesB[1][j];
			distanceVectors[2][i][j] = zA - valuesB[2][j];
		}
	}

	for (int i = 0; i < numValuesA; i++) {
#ifdef __INTEL_COMPILER
		#pragma ivdep
#endif
		for (int j = 0; j < numValuesB; j++) {
			distances[i][j] = distanceVectors[0][i][j] * distanceVectors[0][i][j] +
			                  distanceVectors[1][i][j] * distanceVectors[1][i][j] +
			                  distanceVectors[2][i][j] * distanceVectors[2][i][j];
		}
	}
}


/**
 * @param valuesA positions of molecules in set A: [x1, y1, z1, x2, y2, z2, ... , xN, yN, zN]; Length = 3 * numValuesA
 * @param valuesB positions of molecules in set B
 * @param distanceVectors: linearized array, innermost x,y,z, then molecules b; length = numValuesA * numValuesB * 3
 *        access distance [a][b][0] as a * numValuesB * 3 + b + 0
 * @param distances length = numValuesA * numValuesB; access distance [a][b] as a * numValuesB + b
 */
/*void calculateDistances( double valuesA[], double const valuesB[], int numValuesA, int numValuesB,
		double * __restrict__ distances, double * __restrict__ distanceVectors ) {

	for (int i = 0; i < numValuesA; i++) {
		double xA = valuesA[i];
		double yA = valuesA[i+1];
		double zA = valuesA[i+2];
		for (int j = 0; j < numValuesB; j+=3) {
			distanceVectors[i * numValuesB + j ] = xA - valuesB[j];
			distanceVectors[i * numValuesB + j + 1] = yA - valuesB[j + 1];
			distanceVectors[i * numValuesB + j + 2] = zA - valuesB[j + 2];
		}
	}

	for (int i = 0; i < numValuesA * numValuesB; i+=3) {
		distances[i] = distanceVectors[i] + distanceVectors[i];
		distances[i+1] = distanceVectors[i+1] + distanceVectors[i+1];
		distances[i+2] = distanceVectors[i+2] + distanceVectors[i+2];
	}
}*/
