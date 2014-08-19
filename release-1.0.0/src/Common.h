// Common.h

#ifndef COMMON_H_
#define COMMON_H_

#include <string>

class Common;

//                                        "%Y-%m-%d_%H-%M-%S"
std::string gettimestring(const char* fmt = "%y%m%dT%H%M%S");

//! @brief Helper returning an aligned number.
//!
std::string aligned_number(int number, int num_digits, char c = ' ');

/**
 * Calculate the full distance Matrix as well as the distance vectors
 * between a set of points A and B.
 *
 * @param valuesA Array of three arrays, containing x, then y, then z coordinates,
 *        i.e. [x1, x2, ... xN], [y1, y2, ... yN], [z1, z2, ... zN]
 * @param valuesB see valuesA
 * @param numValuesA size of point set A (indicating the length of the array valuesA[i])
 * @param numValuesB size of point set B
 * @param distances twodimensional array of length [A][B], to which the distances will be stored
 * @param distanceVectors: containing the distance vectors, with size 3 x numValuesA x numValuesB
 */
void calculateDistances( float * __restrict__ valuesA[3], float* __restrict__ const valuesB[3], int numValuesA, int numValuesB,
		float ** __restrict__ distances, float *** __restrict__ distanceVectors );

/**
 * @param valuesA positions of molecules in set A: [x1, y1, z1, x2, y2, z2, ... , xN, yN, zN]; Length = 3 * numValuesA
 * @param valuesB positions of molecules in set B
 * @param distanceVectors: linearized array, innermost x,y,z, then molecules b; length = numValuesA * numValuesB * 3
 *        access distance [a][b][0] as a * numValuesB * 3 + b + 0
 * @param distances length = numValuesA * numValuesB; access distance [a][b] as a * numValuesB + b
 */
/*void calculateDistances( double valuesA[], double const valuesB[], int numValuesA, int numValuesB,
		double * __restrict__ distances, double * __restrict__ distanceVectors ); */

#endif /*COMMON_H_*/
