/*
 * randNum.cpp
 *
 *  Created on: 20.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 *
 */

/*
 * @brief: a tiny function generating 4 digit random numbers in the range [0 , 0.9999].
 */

#include"RandomNumber.h"

//@brief: rand() is initialized by srand(time). srand() is called in main.cpp.

double RandomNumber::randNum(){
	unsigned precision = 10000;		// precision accounts for a random number within the range [0,1[ and the number of digits == precision
	int originalRandom = rand()%precision;
	double random = (double)originalRandom/precision;
	return random;
}
