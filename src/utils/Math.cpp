/*
 * Math.cpp
 *
 *  Created on: 9 Dec 2024
 *      Author: amartyads
 */

#include "Math.h"

bool isNearRel(double a, double b, double maxRelativeDifference) 
{
	const auto greaterNumber = std::max(std::abs(a), std::abs(b));
	const auto absoluteDifference = maxRelativeDifference * greaterNumber;
	const auto diff = std::abs(a - b);
	return (diff <= absoluteDifference);
}

short int findSign(int n) { return n < 0 ? -1 : 1; }
