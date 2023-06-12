/*
 * FFTOrderReduction.h
 *
 *  Created on: Mar 20, 2016
 *  Author: gallardjm
 */
#ifndef FFTORDERRED_H_
#define FFTORDERRED_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**
 * Static class used by the Order Reduction scheme
 * Provides a minimum M2L local order to use for this specific M2L operation
 * using the distance between the source and the target and the expansions' order
 */
class FFTOrderReduction {

public:

	// return the max p (=order+1) that the expansion should use -> for 2,0,0,15 => 16
	static int getM2LOrder(int x, int y, int z, int order);

private:
	static int _ReducedOrder[7];
	static int _order;

	static void initialize(int order);

};

#endif
