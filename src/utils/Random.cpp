#include "Random.h"

Random::Random(int seed) {
	init(seed);
}

void Random::init(int seed) {
	this->ix = (seed ^ (int) 888889999) | (int) 1;
	this->iy = seed ^ (int) 777755555;
	// Calculate normalization factor
	this->am = 2.0 / (1.0 + (unsigned) ((int) -1));
}

float Random::rnd() {
	float rnd;
	const int IA = 16807;
	const int IM = 2147483647;
	const int IQ = 127773;
	const int IR = 2836;

	/* Marsaglia shift sequence */
	ix ^= (ix >> 13);
	ix ^= (ix << 17);
	ix ^= (ix >> 5);

	int k = iy / IQ;
	iy = IA * (iy - k * IQ) - IR * k;
	if (iy < 0)
		iy += IM;
	rnd = am * ((IM & (ix ^ iy)) | (int) 1);
	return rnd;
}

// by Stefan Becker
double Random::gaussDeviate(double stdDeviation) {
	/** Method generates a gaussian distributed deviate with mean = 0 and the specified standard deviation
	 * borrowed from "Numerical Recipes in C++" by W.H. Press et al., Cambridge University Press*/
	static bool iset = true;// flag decides wheter or not the reserve is returned
	static double gset; // normal distributed deviate as a reserve for next call --> speed up
	double fac, rsq, v1, v2;

	//if(randFlag < 0) iset = true; // every 2nd time the reserve value is returned
	if (iset) {
		do {
			v1 = 2.0 * this->rnd() - 1.0;
			v2 = 2.0 * this->rnd() - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		/* Box-Muller transformation to get two normal distributed deviates.
		 * One is returned and gset is saved for the next time (reserve)*/
		fac = sqrt(-2.0 * log(rsq) / rsq);

		/* perform transformation from a standard normal distributed random number v1*fac
		 *to a normal distributed number with specified standrard deviation (mean is 0) by multiplying the stdandard deviation*/
		gset = v1 * fac * stdDeviation;
		iset = false;
		return v2 * fac * stdDeviation;
	} else {
		iset = true;
		return gset;
	}
}

