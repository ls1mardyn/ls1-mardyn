/*
 * DummyExpansion.h
 *
 *  Created on: Mar 31, 2016
 *      Author: gallardjm
 */

#ifndef DUMMYEXPANSION_H_
#define DUMMYEXPANSION_H_
#include <math.h>
#include "bhfmm/fft/FFTAccelerableExpansion.h"

/*
 * Dummy SH_Expansion used to compute the transfer function from the vector
 *
 * Use a DummyStorage (truncated version of SH_Storage)
 */

/**
 * Truncated version of SH_Storage used by DummyExpansion as storage
 */
class DummyStorage {
public:
	/* constructors */
	DummyStorage(unsigned N) :
			numRows(N), totalNumEntries(N * (N + 1) / 2) {
		entries = new double[totalNumEntries];
	}

	/* destructor */
	~DummyStorage() {
		delete[] entries;
	}

	inline double & access(unsigned i, unsigned j) {
		return entries[index(i, j)];
	}
	inline const double & access_const(unsigned i, unsigned j) const {
		return entries[index(i, j)];
	}

	inline unsigned index(unsigned i, unsigned j) const {
		return i * (i + 1) / 2 + j;
	}

private:
	/** the NUMBER OF rows/columns */
	unsigned numRows;
	unsigned totalNumEntries;

	double * entries;
};

/**
 * Truncated version of SH_Expansion used by to compute the transferfunction
 * from the vector
 *
 * Extend FFTAccelerableExpansion so it can be used by a FFTAcceleration
 * evaluate_M_at_r(double X, double Y, double Z) convert the DummyExpansion to a Transferfunction
 * expansion of vector [X,Y,Z] (carthesian coordinates)
 */
class DummyExpansion: public FFTAccelerableExpansion {
public:
	/* constructors */
	DummyExpansion(unsigned ord) :
			FFTAccelerableExpansion(), order(ord), C(ord + 1), S(ord + 1) {
	}

	/* destructor */
	~DummyExpansion() {
	}

	// math operators
	void evaluate_M_at_r(double X, double Y, double Z);

	double & get_C(unsigned l, unsigned m) {
		return acc_C(l, m);
	}
	double & get_S(unsigned l, unsigned m) {
		return acc_S(l, m);
	}

private:
	//direct accessor
	inline double & acc_C(unsigned l, unsigned m) {
		return C.access(l, m);
	}
	inline double & acc_S(unsigned l, unsigned m) {
		return S.access(l, m);
	}

	//const accesssor
	inline double acc_c_C(unsigned l, unsigned m) const {
		return C.access_const(l, m);
	}
	inline double acc_c_S(unsigned l, unsigned m) const {
		return S.access_const(l, m);
	}

	// Fields:
	unsigned order;

	DummyStorage C; // C-terms: real part
	DummyStorage S; // S-terms: imaginary or (-imaginary) part, see 13.4.8
};

#endif
