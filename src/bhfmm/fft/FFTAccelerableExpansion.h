/*
 * FFTAccelerableExpansion.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTACCELERABLEEXP_H_
#define FFTACCELERABLEEXP_H_

#include "bhfmm/fft/FFTDataContainer.h"
#include <cstddef>

/**
 * Interface of the expansion used by the FFTAcceleration,
 * Expansions have to inherit it
 *
 * Requires accessor to the expansion's values and set a FFTDataContainer
 * storage (abstract class by itself)
 * Delete the storage by itself in its virtual destructor.
 */
class FFTAccelerableExpansion {

public:

	//!Constructor, set the storage pointer to NULL by default
	FFTAccelerableExpansion() :
			_FFTData(NULL) {
	}
	//!virtual destructor to ensure correct destructor sequence, delete the _FFTData if set
	virtual ~FFTAccelerableExpansion() {
		if (issetFFTData())
			delete _FFTData;
	}

	//Accessor to the expansions, to be defined by the expansion implementation
	/**
	 * Accessor to the real part of an expansion.
	 * Only require access to the positive m part,
	 * the negative m will be reconstructed using the symmetry
	 *
	 * @param unsigned int l, 0 <= l < max_order
	 * @param unsigned int m, 0 <= m <= l
	 * @return double & C[l][m]
	 */
	virtual double & get_C(unsigned l, unsigned m) =0;

	/**
	 * Accessor to the imaginary part of an expansion.
	 * Only require access to the positive m part,
	 * the negative m will be reconstructed using the symmetry
	 *
	 * @param unsigned int l, 0 <= l < max_order
	 * @param unsigned int m, 0 <= m <= l
	 * @return double & C[l][m]
	 */
	virtual double & get_S(unsigned l, unsigned m) =0;

	/**
	 * FFTDataContainer pointer, abstract class whose implementations
	 * store all required  data to the M2L FFT acceleration by
	 * an FFTAcceleration's implementation
	 *
	 * Will need to be copied by a copy constructor using _FFTData->copyContainer()
	 * to get a copy (see FFTDataContainer.h)
	 */
	FFTDataContainer* _FFTData;

	//! inline isset test
	inline bool issetFFTData() {
		return (_FFTData != NULL);
	}
};

#endif
