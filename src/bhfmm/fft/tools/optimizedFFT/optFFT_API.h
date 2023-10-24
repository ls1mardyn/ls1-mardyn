/*
 * optFFT_API.h
 *
 *  Created on: Mar 06, 2015
 *      Author: gallardjm
 */

#ifndef OPTFFT_API_H_
#define OPTFFT_API_H_

#include "bhfmm/fft/FFTSettings_preprocessor.h"

/**
 * API of the optimized FFT
 *
 * Perform all operations in place
 */
class optFFT_API {

public:

	/**
	 * Perform an optimized FFT using the symmetry of the expansions
	 *
	 * Assume that
	 *  the matrices are the real and imag part of the expansion
	 *  for Expansion(l,m), l is the first coordinate of the matrices
	 *  The symmetry Expansion(l,-m) = (-1)^m conj(Expansion(l,m)) is valid
	 *
	 * The resulting FFT is done on size_x,2*size_y matrices that would represent on the line i
	 * Ex(i,0),.., Ex(i,size_x), | 0, Ex(i,-size_x),.., Ex(i,-1)
	 * (the part after | is reconstruct for the FFT and doesn't appear outside
	 * of it or in the result, Ex(i,j)=0 if j<i or j>i)
	 *
	 * The output is the first half of the reconstructed matrix since it verifies
	 * M(i,j) = conj(M(i,j+size_y))
	 * The output is safe to use for entrywise default operation with other
	 * matrices converted with this function
	 *
	 * @param FFT_precision** & Real, real part
	 * @param FFT_precision** & Imag, imag part
	 * @param const int size_x, size of the first coord of the matrices
	 * @param const int size_y, size of the second coord of the matrices
	 */
	virtual void optimizedFFT(FFT_precision** & Real, FFT_precision** & Imag,
			const int size_x, const int size_y) =0;

	/** Inverse FFT
	 *  Result is assumed to be rescaled
	 *
	 * Take a matrix couples converted with optimizedFFT as input and assume the symmetry M(i,j) = conj(M(i,j+size_y))
	 * on the reconstructed full matrix
	 *
	 * Output the first half of the reconstructed matrix IFFT
	 */
	virtual void optimizedIFFT(FFT_precision** & Real, FFT_precision** & Imag,
			const int size_x, const int size_y) =0;

	//! virtual destructor for inheritance
	virtual ~optFFT_API() {
	}
	;
};

#endif
