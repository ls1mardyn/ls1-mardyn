/*
 * optFFT_API_Factory.h
 *
 *  Created on: Mar 06, 2016
 *      Author: gallardjm
 */
#ifndef OPTFFT_API_FACTORY_H_
#define OPTFFT_API_FACTORY_H_

#include "optFFT_API.h"
#include "FakedOptFFT.h"
#include "BasicOptFFT.h"

class optFFT_API_Factory {
public:

	/**
	 * Return an OptFFT_API
	 * BasicOptFFT is an API for Kurzak's code, it require a square abstracted matrix
	 * (=> input matrix have size_x = 2*size_y) and an order between 5 and 15
	 * FakedOptFFT use the same scheme as a true optFFT but perform the FFT and IFFT
	 * with FFTW by reconstructing the whole matrix and using a default FFT alg., it
	 * works for any matrix.
	 *
	 * @param const int  order, order of the expansion
	 * @param const bool square_matrix, == (size_x == 2*size_y)
	 * @return optFFT_API* pointer to an implementation of the OptFFT_API
	 */
	static optFFT_API* getOptFFT_API(const int order,
			const bool square_matrix) {
		if (square_matrix && order > 4 && order < 16)
			return new BasicOptFFT();
		else
			return new FakedOptFFT();
	}
};

#endif
