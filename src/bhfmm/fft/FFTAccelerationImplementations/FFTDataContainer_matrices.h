/*
 * FFTDataContainer_matrices.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTDATA_MAT_H_
#define FFTDATA_MAT_H_

#include "bhfmm/fft/FFTSettings_preprocessor.h" //tmp include for the typedef FFT_precision
#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/tools/fft_utils.h"

class FFTDataContainer_matrices: public FFTDataContainer {
public:
	FFT_precision** Re;
	FFT_precision** Im;
	int _nx;
	int _ny;

	FFTDataContainer_matrices(int nx, int ny) :
			_nx(nx), _ny(ny) {
	}

	~FFTDataContainer_matrices() {
		delete_matrix(Re);
		delete_matrix(Im);
	}

	FFTDataContainer* copyContainer() {
		FFTDataContainer_matrices* copy = new FFTDataContainer_matrices(_nx,
				_ny);
		copy->Re = copy_matrix(Re, _nx, _ny);
		copy->Im = copy_matrix(Im, _nx, _ny);

		return copy;
	}
};

#endif
