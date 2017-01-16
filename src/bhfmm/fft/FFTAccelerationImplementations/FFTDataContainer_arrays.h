/*
 * FFTDataContainer_arrays.h
 *
 *  Created on: Feb 09, 2016
 *      Author: gallardjm
 */
#ifndef FFTDATA_ARR_H_
#define FFTDATA_ARR_H_

#include "bhfmm/fft/FFTSettings_preprocessor.h" //tmp include for the typedef FFT_precision
#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/tools/fft_utils.h"

class FFTDataContainer_arrays: public FFTDataContainer {

public:
	FFT_precision* Re;
	FFT_precision* Im;
	int _size;

	FFTDataContainer_arrays(int size) :
			_size(size) {
	}

	~FFTDataContainer_arrays() {
		delete_aligned_array(Re);
		delete_aligned_array(Im);
	}

	FFTDataContainer* copyContainer() {
		FFTDataContainer_arrays* copy = new FFTDataContainer_arrays(_size);
		copy->Re = copy_aligned_array(Re, _size);
		copy->Im = copy_aligned_array(Im, _size);

		return copy;
	}
};

#endif
