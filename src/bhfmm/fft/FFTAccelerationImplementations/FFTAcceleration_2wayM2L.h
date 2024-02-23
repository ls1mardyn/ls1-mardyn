/*
 * FFTAcceleration_2wayM2L.h
 *
 *  Created on: Feb 10, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_ARR_H_
#define FFTACC_ARR_H_

#include "bhfmm/fft/FFTAccelerationImplementations/FFTDataContainer_arrays.h"
#include "bhfmm/fft/FFTAccelerationAPI_extensions.h"
#include "bhfmm/fft/FFTSettings_preprocessor.h"
#include "bhfmm/fft/FFTSettings.h"

/**
 * Abstract class with common code for basic fft and fftw implementation
 *
 * Use FFTDataContainer_arrays as Data container (2 arrays)
 *
 * Don't define FFT_M2L and FFT_M2L_vec, but define FFT_M2L_2way
 * the code using the acceleration should downcast its FFTAcceleration
 * to the FFTAcceleration_2Way class and use FFT_M2L_2way instead of FFT_M2L
 */
class FFTAcceleration_2wayM2L: public FFTAccelerationAPI_2Way {

public:

	//child class will be upcasted, virtual destructor required to call the right child class destructor
	virtual ~FFTAcceleration_2wayM2L() {
	}

	void FFT_initialize_Target(FFTAccelerableExpansion & Expansion);

	//Not defined since they shouldn't be used
	void FFT_M2L(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction) {
	}
	void FFT_M2L_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction) {
	}

	//2way M2L mathods from FFTAcceleration_2Way
	void FFT_M2L_2way(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction);
	void FFT_M2L_2way_vec(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction);

protected:
	FFTDataContainer_arrays* getFFTData(FFTAccelerableExpansion & Expansion);

	int _nbZeroes; // number of zeroes used for the padding (see bhfmm/fft/tools/fft_utils)
	int _totalSize; //_fft_nx * _fft_ny + _nbZeroes * 2, is even as _fft_nx must be even
};

#endif
