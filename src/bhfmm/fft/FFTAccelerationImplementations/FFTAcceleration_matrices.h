/*
 * FFTAcceleration_matrices.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_MAT_H_
#define FFTACC_MAT_H_

#include "bhfmm/fft/FFTAccelerationImplementations/FFTDataContainer_matrices.h"
#include "bhfmm/fft/FFTAccelerationAPI.h"
#include "bhfmm/fft/FFTSettings_preprocessor.h"
#include "bhfmm/fft/FFTSettings.h"

/**
 * Abstract class with common code for basic fft and fftw implementation
 *
 * Use FFTDataContainer_matrices as Data container (2 matrices)
 */
class FFTAcceleration_matrices: public FFTAccelerationAPI {

public:

	//child class will be upcasted, virtual destructor required to call the right child class destructor
	virtual ~FFTAcceleration_matrices() {
	}

	void FFT_initialize_Target(FFTAccelerableExpansion & Expansion);

	void FFT_M2L(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction);
	void FFT_M2L_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction);

protected:
	FFTDataContainer_matrices* getFFTData(FFTAccelerableExpansion & Expansion);
};

#endif
