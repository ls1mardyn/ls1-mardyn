/*
 * FFTAcceleration_matrices_optFFT.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_MAT_OPTFFT_H_
#define FFTACC_MAT_OPTFFT_H_

#include "bhfmm/fft/FFTAccelerationImplementations/FFTDataContainer_matrices.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_matrices.h"
#include "bhfmm/fft/tools/optimizedFFT/optFFT_API_Factory.h"

/*
 * Basic implemtation using Kurzak's optimized FFT
 *
 * Use FFTDataContainer_matrices as Data container (2 matrices)
 * initialize_target, M2L and protected function to get the DataContainer defined in abstract FFTAcceleration_matrices
 */
class FFTAcceleration_matrices_optFFT: public FFTAcceleration_matrices {

public:

	FFTAcceleration_matrices_optFFT(int order);
	~FFTAcceleration_matrices_optFFT() {
		delete _optFFT_API;
	}

	void FFT_initialize_Source(FFTAccelerableExpansion & Expansion,
			double radius);
	void FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion);

	void FFT_finalize_Target(FFTAccelerableExpansion & Expansion,
			double radius);

protected:
	optFFT_API* _optFFT_API;
};

#endif
