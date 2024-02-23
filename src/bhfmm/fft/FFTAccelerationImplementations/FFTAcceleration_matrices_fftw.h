/*
 * FFTAcceleration_matrices_fftw.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_MAT_FFTW_H_
#define FFTACC_MAT_FFTW_H_

#include "bhfmm/fft/FFTAccelerationImplementations/FFTDataContainer_matrices.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_matrices.h"
#include "bhfmm/fft/tools/FFTW_Helper.h"

/*
 * Basic implemtation using the FFTW library
 *
 * Use FFTDataContainer_matrices as Data container (2 matrices)
 * initialize_target, M2L and protected function to get the DataContainer defined in abstract FFTAcceleration_matrices
 */
class FFTAcceleration_matrices_fftw: public FFTAcceleration_matrices {

public:

	FFTAcceleration_matrices_fftw(int order);
	~FFTAcceleration_matrices_fftw() {
		delete _fftw_h;
	}

	void FFT_initialize_Source(FFTAccelerableExpansion & Expansion,
			double radius);
	void FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion);

	void FFT_finalize_Target(FFTAccelerableExpansion & Expansion,
			double radius);

protected:
	FFTW_Helper* _fftw_h;

};

#endif
