/*
 * FFTAcceleration_2wayM2L_fftw.h
 *
 *  Created on: Feb 10, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_2WAY_FFTW_H_
#define FFTACC_2WAY_FFTW_H_

#include "bhfmm/fft/FFTAccelerationImplementations/FFTDataContainer_arrays.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_2wayM2L.h"
#include "bhfmm/fft/tools/FFTW_Helper.h"
#include "bhfmm/fft/tools/fft_utils.h"

/*
 * Basic implemtation using the FFTW library
 *
 * Use FFTDataContainer_arrays as Data container (2 arrays, here aligned)
 * initialize_target, M2L and protected function to get the DataContainer defined in abstract FFTAcceleration_2wayM2L
 */
class FFTAcceleration_2wayM2L_fftw: public FFTAcceleration_2wayM2L {

public:

	FFTAcceleration_2wayM2L_fftw(int order);
	~FFTAcceleration_2wayM2L_fftw() {
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
