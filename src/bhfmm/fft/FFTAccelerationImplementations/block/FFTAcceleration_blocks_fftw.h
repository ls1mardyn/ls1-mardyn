/*
 * FFTAcceleration_blocks_fftw.h
 *
 *  Created on: Feb 12, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_BLOCK_FFTW_H_
#define FFTACC_BLOCK_FFTW_H_

#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTDataContainer_blocks.h"
#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTAcceleration_blocks.h"
#include "bhfmm/fft/tools/FFTW_API.h"
#include <stdexcept>

/*
 * Basic implementation using the FFTW library
 *
 * Use FFTDataContainer_blocks as Data container (2 matrices)
 * initialize_target, M2L and protected function to get the DataContainer defined in abstract FFTAcceleration_blocks
 */
class FFTAcceleration_blocks_fftw: public FFTAcceleration_blocks {

public:

	FFTAcceleration_blocks_fftw(int order);
	~FFTAcceleration_blocks_fftw() {
		delete _fftw_api;
	}

	void FFT_initialize_Source(FFTAccelerableExpansion & Expansion,
			double radius);
	void FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion);

	void FFT_finalize_Target(FFTAccelerableExpansion & Expansion,
			double radius);

protected:
	FFTW_API* _fftw_api;

};

#endif
