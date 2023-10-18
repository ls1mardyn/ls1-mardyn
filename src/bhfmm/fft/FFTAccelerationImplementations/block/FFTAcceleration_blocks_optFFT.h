/*
 * FFTAcceleration_blocks_optFFT.h
 *
 *  Created on: Mar 08, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_BLOCK_OPTFFT_H_
#define FFTACC_BLOCK_OPTFFT_H_

#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTDataContainer_blocks.h"
#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTAcceleration_blocks.h"
#include "bhfmm/fft/tools/optimizedFFT/optFFT_API_Factory.h"
#include "bhfmm/fft/FFTSettings.h"
#include <stdexcept>

/*
 * Basic implementation using the optFFT
 *
 * Use FFTDataContainer_blocks as Data container (2 matrices)
 * initialize_target, M2L and protected function to get the DataContainer defined in abstract FFTAcceleration_blocks
 */
class FFTAcceleration_blocks_optFFT: public FFTAcceleration_blocks {

public:

	FFTAcceleration_blocks_optFFT(int order);
	~FFTAcceleration_blocks_optFFT() {
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
