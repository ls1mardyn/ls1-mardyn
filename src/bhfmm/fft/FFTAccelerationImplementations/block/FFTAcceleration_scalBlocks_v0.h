/*
 * FFTAcceleration_scalBlocks_v0.h
 *
 *  Created on: Mar 15, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_SCALBLOCKS_V0_H_
#define FFTACC_SCALBLOCKS_V0_H_

#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTDataContainer_scalBlocks_v0.h"
#include "bhfmm/fft/FFTAccelerationAPI.h"
#include "bhfmm/fft/FFTSettings_preprocessor.h"
#include "bhfmm/fft/FFTSettings.h"
#include "bhfmm/fft/tools/optimizedFFT/optFFT_API_Factory.h"

/*
 * v0 version of scaling block
 *
 * Use FFTDataContainer_scalBlocks_v0 as Data container (2 arrays of matrices)
 */
class FFTAcceleration_scalBlocks_v0: public FFTAccelerationAPI {

public:

	FFTAcceleration_scalBlocks_v0(int order);
	~FFTAcceleration_scalBlocks_v0() {
		delete _optFFT_API;
		delete[] _blockSize;
	}

	void FFT_initialize_Target(FFTAccelerableExpansion & Expansion);
	void FFT_initialize_Source(FFTAccelerableExpansion & Expansion,
			double radius);
	void FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion);

	void FFT_finalize_Target(FFTAccelerableExpansion & Expansion,
			double radius);

	void FFT_M2L(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction);
	void FFT_M2L_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction);

protected:
	FFTDataContainer_scalBlocks_v0* getFFTData(
			FFTAccelerableExpansion & Expansion);
	FFTDataContainer_scalBlocks_v0* getFFTData_scal(
			FFTAccelerableExpansion & Expansion);
	optFFT_API* _optFFT_API;

	int _nbBlocks;
	int _nbLinePerBlock; //number of usefull line per block
	int* _blockSize;     //fft_ny size of each block

};

#endif
