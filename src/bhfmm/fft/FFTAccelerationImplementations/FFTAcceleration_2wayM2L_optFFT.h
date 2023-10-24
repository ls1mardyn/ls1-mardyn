/*
 * FFTAcceleration_2wayM2L_optFFT.h
 *
 *  Created on: Feb 10, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_2WAY_OPTFFT_H_
#define FFTACC_2WAY_OPTFFT_H_

#include "bhfmm/fft/FFTAccelerationImplementations/FFTDataContainer_arrays.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_2wayM2L.h"
#include "bhfmm/fft/tools/optimizedFFT/optFFT_API_Factory.h"
#include "bhfmm/fft/tools/fft_utils.h"

/*
 * Basic implemtation using optimized FFT
 *
 * Use FFTDataContainer_matrices as Data container (2 matrices)
 * initialize_target, M2L and protected function to get the DataContainer
 * defined in parent class
 */
class FFTAcceleration_2wayM2L_optFFT: public FFTAcceleration_2wayM2L {

public:

	FFTAcceleration_2wayM2L_optFFT(int order);
	~FFTAcceleration_2wayM2L_optFFT() {
		for (int i = 0; i < mardyn_get_max_threads(); ++i) {
			delete_matrix(Re[i]);
			delete_matrix(Im[i]);
		}
	}

	void FFT_initialize_Source(FFTAccelerableExpansion & Expansion,
			double radius);
	void FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion);

	void FFT_finalize_Target(FFTAccelerableExpansion & Expansion,
			double radius);

protected:
	//Matrices used as tmp data storage to perform optFFT (works only on matrices)
	FFT_precision*** Re;
	FFT_precision*** Im;
	optFFT_API* _optFFT_API;

	//convert matrix for the optFFT_API to array for the FFTData
	void matrix2array(FFT_precision** & m, FFT_precision* & a);
	void array2matrix(FFT_precision* & a, FFT_precision** & m);
};

#endif
