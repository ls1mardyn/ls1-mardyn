/*
 * FFTAcceleration_scalBlocks_optFFT.h
 *
 *  Created on: Mar 16, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_SCALBLOCKS_OPTFFT_H_
#define FFTACC_SCALBLOCKS_OPTFFT_H_

#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTDataContainer_scalBlocks.h"
#include "bhfmm/fft/FFTAccelerationAPI_extensions.h"
#include "bhfmm/fft/FFTSettings_preprocessor.h"
#include "bhfmm/fft/FFTSettings.h"
#include "bhfmm/fft/tools/optimizedFFT/optFFT_API_Factory.h"

/*
 * Use FFTDataContainer_scalBlocks as Data container (2 arrays of matrices)
 */
class FFTAcceleration_scalBlocks_optFFT: public FFTAccelerationAPI_full {

public:

	FFTAcceleration_scalBlocks_optFFT(int order);
	~FFTAcceleration_scalBlocks_optFFT() {
		delete _optFFT_API;
		delete[] _blockSize;
        for (int i = 0; i < mardyn_get_max_threads(); ++i) {
            delete_matrix(_Re_tmp[i]);
            delete_matrix(_Im_tmp[i]);
        }
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

	//From FFTAcceleration_full
	void FFT_M2L_OrderReduction(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction, int order);
	void FFT_M2L_OrderReduction_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction, int order);

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

	void FFT_M2L_2way_ORed(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction, int order);
	void FFT_M2L_2way_ORed_vec(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction, int order);

protected:
	FFTDataContainer_scalBlocks* getFFTData(
			FFTAccelerableExpansion & Expansion);
	FFTDataContainer_scalBlocks* getFFTData_scal(
			FFTAccelerableExpansion & Expansion);
	optFFT_API* _optFFT_API;

	int _nbBlocks;
	int _nbLinePerBlock; //number of usefull line per block
	int* _blockSize;     //fft_ny size of each block
    //tmp matrices for FFT for each thread
	FFT_precision*** _Re_tmp;
	FFT_precision*** _Im_tmp;

	template<bool Vect, bool OrderRed>
	void FFT_M2L_template(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction, int order);
	template<bool Vect, bool OrderRed>
	void FFT_M2L_2way_template(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction, int order);

};

#endif
