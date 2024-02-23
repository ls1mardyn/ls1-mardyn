/*
 * FFTAcceleration_blocks.h
 *
 *  Created on: Feb 12, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_BLOCK_H_
#define FFTACC_BLOCK_H_

#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTDataContainer_blocks.h"
#include "bhfmm/fft/FFTAccelerationAPI_extensions.h"
#include "bhfmm/fft/FFTSettings_preprocessor.h"
#include "bhfmm/fft/FFTSettings.h"

/*
 * Abstract class with common code for basic fft and fftw implementation
 *
 * Use FFTDataContainer_blocks as Data container (2 arrays of matrices)
 */
class FFTAcceleration_blocks: public FFTAccelerationAPI_full {

public:

	//child class will be upcasted, virtual destructor required to call the right child class destructor
	virtual ~FFTAcceleration_blocks() {
	}

	void FFT_initialize_Target(FFTAccelerableExpansion & Expansion);

	void FFT_M2L(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction);
	void FFT_M2L_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction);

	//From FFTAcceleration_OrderReduction
	void FFT_M2L_OrderReduction(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction, int order);
	void FFT_M2L_OrderReduction_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction, int order);

	//2way M2L, require _fft_nx = 2 * _p
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
	FFTDataContainer_blocks* getFFTData(FFTAccelerableExpansion & Expansion);

	int _nbBlocks;
	int _nbLinePerBlock; //number of usefull line per block

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
