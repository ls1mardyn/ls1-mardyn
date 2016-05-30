/*
 * FFTAcceleration_extensions.h
 *
 *  Created on: Mar 20, 2016
 *  Author: gallardjm
 */
#ifndef FFTACC_EXTS_H_
#define FFTACC_EXTS_H_

#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/FFTAccelerableExpansion.h"
#include "bhfmm/fft/FFTAccelerationAPI.h"

/**
 * Various additionnal M2L schemes
 * FMM algorithm will have to static_cast to it to use the order reduction in the M2L
 * All extends FFTAcceleration so FFTAcceleration's implementation can directly inherit them
 */

//!API for 2 way M2L
class FFTAccelerationAPI_2Way: public FFTAccelerationAPI {

public:

	//! destructor, child class will be upcasted, virtual destructor required to call the right child class destructor
	virtual ~FFTAccelerationAPI_2Way() {
	}

	//! M2L using 2way scheme (see doc/2wayM2L)
	virtual void FFT_M2L_2way(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction) =0;
	virtual void FFT_M2L_2way_vec(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction) =0;
};

//! API for OrderReduction scheme (extends 2way M2L)
class FFTAccelerationAPI_full: public FFTAccelerationAPI_2Way {

public:

	//! destructor, child class will be upcasted, virtual destructor required to call the right child class destructor
	virtual ~FFTAccelerationAPI_full() {
	}

	//! M2L using Order Reduction
	virtual void FFT_M2L_OrderReduction(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction, int order) =0;
	virtual void FFT_M2L_OrderReduction_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction, int order) =0;

	//M2L using 2way from inheritance

	//! M2L using both Order Reduction and 2way
	virtual void FFT_M2L_2way_ORed(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction, int order) =0;
	virtual void FFT_M2L_2way_ORed_vec(FFTAccelerableExpansion & Source1,
			FFTAccelerableExpansion & Source2,
			FFTAccelerableExpansion & Target1,
			FFTAccelerableExpansion & Target2,
			FFTDataContainer* TransferFunction, int order) =0;

};

#endif
