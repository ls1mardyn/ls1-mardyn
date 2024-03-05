/*
 * FFTFactory.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTFACTORY_H_
#define FFTFACTORY_H_

#include "bhfmm/fft/FFTSettings.h"
#include "bhfmm/fft/FFTAccelerationAPI.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_matrices_optFFT.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_matrices_fftw.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_2wayM2L_fftw.h"
#include "bhfmm/fft/FFTAccelerationImplementations/FFTAcceleration_2wayM2L_optFFT.h"
#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTAcceleration_blocks_fftw.h"
#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTAcceleration_blocks_optFFT.h"
#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTAcceleration_scalBlocks_v0.h"
#include "bhfmm/fft/FFTAccelerationImplementations/block/FFTAcceleration_scalBlocks_optFFT.h"
#include "bhfmm/fft/TransferFunctionManagerAPI.h"
#include "bhfmm/fft/transferFunctionManager/TransferFunctionManager.h"
#include "bhfmm/fft/transferFunctionManager/TransferFunctionManager_UniformGrid.h"
#include <stdexcept>

/**
 * Factory static class to get an FFTAcceleration and a TransferFunctionManager
 * based on the current settings defined in the static class FFTSettings
 */
class FFTFactory {

public:

	static FFTAccelerationAPI* getFFTAccelerationAPI(int order) {
		if (FFTSettings::USE_BLOCK) {
			if (FFTSettings::USE_FFTW) {
				if (FFTSettings::USE_ADVBLOCK)
					throw std::invalid_argument(
							"No advblock fftw implementation, use fft instead");
				else
					return new FFTAcceleration_blocks_fftw(order);
			} else {
				if (FFTSettings::USE_ADVBLOCK)
					return new FFTAcceleration_scalBlocks_optFFT(order); //return new FFTAcceleration_scalBlocks_v0(order);
				else
					return new FFTAcceleration_blocks_optFFT(order);
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_FFTW)
					return new FFTAcceleration_2wayM2L_fftw(order);
				else
					return new FFTAcceleration_2wayM2L_optFFT(order);
			} else {
				if (FFTSettings::USE_FFTW)
					return new FFTAcceleration_matrices_fftw(order);
				else
					return new FFTAcceleration_matrices_optFFT(order);
			}
		}
	}

	static TransferFunctionManagerAPI* getTransferFunctionManagerAPI(int ord,
			FFTAccelerationAPI* FFTA) {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID)
			return new TransferFunctionManager_UniformGrid(ord, FFTA,
					FFTSettings::TFMANAGER_VERBOSE);
		else
			return new TransferFunctionManager(ord, FFTA,
					FFTSettings::TFMANAGER_VERBOSE);
	}

};

#endif
