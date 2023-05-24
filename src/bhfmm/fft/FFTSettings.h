/*
 * FFTSettings.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTSETTINGS_H_
#define FFTSETTINGS_H_

#include <stdio.h>
#include <vector>
#include <string>
#include "bhfmm/fft/FFTSettings_preprocessor.h"

/**
 * Static class containing all the FFT Acceleration's settings
 * 
 * All settings are public and can be changed dynamically, the FFTFactory
 * take them into account to provide the right FFTAcceleration and the FMM
 * algorithm should also take them into account
 * Settings shouldn't be changed once a FFTAcceleration has been created
 */
class FFTSettings {
public:

	//! Seetings, default values in FFTSettings.cpp

	//FFT or FFTW (if none use the no FFT acceleration)
	static bool USE_FFT;  // set to use any FFT acceleration (default optFFT)
	static bool USE_FFTW; //set to use fftw instead of optFFT

	//TransferFunctionManager settings
	static bool USE_TFMANAGER_UNIFORMGRID; //set to use memoized transfer function (require uniform grid)
	//WARNING: if USE_TFMANAGER_UNIFORMGRID is set the tf should not be freed after an M2L, else it should be freed to avoid memory leaks

	static bool TFMANAGER_VERBOSE; //set to print the TFManager stats at its destruct

	static bool USE_VECTORIZATION; //set to use vectorized FFT (/!\ can be faked, see FFTSettings_preprocessor.h)
	static bool USE_2WAY_M2L; //set to use a 2way M2L (best on non uniform grid, requires vectorization)

	static bool USE_BLOCK; //set to use a block decomposition
	static bool USE_ADVBLOCK; //set to use a scaling block decomposition (requires USE_BLOCK)

	static bool USE_ORDER_REDUCTION; //Order Reduction scheme, only with blocks

	//sets Settings to the best default value for a given order
	static void autoSetting(int order);

	static bool issetFFTAcceleration() {
		return USE_FFT;
	}

	//Used in the bhfmm code for command line interface
	static void setOptions(std::string option);
	static std::vector<std::string> getAvailableOptions();
	static void printCurrentOptions();

};

#endif
