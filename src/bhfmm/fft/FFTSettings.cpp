/*
 * FFTSettings.cpp
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#include "FFTSettings.h"

//! default settings value
bool FFTSettings::USE_FFT = true; // set to use FFT acceleration
bool FFTSettings::USE_FFTW = false; //set to use fftw instead of optFFT
bool FFTSettings::USE_TFMANAGER_UNIFORMGRID = true; //set to use memoized transfer function (require uniform grid)
bool FFTSettings::TFMANAGER_VERBOSE = true;
bool FFTSettings::USE_VECTORIZATION = false; //set to use vectorized FFT
bool FFTSettings::USE_2WAY_M2L = false; //set to use a 2way M2L (best on non uniform grid, requires vectorization)  
bool FFTSettings::USE_BLOCK = false;
bool FFTSettings::USE_ADVBLOCK = false;
bool FFTSettings::USE_ORDER_REDUCTION = false;

/*
 * Default settings:
 * 
 bool FFTSettings::USE_FFT = true; // set to use FFT acceleration
 bool FFTSettings::USE_FFTW = false; //set to use fftw instead of optFFT
 bool FFTSettings::USE_TFMANAGER_UNIFORMGRID = true; //set to use memoized transfer function (require uniform grid)
 bool FFTSettings::TFMANAGER_VERBOSE = true;
 bool FFTSettings::USE_VECTORIZATION = false; //set to use vectorized FFT
 bool FFTSettings::USE_2WAY_M2L = false; //set to use a 2way M2L (best on non uniform grid, requires vectorization)
 bool FFTSettings::USE_BLOCK = false;
 *
 */

void FFTSettings::autoSetting(int order) {
	if (order < 16) {
		FFTSettings::USE_FFT = true;
		FFTSettings::USE_FFTW = false;
		FFTSettings::USE_TFMANAGER_UNIFORMGRID = true;
		FFTSettings::USE_VECTORIZATION = true;
		FFTSettings::USE_2WAY_M2L = false;
		FFTSettings::USE_BLOCK = false;
		FFTSettings::USE_ADVBLOCK = false;
		FFTSettings::USE_ORDER_REDUCTION = false;
	} else {
		FFTSettings::USE_FFT = true;
		FFTSettings::USE_FFTW = false;
		FFTSettings::USE_TFMANAGER_UNIFORMGRID = true;
		FFTSettings::USE_VECTORIZATION = true;
		FFTSettings::USE_2WAY_M2L = false;
		FFTSettings::USE_BLOCK = true;
		FFTSettings::USE_ADVBLOCK = true;
		FFTSettings::USE_ORDER_REDUCTION = true;
	}
}

void FFTSettings::setOptions(std::string option) {
	if (option == "off") {
		FFTSettings::USE_FFT = false;
	} else if (option == "fft") {
		FFTSettings::USE_FFT = true;
		FFTSettings::USE_FFTW = false;
	} else if (option == "fftw") {
		FFTSettings::USE_FFT = true;
		FFTSettings::USE_FFTW = true;
	} else if (option == "TFMem") {
		FFTSettings::USE_TFMANAGER_UNIFORMGRID = true;
	} else if (option == "noTFMem") {
		FFTSettings::USE_TFMANAGER_UNIFORMGRID = false;
	} else if (option == "vec") {
		FFTSettings::USE_VECTORIZATION = true;
	} else if (option == "novec") {
		FFTSettings::USE_VECTORIZATION = false;
	} else if (option == "2way") {
		FFTSettings::USE_2WAY_M2L = true;
	} else if (option == "no2way") {
		FFTSettings::USE_2WAY_M2L = false;
	} else if (option == "block") {
		FFTSettings::USE_BLOCK = true;
		FFTSettings::USE_ADVBLOCK = false;
	} else if (option == "advblock") {
		FFTSettings::USE_BLOCK = true;
		FFTSettings::USE_ADVBLOCK = true;
	} else if (option == "noblock") {
		FFTSettings::USE_BLOCK = false;
		FFTSettings::USE_ADVBLOCK = false;
	} else if (option == "ORed") {
		FFTSettings::USE_ORDER_REDUCTION = true;
	} else if (option == "noORed") {
		FFTSettings::USE_ORDER_REDUCTION = false;
	} else {
		printf("Unrecognized FFT Option!!! \nUsing default parameters:\n");
		FFTSettings::printCurrentOptions();
	}
}

std::vector<std::string> FFTSettings::getAvailableOptions() {
	std::vector<std::string> options;
	options.push_back("off:     disable fft acceleration");
	options.push_back("fft:     use opt. fft acceleration");
	options.push_back("fftw:    use FFTW acceleration");
	options.push_back(
			"TFMem:   enable TransferFunction precomputation (requires uniform grid)");
	options.push_back("noTFMem: disable TransferFunction precomputation");
	options.push_back(
			"vec:     enable vectorization (must have been compiled with vec support, see FFTSettings_preprocessor.h)");
	options.push_back("novec:   disable vectorization");
	options.push_back(
			"2way:    enable 2way M2L (better on non uniform grid), requires vec");
	options.push_back("no2way:  disable 2way M2L");
	options.push_back("block:   enable block FFT (only with fftw)");
	options.push_back("advblock:enable advanced block FFT (only with fftw)");
	options.push_back("noblock: disable block FFT");
	options.push_back(
			"ORed:    enable order reduction (block and advblock only)");
	options.push_back("noORed:  disable order reduction");
	return options;
}

void FFTSettings::printCurrentOptions() {
	if (FFTSettings::USE_FFT) {
		if (FFTSettings::USE_FFTW)
			printf("FFTW Acceleration\n");
		else
			printf("opt. FFT Acceleration\n");

		if (FFTSettings::USE_BLOCK) {
			if (FFTSettings::USE_ADVBLOCK)
				printf("Adv. Block FFT enabled\n");
			else
				printf("Block FFT enabled\n");
		}

		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID)
			printf(
					"TransferFunctionManager_UniformGrid (precomputation of TF, requires an uniform grid)\n");
		else
			printf("default TransferFunctionManager (no precomputation)\n");

		if (FFTSettings::USE_VECTORIZATION)
#if defined (__TEST_FAKE_VECTORIZATION__)
			printf("/!\\ FAKED Vectorization. Alignement: %i bits\n", __FFT_MATRIX_ALIGNMENT__);
#else
			printf("Vectorized M2L enabled. Alignement: %i bits\n",
					__FFT_MATRIX_ALIGNMENT__);
#endif

		if (FFTSettings::USE_2WAY_M2L)
			printf("2way M2L enabled\n");

		if (FFTSettings::USE_ORDER_REDUCTION)
			printf("Order Reduction enabled\n");

	} else
		printf("No FFT Acceleration\n");

}
