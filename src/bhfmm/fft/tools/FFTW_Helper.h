/*
 * FFTW_Helper.h
 *
 *  Created on: Sep 14, 2015
 *      Author: gallardjm
 */

#ifndef FFTW_HELPER_H_
#define FFTW_HELPER_H_

#include "bhfmm/fft/FFTSettings_preprocessor.h"
#include "bhfmm/fft/FFTAccelerableExpansion.h"
#include <fftw3.h>
#include <iostream>

/**
 * Used by fftw implementations of FFTAcceleration
 */
class FFTW_Helper {
public:
	FFTW_Helper(const int p, const int nx, const int ny);
	~FFTW_Helper();

	void FFT2Local(FFTAccelerableExpansion & Expansion, double radius);

	inline void execute_S2FFT() {
#if defined(__SINGLE_PRECISION_FFT__)
		fftwf_execute(S2FFT);
#else
		fftw_execute(S2FFT);
#endif
	}

	inline void execute_T2FFT() {
#if defined(__SINGLE_PRECISION_FFT__)
		fftwf_execute(T2FFT);
#else
		fftw_execute(T2FFT);
#endif
	}

	inline void execute_FFT2L() {
#if defined(__SINGLE_PRECISION_FFT__)
		fftwf_execute(FFT2L);
#else
		fftw_execute(FFT2L);
#endif
	}

#if defined(__SINGLE_PRECISION_FFT__)

	fftwf_complex* Source2FFT(FFTAccelerableExpansion & Expansion, double radius);
	fftwf_complex* TransferFunction2FFT(FFTAccelerableExpansion & Expansion);
	void getInLocal(fftwf_complex* & in) {in = FFT2L_in;}

private:
	int _p;
	int _nx;
	int _ny;
	fftwf_plan S2FFT;
	fftwf_complex *S2FFT_in;
	fftwf_complex *S2FFT_out;
	fftwf_plan T2FFT;
	fftwf_complex *T2FFT_in;
	fftwf_complex *T2FFT_out;
	fftwf_plan FFT2L;
	fftwf_complex *FFT2L_in;
	fftwf_complex *FFT2L_out;

#else

	fftw_complex* Source2FFT(FFTAccelerableExpansion & Expansion,
			double radius);
	fftw_complex* TransferFunction2FFT(FFTAccelerableExpansion & Expansion);
	void getInLocal(fftw_complex* & in) {
		in = FFT2L_in;
	}

private:
	int _p;
	int _nx;
	int _ny;
	fftw_plan S2FFT;
	fftw_complex *S2FFT_in;
	fftw_complex *S2FFT_out;
	fftw_plan T2FFT;
	fftw_complex *T2FFT_in;
	fftw_complex *T2FFT_out;
	fftw_plan FFT2L;
	fftw_complex *FFT2L_in;
	fftw_complex *FFT2L_out;

#endif

};

#endif
