/*
 * FFTW_API.h
 *
 *  Created on: Sep 14, 2015
 *      Author: gallardjm
 */

#ifndef FFTW_API_H_
#define FFTW_API_H_

#include <stddef.h> //NULL
#include <fftw3.h>
#include "bhfmm/fft/FFTSettings_preprocessor.h"

/**
 * Basic API to perform 2d FFT and IFFT using the fftw lib
 *
 * Can perform 1d FFT
 * Uses double or single precision depending on FFTSettings_preprocessor
 *
 * getIn_Forward() to get the input matrices to fill with the value
 * FFTAndGetOutput_Forward() to perform the FFT and get the output
 * Backward for IFFT
 */
class FFTW_API {

public:

	/**
	 * Constructor, give ny=1 to perform 1d FFT
	 *
	 * @param int nx
	 * @param int ny
	 */
	FFTW_API(int nx, int ny) :
			_nx(nx), _ny(ny) {
#if defined(__SINGLE_PRECISION_FFT__)
		_f_in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
		_f_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
		_b_in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
		_b_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);

		_forward = fftwf_plan_dft_2d(nx, ny, _f_in, _f_out, FFTW_FORWARD, FFTW_ESTIMATE);
		_backward = fftwf_plan_dft_2d(nx, ny, _b_in, _b_out, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
		_f_in = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * nx * ny);
		_f_out = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * nx * ny);
		_b_in = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * nx * ny);
		_b_out = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * nx * ny);

		//if(ny>1) {
		_forward = fftw_plan_dft_2d(nx, ny, _f_in, _f_out, FFTW_FORWARD,
				FFTW_ESTIMATE);
		_backward = fftw_plan_dft_2d(nx, ny, _b_in, _b_out, FFTW_BACKWARD,
				FFTW_ESTIMATE);
		//} else {
		//  _forward  = fftw_plan_dft_1d(nx, _f_in, _f_out, FFTW_FORWARD,  FFTW_ESTIMATE);
		//  _backward = fftw_plan_dft_1d(nx, _b_in, _b_out, FFTW_BACKWARD, FFTW_ESTIMATE);
		//}
#endif
	}

	//! destructor, clean plans and frees allocated memory
	~FFTW_API() {
#if defined(__SINGLE_PRECISION_FFT__)
		fftwf_destroy_plan(_forward);
		fftwf_destroy_plan(_backward);
		fftwf_free(_f_in);
		fftwf_free(_f_out);
		fftwf_free(_b_in);
		fftwf_free(_b_out);
#else
		fftw_destroy_plan(_forward);
		fftw_destroy_plan(_backward);
		fftw_free(_f_in);
		fftw_free(_f_out);
		fftw_free(_b_in);
		fftw_free(_b_out);

		fftw_cleanup();
#endif
	}

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* getIn_Forward() {return _f_in;}
	fftwf_complex* getIn_Backward() {return _b_in;}

	fftwf_complex* FFTAndGetOutput_Forward() {fftwf_execute(_forward); return _f_out;}
	fftwf_complex* FFTAndGetOutput_Backward() {fftwf_execute(_backward); return _b_out;}
#else
	fftw_complex* getIn_Forward() {
		return _f_in;
	}
	fftw_complex* getIn_Backward() {
		return _b_in;
	}

	fftw_complex* FFTAndGetOutput_Forward() {
		fftw_execute(_forward);
		return _f_out;
	}
	fftw_complex* FFTAndGetOutput_Backward() {
		fftw_execute(_backward);
		return _b_out;
	}
#endif

private:
	int _nx;
	int _ny;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_plan _forward;
	fftwf_complex* _f_in;
	fftwf_complex* _f_out;
	fftwf_plan _backward;
	fftwf_complex* _b_in;
	fftwf_complex* _b_out;
#else
	fftw_plan _forward;
	fftw_complex* _f_in;
	fftw_complex* _f_out;
	fftw_plan _backward;
	fftw_complex* _b_in;
	fftw_complex* _b_out;
#endif
};

#endif
