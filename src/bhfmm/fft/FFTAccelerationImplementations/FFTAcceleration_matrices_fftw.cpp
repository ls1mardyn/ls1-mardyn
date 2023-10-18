/*
 * FFTAcceleration_matrices_fftw.cpp
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_matrices_fftw.h"

FFTAcceleration_matrices_fftw::FFTAcceleration_matrices_fftw(int order) {
	_p = order + 1;
	_fft_nx = 2 * _p;
	_fft_ny = 2 * _p;

	_fftw_h = new FFTW_Helper(_p, _fft_nx, _fft_ny);
}

void FFTAcceleration_matrices_fftw::FFT_initialize_Source(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_matrices* FFTData = getFFTData(Expansion);
	FFT_precision** & Re = FFTData->Re;
	FFT_precision** & Im = FFTData->Im;
	int n, m;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* out = _fftw_h->Source2FFT(Expansion, radius);
#else
	fftw_complex* out = _fftw_h->Source2FFT(Expansion, radius);
#endif
	for (n = 0; n < _fft_nx; n++)
		for (m = 0; m < _fft_ny; m++) {
			Re[n][m] = out[n * _fft_ny + m][0];
			Im[n][m] = out[n * _fft_ny + m][1];
		}

}

void FFTAcceleration_matrices_fftw::FFT_initialize_TransferFunction(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_matrices* FFTData = getFFTData(Expansion);
	FFT_precision** & Re = FFTData->Re;
	FFT_precision** & Im = FFTData->Im;

	int n, m;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* out = _fftw_h->TransferFunction2FFT(Expansion);
#else
	fftw_complex* out = _fftw_h->TransferFunction2FFT(Expansion);
#endif

	for (n = 0; n < _fft_nx; n++)
		for (m = 0; m < _fft_ny; m++) {
			Re[n][m] = out[n * _fft_ny + m][0];
			Im[n][m] = out[n * _fft_ny + m][1];
		}

}

void FFTAcceleration_matrices_fftw::FFT_finalize_Target(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_matrices* FFTData = getFFTData(Expansion);
	FFT_precision** & Re = FFTData->Re;
	FFT_precision** & Im = FFTData->Im;

	int n, m;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* in;
#else
	fftw_complex* in;
#endif
	_fftw_h->getInLocal(in);

	for (n = 0; n < _fft_nx; n++)
		for (m = 0; m < _fft_ny; m++) {
			in[n * _fft_ny + m][0] = Re[n][m];
			in[n * _fft_ny + m][1] = Im[n][m];
		}

	_fftw_h->FFT2Local(Expansion, radius);
}

