/*
 * FFTAcceleration_2wayM2L_fftw.cpp
 *
 *  Created on: Feb 10, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_2wayM2L_fftw.h"

FFTAcceleration_2wayM2L_fftw::FFTAcceleration_2wayM2L_fftw(int order) {
	_p = order + 1;
	_fft_nx = 2 * _p;
	_fft_ny = 2 * _p - 1;

	_nbZeroes = nb_zeroes_2wayFFT(_fft_nx, _fft_ny);
	_totalSize = _fft_nx * _fft_ny + 2 * _nbZeroes;

	_fftw_h = new FFTW_Helper(_p, _fft_nx, _fft_ny);
}

void FFTAcceleration_2wayM2L_fftw::FFT_initialize_Source(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_arrays* FFTData = getFFTData(Expansion);
	FFT_precision* Re = FFTData->Re;
	FFT_precision* Im = FFTData->Im;

	int n;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* out = _fftw_h->Source2FFT(Expansion, radius);
#else
	fftw_complex* out = _fftw_h->Source2FFT(Expansion, radius);
#endif

	for (n = 0; n < _fft_nx * _fft_ny / 2; n++) {
		Re[n] = out[n][0];
		Im[n] = out[n][1];
	}
	for (; n < _fft_nx * _fft_ny; n++) {
		Re[n + _nbZeroes] = out[n][0];
		Im[n + _nbZeroes] = out[n][1];
	}
}

void FFTAcceleration_2wayM2L_fftw::FFT_initialize_TransferFunction(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_arrays* FFTData = getFFTData(Expansion);
	FFT_precision* Re = FFTData->Re;
	FFT_precision* Im = FFTData->Im;

	int n;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* out = _fftw_h->TransferFunction2FFT(Expansion);
#else
	fftw_complex* out = _fftw_h->TransferFunction2FFT(Expansion);
#endif

	for (n = 0; n < _fft_nx * _fft_ny / 2; n++) {
		Re[n] = out[n][0];
		Im[n] = out[n][1];
	}
	for (; n < _fft_nx * _fft_ny; n++) {
		Re[n + _nbZeroes] = out[n][0];
		Im[n + _nbZeroes] = out[n][1];
	}

}

void FFTAcceleration_2wayM2L_fftw::FFT_finalize_Target(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_arrays* FFTData = getFFTData(Expansion);
	FFT_precision* Re = FFTData->Re;
	FFT_precision* Im = FFTData->Im;

	int n;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* in;
#else
	fftw_complex* in;
#endif
	_fftw_h->getInLocal(in);

	for (n = 0; n < _fft_nx * _fft_ny / 2; n++) {
		in[n][0] = Re[n];
		in[n][1] = Im[n];
	}
	for (; n < _fft_nx * _fft_ny; n++) {
		in[n][0] = Re[n + _nbZeroes];
		in[n][1] = Im[n + _nbZeroes];
	}

	_fftw_h->FFT2Local(Expansion, radius);
}

