/*
 * FFTAcceleration_matrices.cpp
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_matrices.h"

FFTDataContainer_matrices* FFTAcceleration_matrices::getFFTData(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_matrices* FFTData = NULL;
	if (!Expansion.issetFFTData()) {
		FFTData = new FFTDataContainer_matrices(_fft_nx, _fft_ny);
		FFTData->Re = alloc_matrix(_fft_nx, _fft_ny);
		FFTData->Im = alloc_matrix(_fft_nx, _fft_ny);
		Expansion._FFTData = FFTData;
	} else {
		FFTData = static_cast<FFTDataContainer_matrices*>(Expansion._FFTData);
	}

	return FFTData;
}

void FFTAcceleration_matrices::FFT_initialize_Target(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_matrices* FFTData = getFFTData(Expansion);

	clear_matrix(FFTData->Re, _fft_nx, _fft_ny);
	clear_matrix(FFTData->Im, _fft_nx, _fft_ny);
}

void FFTAcceleration_matrices::FFT_M2L(FFTAccelerableExpansion & Source,
		FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction) {
	FFTDataContainer_matrices* S_FFTData = getFFTData(Source);
	FFT_precision** & S_Re = S_FFTData->Re;
	FFT_precision** & S_Im = S_FFTData->Im;

	FFTDataContainer_matrices* T_FFTData = getFFTData(Target);
	FFT_precision** & T_Re = T_FFTData->Re;
	FFT_precision** & T_Im = T_FFTData->Im;

	FFTDataContainer_matrices* Tf_FFTData =
			static_cast<FFTDataContainer_matrices*>(TransferFunction);
	FFT_precision** & Tf_Re = Tf_FFTData->Re;
	FFT_precision** & Tf_Im = Tf_FFTData->Im;

	int i;
	const int end_i = _fft_nx * _fft_ny;
	const FFT_precision * const s_re = &(S_Re[0][0]);
	const FFT_precision * const s_im = &(S_Im[0][0]);
	const FFT_precision * const tr_re = &(Tf_Re[0][0]);
	const FFT_precision * const tr_im = &(Tf_Im[0][0]);
	FFT_precision * const t_re = &(T_Re[0][0]);
	FFT_precision * const t_im = &(T_Im[0][0]);

	for (i = 0; i < end_i; ++i) {
		t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
		t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
	}

}

void FFTAcceleration_matrices::FFT_M2L_vec(FFTAccelerableExpansion & Source,
		FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction) {
	FFTDataContainer_matrices* S_FFTData = getFFTData(Source);
	FFT_precision** & S_Re = S_FFTData->Re;
	FFT_precision** & S_Im = S_FFTData->Im;

	FFTDataContainer_matrices* T_FFTData = getFFTData(Target);
	FFT_precision** & T_Re = T_FFTData->Re;
	FFT_precision** & T_Im = T_FFTData->Im;

	FFTDataContainer_matrices* Tf_FFTData =
			static_cast<FFTDataContainer_matrices*>(TransferFunction);
	FFT_precision** & Tf_Re = Tf_FFTData->Re;
	FFT_precision** & Tf_Im = Tf_FFTData->Im;

	int i;
	const int end_i = _fft_nx * _fft_ny;
	const FFT_precision * const s_re = &(S_Re[0][0]);
	const FFT_precision * const s_im = &(S_Im[0][0]);
	const FFT_precision * const tr_re = &(Tf_Re[0][0]);
	const FFT_precision * const tr_im = &(Tf_Im[0][0]);
	FFT_precision * const t_re = &(T_Re[0][0]);
	FFT_precision * const t_im = &(T_Im[0][0]);

#pragma omp simd aligned (t_re, t_im, s_re, s_im: __FFT_MATRIX_ALIGNMENT__)
	for (i = 0; i < end_i; ++i) {
		t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
		t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
	}
}

