/*
 * FFTAcceleration_blocks_fftw.cpp
 *
 *  Created on: Feb 12, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_blocks_fftw.h"

FFTAcceleration_blocks_fftw::FFTAcceleration_blocks_fftw(int order) {
	_p = order + 1;
	_nbLinePerBlock = 4;
	if (_p % _nbLinePerBlock != 0) {
		throw std::invalid_argument(
				"FFTAcceleration: _p % _nbLinePerBlock != 0");
	}

	_nbBlocks = _p / _nbLinePerBlock;
	_fft_nx = 2 * _nbLinePerBlock - 1; // litterature uses 2*_nbLinePerBlock (no -1)
	_fft_ny = 2 * _p - 1;

	_fftw_api = new FFTW_API(_fft_nx, _fft_ny);
}

void FFTAcceleration_blocks_fftw::FFT_initialize_Source(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_blocks* FFTData = getFFTData(Expansion);
	FFT_precision*** & Re = FFTData->Re;
	FFT_precision*** & Im = FFTData->Im;

	int b, i, j, trueI;
	FFT_precision minus_one_power_i = 1.0;
	FFT_precision minus_one_neg_j = 1.0;
	FFT_precision r = (FFT_precision) radius;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* in = _fftw_api->getIn_Forward();
	fftwf_complex* out;
#else
	fftw_complex* in = _fftw_api->getIn_Forward();
	fftw_complex* out;
#endif

	for (b = 0; b < _nbBlocks; b++) {
		minus_one_power_i = 1.0;

		for (i = 0; i < _nbLinePerBlock; i++) {
			trueI = i + _nbLinePerBlock * b;
			for (j = 0; j < _fft_ny; j++) {
				if (_p - 1 <= trueI + j && _p - 1 >= j - trueI) {
					if (j - _p + 1 < 0) {
						minus_one_neg_j = 1.0
								- 2.0 * (FFT_precision) ((_p - 1 - j) % 2);
						in[i * _fft_ny + j][0] =
								(FFT_precision) Expansion.get_C(trueI,
										_p - 1 - j) * minus_one_neg_j;
						in[i * _fft_ny + j][1] =
								(FFT_precision) Expansion.get_S(trueI,
										_p - 1 - j) * -minus_one_neg_j;
					} else {
						in[i * _fft_ny + j][0] =
								(FFT_precision) Expansion.get_C(trueI,
										j - _p + 1);
						in[i * _fft_ny + j][1] =
								(FFT_precision) Expansion.get_S(trueI,
										j - _p + 1);
					}

					in[i * _fft_ny + j][0] = minus_one_power_i
							* in[i * _fft_ny + j][0] / r;
					in[i * _fft_ny + j][1] = -minus_one_power_i
							* in[i * _fft_ny + j][1] / r; //we want to use use conj(Source)

				} else {
					in[i * _fft_ny + j][0] = 0.0;
					in[i * _fft_ny + j][1] = 0.0;
				}
			}

			minus_one_power_i *= -1.0;
			r *= (FFT_precision) radius;

		}

		for (; i < _fft_nx; i++)
			for (j = 0; j < _fft_ny; j++) {
				in[i * _fft_ny + j][0] = 0.0;
				in[i * _fft_ny + j][1] = 0.0;
			}

		out = _fftw_api->FFTAndGetOutput_Forward();

		for (i = 0; i < _fft_nx; i++)
			for (j = 0; j < _fft_ny; j++) {
				Re[b][i][j] = out[i * _fft_ny + j][0];
				Im[b][i][j] = out[i * _fft_ny + j][1];
			}
	}

}

void FFTAcceleration_blocks_fftw::FFT_initialize_TransferFunction(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_blocks* FFTData = getFFTData(Expansion);
	FFT_precision*** & Re = FFTData->Re;
	FFT_precision*** & Im = FFTData->Im;

	int b, i, j, trueI;
	FFT_precision minus_one_neg_j = 1.0;

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* in = _fftw_api->getIn_Forward();
	fftwf_complex* out;
#else
	fftw_complex* in = _fftw_api->getIn_Forward();
	fftw_complex* out;
#endif

	for (b = 0; b < _nbBlocks; b++) {

		for (i = 0; i < _nbLinePerBlock; i++) {
			trueI = i + _nbLinePerBlock * b;
			for (j = 0; j < _fft_ny; j++) {
				if (_p - 1 <= j + trueI && _p - 1 >= j - trueI) {
					if (_p - 1 - j < 0) {
						minus_one_neg_j = 1.0
								- 2.0 * (FFT_precision) ((j + 1 - _p) % 2);
						in[(_nbLinePerBlock - i - 1) * _fft_ny + j][0] =
								(FFT_precision) Expansion.get_C(trueI,
										j + 1 - _p) * minus_one_neg_j;
						in[(_nbLinePerBlock - i - 1) * _fft_ny + j][1] =
								(FFT_precision) Expansion.get_S(trueI,
										j + 1 - _p) * -minus_one_neg_j;
					} else {
						in[(_nbLinePerBlock - i - 1) * _fft_ny + j][0] =
								(FFT_precision) Expansion.get_C(trueI,
										_p - 1 - j);
						in[(_nbLinePerBlock - i - 1) * _fft_ny + j][1] =
								(FFT_precision) Expansion.get_S(trueI,
										_p - 1 - j);
					}

				} else {
					in[(_nbLinePerBlock - i - 1) * _fft_ny + j][0] = 0.0;
					in[(_nbLinePerBlock - i - 1) * _fft_ny + j][1] = 0.0;
				}
			}
		}

		for (; i < _fft_nx; i++)
			for (j = 0; j < _fft_ny; j++) {
				in[i * _fft_ny + j][0] = 0.0;
				in[i * _fft_ny + j][1] = 0.0;
			}

		out = _fftw_api->FFTAndGetOutput_Forward();

		for (i = 0; i < _fft_nx; i++)
			for (j = 0; j < _fft_ny; j++) {
				Re[b][i][j] = out[i * _fft_ny + j][0];
				Im[b][i][j] = out[i * _fft_ny + j][1];
			}
	}

}

void FFTAcceleration_blocks_fftw::FFT_finalize_Target(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_blocks* FFTData = getFFTData(Expansion);
	FFT_precision*** & Re = FFTData->Re;
	FFT_precision*** & Im = FFTData->Im;

	int b, i, j, trueI;
	FFT_precision scaling = 1.0 / (FFT_precision) (_fft_nx * _fft_ny);
	FFT_precision scaling2 = 1.0 / (FFT_precision) (_fft_nx * _fft_ny);

#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* in = _fftw_api->getIn_Backward();
	fftwf_complex* out;
#else
	fftw_complex* in = _fftw_api->getIn_Backward();
	fftw_complex* out;
#endif

	for (b = 0; b < _nbBlocks; b++) {
		for (i = 0; i < _fft_nx; i++)
			for (j = 0; j < _fft_ny; j++) {
				in[i * _fft_ny + j][0] = Re[b][i][j];
				in[i * _fft_ny + j][1] = Im[b][i][j];
			}

		out = _fftw_api->FFTAndGetOutput_Backward();

		for (i = 0; i < _nbLinePerBlock; i++) {
			trueI = i + _nbLinePerBlock * b;
			for (j = 0; j <= trueI; j++) {
				Expansion.get_C(trueI, j) += (double) (out[(_nbLinePerBlock - i
						- 1) * _fft_ny + (_fft_ny - j - 1)][0] * scaling);
				Expansion.get_S(trueI, j) += (double) (out[(_nbLinePerBlock - i
						- 1) * _fft_ny + (_fft_ny - j - 1)][1] * scaling);
			}
			scaling /= (FFT_precision) radius;
		}

		if (b != 0)
			for (i = 1; i < _nbLinePerBlock; i++) {
				trueI = i + _nbLinePerBlock * (b - 1);
				for (j = 0; j <= trueI; j++) {
					Expansion.get_C(trueI, j) += (double) (out[(_nbLinePerBlock
							* 2 - i - 1) * _fft_ny + (_fft_ny - j - 1)][0]
							* scaling2);
					Expansion.get_S(trueI, j) += (double) (out[(_nbLinePerBlock
							* 2 - i - 1) * _fft_ny + (_fft_ny - j - 1)][1]
							* scaling2);
				}
				scaling2 /= (FFT_precision) radius;
			}
		scaling2 /= (FFT_precision) radius;

	}
}

