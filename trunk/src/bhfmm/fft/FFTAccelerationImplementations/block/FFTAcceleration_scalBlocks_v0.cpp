/*
 * FFTAcceleration_scalBlocks_v0.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_scalBlocks_v0.h"

FFTAcceleration_scalBlocks_v0::FFTAcceleration_scalBlocks_v0(int order) {
	printf("scal Block\n");
	_p = order + 1;
	_nbLinePerBlock = 4;
	if (_p % _nbLinePerBlock != 0) {
		throw std::invalid_argument(
				"FFTAcceleration: _p % _nbLinePerBlock != 0");
	}

	_nbBlocks = _p / _nbLinePerBlock;
	_fft_nx = 2 * _nbLinePerBlock - 1; // litterature uses 2*_nbLinePerBlock (no -1)
	_fft_ny = _p;

	_optFFT_API = optFFT_API_Factory::getOptFFT_API(order, false);
	_blockSize = new int[_nbBlocks];

	//initializing _blockSize
	int repetition, i;

	for (i = 0; i < _nbBlocks / 2; i++) { //first half should have smaller blocks
		repetition = _fft_ny / (_nbLinePerBlock * (i + 1)); //will be at least 2
		while (_fft_ny % repetition != 0)
			repetition--;
		_blockSize[i] = _fft_ny / repetition;
	}
	for (; i < _nbBlocks; i++)
		_blockSize[i] = _fft_ny;
}

FFTDataContainer_scalBlocks_v0* FFTAcceleration_scalBlocks_v0::getFFTData(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_scalBlocks_v0* FFTData = NULL;
	if (!Expansion.issetFFTData()) {
		FFTData = new FFTDataContainer_scalBlocks_v0(_nbBlocks, _fft_nx,
				_fft_ny, false, _blockSize);
		FFTData->Re = alloc_blocks(_nbBlocks, _fft_nx, _fft_ny);
		FFTData->Im = alloc_blocks(_nbBlocks, _fft_nx, _fft_ny);
		Expansion._FFTData = FFTData;
	} else {
		FFTData =
				static_cast<FFTDataContainer_scalBlocks_v0*>(Expansion._FFTData);
	}

	return FFTData;
}

FFTDataContainer_scalBlocks_v0* FFTAcceleration_scalBlocks_v0::getFFTData_scal(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_scalBlocks_v0* FFTData = NULL;
	if (!Expansion.issetFFTData()) {
		FFTData = new FFTDataContainer_scalBlocks_v0(_nbBlocks, _fft_nx,
				_fft_ny, true, _blockSize);
		FFTData->Re = alloc_scalBlocks(_nbBlocks, _fft_nx, _blockSize);
		FFTData->Im = alloc_scalBlocks(_nbBlocks, _fft_nx, _blockSize);
		Expansion._FFTData = FFTData;
	} else {
		FFTData =
				static_cast<FFTDataContainer_scalBlocks_v0*>(Expansion._FFTData);
	}

	return FFTData;
}

void FFTAcceleration_scalBlocks_v0::FFT_initialize_Target(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_scalBlocks_v0* FFTData = getFFTData(Expansion);

	clear_blocks(FFTData->Re, _nbBlocks, _fft_nx, _fft_ny);
	clear_blocks(FFTData->Im, _nbBlocks, _fft_nx, _fft_ny);
}

void FFTAcceleration_scalBlocks_v0::FFT_initialize_Source(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_scalBlocks_v0* FFTData = getFFTData(Expansion);
	FFT_precision*** & Re = FFTData->Re;
	FFT_precision*** & Im = FFTData->Im;

	int b, i, j, trueI;
	FFT_precision scaling = 1.0 / (FFT_precision) radius; //include radius scaling and -1^i

	for (b = 0; b < _nbBlocks; b++) {
		for (i = 0; i < _nbLinePerBlock; i++) {
			trueI = i + _nbLinePerBlock * b;
			for (j = 0; j <= trueI; j++) {
				Re[b][i][j] = (FFT_precision) Expansion.get_C(trueI, j)
						* scaling;
				Im[b][i][j] = -(FFT_precision) Expansion.get_S(trueI, j)
						* scaling; //we want to use use conj(Source)
			}
			for (; j < _fft_ny; j++) {
				Re[b][i][j] = 0.0;
				Im[b][i][j] = 0.0;
			}
			scaling /= -(FFT_precision) radius;
		}

		for (; i < _fft_nx; i++)
			for (j = 0; j < _fft_ny; j++) {
				Re[b][i][j] = 0.0;
				Im[b][i][j] = 0.0;
			}

		_optFFT_API->optimizedFFT(Re[b], Im[b], _fft_nx, _fft_ny);
	}
}

void FFTAcceleration_scalBlocks_v0::FFT_finalize_Target(
		FFTAccelerableExpansion & Expansion, double radius) {
	FFTDataContainer_scalBlocks_v0* FFTData = getFFTData(Expansion);
	FFT_precision*** & Re = FFTData->Re;
	FFT_precision*** & Im = FFTData->Im;

	int b, i, j, trueI;
	double scaling = 1.0; // / (FFT_precision)(_fft_nx * _fft_ny);
	double scaling2 = 1.0; // / (FFT_precision)(_fft_nx * _fft_ny);
	double minus_one_power_j;

	for (b = 0; b < _nbBlocks; b++) {

		_optFFT_API->optimizedIFFT(Re[b], Im[b], _fft_nx, _fft_ny);

		for (i = 0; i < _nbLinePerBlock; i++) {
			trueI = i + _nbLinePerBlock * b;
			minus_one_power_j = 1.0;
			for (j = 0; j <= trueI; j++) {
				Expansion.get_C(trueI, j) += minus_one_power_j
						* (double) (Re[b][_nbLinePerBlock - i - 1][j])
						* scaling;
				Expansion.get_S(trueI, j) += -minus_one_power_j
						* (double) (Im[b][_nbLinePerBlock - i - 1][j])
						* scaling;
				minus_one_power_j *= -1.0;
			}
			scaling /= radius;
		}

		if (b != 0) {
			for (i = 1; i < _nbLinePerBlock; i++) {
				trueI = i + _nbLinePerBlock * (b - 1);
				minus_one_power_j = 1.0;
				for (j = 0; j <= trueI; j++) {
					Expansion.get_C(trueI, j) += minus_one_power_j
							* (double) (Re[b][_nbLinePerBlock * 2 - i - 1][j])
							* scaling2;
					Expansion.get_S(trueI, j) += -minus_one_power_j
							* (double) (Im[b][_nbLinePerBlock * 2 - i - 1][j])
							* scaling2;
					minus_one_power_j *= -1.0;
				}
				scaling2 /= radius;
			}
		}
		scaling2 /= (FFT_precision) radius;
	}
}

void FFTAcceleration_scalBlocks_v0::FFT_M2L(FFTAccelerableExpansion & Source,
		FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction) {
	FFTDataContainer_scalBlocks_v0* S_FFTData = getFFTData(Source);
	FFT_precision*** & S_Re = S_FFTData->Re;
	FFT_precision*** & S_Im = S_FFTData->Im;

	FFTDataContainer_scalBlocks_v0* T_FFTData = getFFTData(Target);
	FFT_precision*** & T_Re = T_FFTData->Re;
	FFT_precision*** & T_Im = T_FFTData->Im;

	FFTDataContainer_scalBlocks_v0* Tf_FFTData =
			static_cast<FFTDataContainer_scalBlocks_v0*>(TransferFunction);
	FFT_precision*** & Tf_Re = Tf_FFTData->Re;
	FFT_precision*** & Tf_Im = Tf_FFTData->Im;

	int i, b;
	int currentBlock, tfBlock;
	int lastBlock = _nbBlocks - 1; // could be reduced for far way interaction
	int repetition;
	int j, j_rep;
	const int end_i = _fft_nx * _fft_ny; //size of a block
	FFT_precision * t_re;
	FFT_precision * t_im;
	FFT_precision * s_re;
	FFT_precision * s_im;
	FFT_precision * tr_re;
	FFT_precision * tr_im;

	for (currentBlock = 0; currentBlock <= lastBlock; currentBlock++) {
		t_re = &(T_Re[currentBlock][0][0]);
		t_im = &(T_Im[currentBlock][0][0]);
		for (b = 0; b <= lastBlock - currentBlock; b++) {
			tfBlock = b + currentBlock;
			if (tfBlock < _nbBlocks / 2) { //scaling blocks scheme
				repetition = _fft_ny / _blockSize[tfBlock];
				for (j = 0; j < _blockSize[tfBlock]; j++) {
					j_rep = repetition * j;
					for (i = 0; i < _fft_nx; i++) {
						T_Re[currentBlock][i][j_rep] += S_Re[b][i][j_rep]
								* Tf_Re[tfBlock][i][j]
								- S_Im[b][i][j_rep] * Tf_Im[tfBlock][i][j];
						T_Im[currentBlock][i][j_rep] += S_Re[b][i][j_rep]
								* Tf_Im[tfBlock][i][j]
								+ S_Im[b][i][j_rep] * Tf_Re[tfBlock][i][j];
					}
				}

			} else { //same size block, default scheme
				s_re = &(S_Re[b][0][0]);
				s_im = &(S_Im[b][0][0]);
				tr_re = &(Tf_Re[tfBlock][0][0]);
				tr_im = &(Tf_Im[tfBlock][0][0]);
				for (i = 0; i < end_i; ++i) {
					t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
					t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
				}
			}
		}
	}

}

void FFTAcceleration_scalBlocks_v0::FFT_M2L_vec(
		FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target,
		FFTDataContainer* TransferFunction) {
	return FFT_M2L(Source, Target, TransferFunction);
} //TODO

void FFTAcceleration_scalBlocks_v0::FFT_initialize_TransferFunction(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_scalBlocks_v0* FFTData = getFFTData_scal(Expansion);
	FFT_precision*** & Re = FFTData->Re;
	FFT_precision*** & Im = FFTData->Im;

	int b, i, j, trueI;
	FFT_precision minus_one_power_j = 1.0;
	FFT_precision repetition;

	for (b = 0; b < _nbBlocks; b++) {
		for (i = 0; i < _nbLinePerBlock; i++) {
			trueI = i + _nbLinePerBlock * b;
			minus_one_power_j = 1.0;
			for (j = 0; j <= trueI; j++) {
				Re[b][_nbLinePerBlock - i - 1][j] =
						(FFT_precision) Expansion.get_C(trueI, j)
								* minus_one_power_j;
				Im[b][_nbLinePerBlock - i - 1][j] =
						(FFT_precision) Expansion.get_S(trueI, j)
								* -minus_one_power_j;
				minus_one_power_j *= -1.0;
			}

			for (; j < _blockSize[b]; j++) {
				Re[b][_nbLinePerBlock - i - 1][j] = 0.0;
				Im[b][_nbLinePerBlock - i - 1][j] = 0.0;
			}
		}

		for (; i < _fft_nx; i++)
			for (j = 0; j < _blockSize[b]; j++) {
				Re[b][i][j] = 0.0;
				Im[b][i][j] = 0.0;
			}

		_optFFT_API->optimizedFFT(Re[b], Im[b], _fft_nx, _blockSize[b]);

		if (b < _nbBlocks / 2) {
			repetition = (FFT_precision) (_fft_ny / _blockSize[b]);
			for (i = 0; i < _fft_nx; i++)
				for (j = 0; j < _blockSize[b]; j++) {
					Re[b][i][j] *= repetition;
					Im[b][i][j] *= repetition;
				}
		}
	}
}

