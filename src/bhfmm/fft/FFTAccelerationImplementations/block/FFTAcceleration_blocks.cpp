/*
 * FFTAcceleration_blocks.cpp
 *
 *  Created on: Feb 12, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_blocks.h"

FFTDataContainer_blocks* FFTAcceleration_blocks::getFFTData(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_blocks* FFTData = NULL;
	if (!Expansion.issetFFTData()) {
		FFTData = new FFTDataContainer_blocks(_nbBlocks, _fft_nx, _fft_ny);
		FFTData->Re = alloc_blocks(_nbBlocks, _fft_nx, _fft_ny);
		FFTData->Im = alloc_blocks(_nbBlocks, _fft_nx, _fft_ny);
		Expansion._FFTData = FFTData;
	} else {
		FFTData = static_cast<FFTDataContainer_blocks*>(Expansion._FFTData);
	}

	return FFTData;
}

void FFTAcceleration_blocks::FFT_initialize_Target(
		FFTAccelerableExpansion & Expansion) {
	FFTDataContainer_blocks* FFTData = getFFTData(Expansion);

	clear_blocks(FFTData->Re, _nbBlocks, _fft_nx, _fft_ny);
	clear_blocks(FFTData->Im, _nbBlocks, _fft_nx, _fft_ny);
}

template<bool Vect, bool OrderRed>
void FFTAcceleration_blocks::FFT_M2L_template(FFTAccelerableExpansion & Source,
		FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction,
		int order) {
	FFTDataContainer_blocks* S_FFTData = getFFTData(Source);
	FFT_precision*** & S_Re = S_FFTData->Re;
	FFT_precision*** & S_Im = S_FFTData->Im;

	FFTDataContainer_blocks* T_FFTData = getFFTData(Target);
	FFT_precision*** & T_Re = T_FFTData->Re;
	FFT_precision*** & T_Im = T_FFTData->Im;

	FFTDataContainer_blocks* Tf_FFTData =
			static_cast<FFTDataContainer_blocks*>(TransferFunction);
	FFT_precision*** & Tf_Re = Tf_FFTData->Re;
	FFT_precision*** & Tf_Im = Tf_FFTData->Im;

	int i, b;
	int currentBlock;

	int lastBlock;

	if (OrderRed) {
		lastBlock = (order / _nbLinePerBlock); //compute up to order
		if (order % _nbLinePerBlock == 0)
			lastBlock--;
		if (lastBlock < 0)
			lastBlock = 0;
		else if (lastBlock > _nbBlocks - 1)
			lastBlock = _nbBlocks - 1;
	} else {
		lastBlock = _nbBlocks - 1;
	}

	const int end_i = _fft_nx * _fft_ny; //size of a block
	FFT_precision * t_re;
	FFT_precision * t_im;
	FFT_precision * s_re;
	FFT_precision * s_im;
	FFT_precision * tf_re;
	FFT_precision * tf_im;

	for (currentBlock = 0; currentBlock <= lastBlock; currentBlock++) {
		t_re = &(T_Re[currentBlock][0][0]);
		t_im = &(T_Im[currentBlock][0][0]);
		for (b = 0; b <= lastBlock - currentBlock; b++) {
			s_re = &(S_Re[b][0][0]);
			s_im = &(S_Im[b][0][0]);
			tf_re = &(Tf_Re[b + currentBlock][0][0]);
			tf_im = &(Tf_Im[b + currentBlock][0][0]);

			if (Vect) {
#pragma omp simd aligned (t_re, t_im, s_re, s_im, tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
				for (i = 0; i < end_i; ++i) {
					t_re[i] += s_re[i] * tf_re[i] - s_im[i] * tf_im[i];
					t_im[i] += s_re[i] * tf_im[i] + s_im[i] * tf_re[i];
				}
			} else {
				for (i = 0; i < end_i; ++i) {
					t_re[i] += s_re[i] * tf_re[i] - s_im[i] * tf_im[i];
					t_im[i] += s_re[i] * tf_im[i] + s_im[i] * tf_re[i];
				}
			}
		}
	}

}

void FFTAcceleration_blocks::FFT_M2L(FFTAccelerableExpansion & Source,
		FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction) {
	return FFT_M2L_template<false, false>(Source, Target, TransferFunction, 0);
}

void FFTAcceleration_blocks::FFT_M2L_vec(FFTAccelerableExpansion & Source,
		FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction) {
	return FFT_M2L_template<true, false>(Source, Target, TransferFunction, 0);
}

void FFTAcceleration_blocks::FFT_M2L_OrderReduction(
		FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target,
		FFTDataContainer* TransferFunction, int order) {
	return FFT_M2L_template<false, true>(Source, Target, TransferFunction,
			order);
}

void FFTAcceleration_blocks::FFT_M2L_OrderReduction_vec(
		FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target,
		FFTDataContainer* TransferFunction, int order) {
	return FFT_M2L_template<true, true>(Source, Target, TransferFunction, order);
}

template<bool Vect, bool OrderRed>
void FFTAcceleration_blocks::FFT_M2L_2way_template(
		FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2,
		FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2,
		FFTDataContainer* TransferFunction, int order) {
	FFTDataContainer_blocks* Source1Data = getFFTData(Source1);
	FFT_precision*** & S1_Re = Source1Data->Re;
	FFT_precision*** & S1_Im = Source1Data->Im;

	FFTDataContainer_blocks* Source2Data = getFFTData(Source2);
	FFT_precision*** & S2_Re = Source2Data->Re;
	FFT_precision*** & S2_Im = Source2Data->Im;

	FFTDataContainer_blocks* Target1Data = getFFTData(Target1);
	FFT_precision*** & T1_Re = Target1Data->Re;
	FFT_precision*** & T1_Im = Target1Data->Im;

	FFTDataContainer_blocks* Target2Data = getFFTData(Target2);
	FFT_precision*** & T2_Re = Target2Data->Re;
	FFT_precision*** & T2_Im = Target2Data->Im;

	FFTDataContainer_blocks* TFData =
			static_cast<FFTDataContainer_blocks*>(TransferFunction);
	FFT_precision*** & Tf_Re = TFData->Re;
	FFT_precision*** & Tf_Im = TFData->Im;

	FFT_precision * t1_re;
	FFT_precision * t1_im;
	FFT_precision * s1_re;
	FFT_precision * s1_im;
	FFT_precision * t2_re;
	FFT_precision * t2_im;
	FFT_precision * s2_re;
	FFT_precision * s2_im;
	FFT_precision * tf_re;
	FFT_precision * tf_im;

	int i, j, b;
	int currentBlock;

	int lastBlock;

	if (OrderRed) {
		lastBlock = (order / _nbLinePerBlock); //compute up to order
		if (order % _nbLinePerBlock == 0)
			lastBlock--;
		if (lastBlock < 0)
			lastBlock = 0;
		else if (lastBlock > _nbBlocks - 1)
			lastBlock = _nbBlocks - 1;
	} else {
		lastBlock = _nbBlocks - 1;
	}

	const int end_i = _fft_nx * _fft_ny; //size of a block
	const int shift_i = end_i / 2;
	/*
	 * Tr(r)[i][j] = -1^(p+1) Tr(-r)[i+max_i/2][j]
	 * see doc/2WayM2L.txt
	 */
	for (currentBlock = 0; currentBlock <= lastBlock; currentBlock++) {
		t1_re = &(T1_Re[currentBlock][0][0]);
		t1_im = &(T1_Im[currentBlock][0][0]);
		t2_re = &(T2_Re[currentBlock][0][0]);
		t2_im = &(T2_Im[currentBlock][0][0]);
		for (b = 0; b <= lastBlock - currentBlock; b++) {
			s1_re = &(S1_Re[b][0][0]);
			s1_im = &(S1_Im[b][0][0]);
			s2_re = &(S2_Re[b][0][0]);
			s2_im = &(S2_Im[b][0][0]);
			tf_re = &(Tf_Re[b + currentBlock][0][0]);
			tf_im = &(Tf_Im[b + currentBlock][0][0]);

			if (Vect) {
				j = shift_i;
#pragma omp simd aligned (t1_re, t1_im, s1_re, s1_im, \
								  t2_re, t2_im, s2_re, s2_im, \
								  tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
				for (i = 0; i < shift_i; ++i) {
					t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
					t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
					t1_im[j] -= s2_re[j] * tf_im[i] + s2_im[j] * tf_re[i];
					t1_re[j] -= s2_re[j] * tf_re[i] - s2_im[j] * tf_im[i];
					++j;
				}
				j = 0;
#pragma omp simd aligned (t1_re, t1_im, s1_re, s1_im, \
								  t2_re, t2_im, s2_re, s2_im, \
								  tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
				for (i = shift_i; i < end_i; ++i) {
					t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
					t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
					t1_im[j] -= s2_re[j] * tf_im[i] + s2_im[j] * tf_re[i];
					t1_re[j] -= s2_re[j] * tf_re[i] - s2_im[j] * tf_im[i];
					++j;
				}
			} else {
				j = shift_i;
				for (i = 0; i < shift_i; ++i) {
					t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
					t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
					t1_im[j] -= s2_re[j] * tf_im[i] + s2_im[j] * tf_re[i];
					t1_re[j] -= s2_re[j] * tf_re[i] - s2_im[j] * tf_im[i];
					++j;
				}
				j = 0;
				for (; i < end_i; ++i) {
					t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
					t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
					t1_im[j] -= s2_re[j] * tf_im[i] + s2_im[j] * tf_re[i];
					t1_re[j] -= s2_re[j] * tf_re[i] - s2_im[j] * tf_im[i];
					++j;
				}
			}
		}
	}
}

void FFTAcceleration_blocks::FFT_M2L_2way(FFTAccelerableExpansion & Source1,
		FFTAccelerableExpansion & Source2, FFTAccelerableExpansion & Target1,
		FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction) {
	return FFT_M2L_2way_template<false, false>(Source1, Source2, Target1,
			Target2, TransferFunction, 0);
}

void FFTAcceleration_blocks::FFT_M2L_2way_vec(FFTAccelerableExpansion & Source1,
		FFTAccelerableExpansion & Source2, FFTAccelerableExpansion & Target1,
		FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction) {
	return FFT_M2L_2way_template<true, false>(Source1, Source2, Target1,
			Target2, TransferFunction, 0);
}

void FFTAcceleration_blocks::FFT_M2L_2way_ORed(
		FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2,
		FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2,
		FFTDataContainer* TransferFunction, int order) {
	return FFT_M2L_2way_template<false, true>(Source1, Source2, Target1,
			Target2, TransferFunction, order);
}

void FFTAcceleration_blocks::FFT_M2L_2way_ORed_vec(
		FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2,
		FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2,
		FFTDataContainer* TransferFunction, int order) {
	return FFT_M2L_2way_template<true, true>(Source1, Source2, Target1, Target2,
			TransferFunction, order);
}

/*
 void FFTAcceleration_blocks::FFT_M2L(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction)
 {
 FFTDataContainer_blocks* S_FFTData = getFFTData(Source);
 FFT_precision*** & S_Re = S_FFTData->Re;
 FFT_precision*** & S_Im = S_FFTData->Im;

 FFTDataContainer_blocks* T_FFTData = getFFTData(Target);
 FFT_precision*** & T_Re = T_FFTData->Re;
 FFT_precision*** & T_Im = T_FFTData->Im;

 FFTDataContainer_blocks* Tf_FFTData = static_cast<FFTDataContainer_blocks*>(TransferFunction);
 FFT_precision*** & Tf_Re = Tf_FFTData->Re;
 FFT_precision*** & Tf_Im = Tf_FFTData->Im;

 int i,b;
 int currentBlock;
 int lastBlock = _nbBlocks-1; // could be reduced for far way interaction
 const int end_i = _fft_nx * _fft_ny; //size of a block
 FFT_precision * t_re;
 FFT_precision * t_im;
 FFT_precision * s_re;
 FFT_precision * s_im;
 FFT_precision * tr_re;
 FFT_precision * tr_im;


 for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
 t_re = &(T_Re[currentBlock][0][0]);
 t_im = &(T_Im[currentBlock][0][0]);
 for(b=0; b<=lastBlock-currentBlock; b++) {
 s_re = &(S_Re[b][0][0]);
 s_im = &(S_Im[b][0][0]);
 tr_re = &(Tf_Re[b+currentBlock][0][0]);
 tr_im = &(Tf_Im[b+currentBlock][0][0]);


 for(i=0; i<end_i; ++i) {
 t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
 t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
 }
 }
 }


 }

 void FFTAcceleration_blocks::FFT_M2L_vec(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction)
 {
 FFTDataContainer_blocks* S_FFTData = getFFTData(Source);
 FFT_precision*** & S_Re = S_FFTData->Re;
 FFT_precision*** & S_Im = S_FFTData->Im;

 FFTDataContainer_blocks* T_FFTData = getFFTData(Target);
 FFT_precision*** & T_Re = T_FFTData->Re;
 FFT_precision*** & T_Im = T_FFTData->Im;

 FFTDataContainer_blocks* Tf_FFTData = static_cast<FFTDataContainer_blocks*>(TransferFunction);
 FFT_precision*** & Tf_Re = Tf_FFTData->Re;
 FFT_precision*** & Tf_Im = Tf_FFTData->Im;

 int i,b;
 int currentBlock;
 int lastBlock = _nbBlocks-1; // could be reduced for far way interaction
 const int end_i = _fft_nx * _fft_ny; //size of a block
 FFT_precision * t_re;
 FFT_precision * t_im;
 FFT_precision * s_re;
 FFT_precision * s_im;
 FFT_precision * tr_re;
 FFT_precision * tr_im;


 for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
 t_re = &(T_Re[currentBlock][0][0]);
 t_im = &(T_Im[currentBlock][0][0]);
 for(b=0; b<=lastBlock-currentBlock; b++) {
 s_re = &(S_Re[b][0][0]);
 s_im = &(S_Im[b][0][0]);
 tr_re = &(Tf_Re[b+currentBlock][0][0]);
 tr_im = &(Tf_Im[b+currentBlock][0][0]);

 #pragma ivdep
 #pragma vector aligned
 for(i=0; i<end_i; ++i) {
 t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
 t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
 }
 }
 }

 }



 void FFTAcceleration_blocks::FFT_M2L_OrderReduction(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction, int order)
 {
 FFTDataContainer_blocks* S_FFTData = getFFTData(Source);
 FFT_precision*** & S_Re = S_FFTData->Re;
 FFT_precision*** & S_Im = S_FFTData->Im;

 FFTDataContainer_blocks* T_FFTData = getFFTData(Target);
 FFT_precision*** & T_Re = T_FFTData->Re;
 FFT_precision*** & T_Im = T_FFTData->Im;

 FFTDataContainer_blocks* Tf_FFTData = static_cast<FFTDataContainer_blocks*>(TransferFunction);
 FFT_precision*** & Tf_Re = Tf_FFTData->Re;
 FFT_precision*** & Tf_Im = Tf_FFTData->Im;

 int i,b;
 int currentBlock;

 int lastBlock = (order/_nbLinePerBlock); //compute up to order
 if(order % _nbLinePerBlock == 0) lastBlock--;
 if(lastBlock < 0) lastBlock = 0;
 else if(lastBlock > _nbBlocks-1) lastBlock = _nbBlocks-1;

 const int end_i = _fft_nx * _fft_ny; //size of a block
 FFT_precision * t_re;
 FFT_precision * t_im;
 FFT_precision * s_re;
 FFT_precision * s_im;
 FFT_precision * tr_re;
 FFT_precision * tr_im;


 for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
 t_re = &(T_Re[currentBlock][0][0]);
 t_im = &(T_Im[currentBlock][0][0]);
 for(b=0; b<=lastBlock-currentBlock; b++) {
 s_re = &(S_Re[b][0][0]);
 s_im = &(S_Im[b][0][0]);
 tr_re = &(Tf_Re[b+currentBlock][0][0]);
 tr_im = &(Tf_Im[b+currentBlock][0][0]);

 #pragma ivdep
 #pragma vector aligned
 for(i=0; i<end_i; ++i) {
 t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
 t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
 }
 }
 }

 }
 */
