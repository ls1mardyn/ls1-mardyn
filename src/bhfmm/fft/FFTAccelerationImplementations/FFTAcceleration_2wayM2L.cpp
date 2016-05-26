/*
 * FFTAcceleration_2wayM2L.cpp
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_2wayM2L.h"


FFTDataContainer_arrays* FFTAcceleration_2wayM2L::getFFTData(FFTAccelerableExpansion & Expansion)
{
  FFTDataContainer_arrays* FFTData = NULL;
  if(! Expansion.issetFFTData()) {
    FFTData = new FFTDataContainer_arrays(_totalSize);
    FFTData->Re = alloc_aligned_array(_totalSize);
    FFTData->Im = alloc_aligned_array(_totalSize);
    Expansion._FFTData = FFTData;
  }
  else
  {
    FFTData = static_cast<FFTDataContainer_arrays*>(Expansion._FFTData);
  }
  
  return FFTData;
}

void FFTAcceleration_2wayM2L::FFT_initialize_Target(FFTAccelerableExpansion & Expansion) 
{
  FFTDataContainer_arrays* FFTData = getFFTData(Expansion);
  
  clear_array(FFTData->Re, _totalSize);
  clear_array(FFTData->Im, _totalSize);
}

void FFTAcceleration_2wayM2L::FFT_M2L_2way(FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2, 
      FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction)
{   
  FFTDataContainer_arrays* Source1Data = getFFTData(Source1);
  FFT_precision* & S1_Re = Source1Data->Re;
  FFT_precision* & S1_Im = Source1Data->Im;
  
  FFTDataContainer_arrays* Source2Data = getFFTData(Source2);
  FFT_precision* & S2_Re = Source2Data->Re;
  FFT_precision* & S2_Im = Source2Data->Im;
  
  FFTDataContainer_arrays* Target1Data = getFFTData(Target1);
  FFT_precision* & T1_Re = Target1Data->Re;
  FFT_precision* & T1_Im = Target1Data->Im;
  
  FFTDataContainer_arrays* Target2Data = getFFTData(Target2);
  FFT_precision* & T2_Re = Target2Data->Re;
  FFT_precision* & T2_Im = Target2Data->Im;

  FFTDataContainer_arrays* TFData = static_cast<FFTDataContainer_arrays*>(TransferFunction);
  FFT_precision* & Tr_Re = TFData->Re;
  FFT_precision* & Tr_Im = TFData->Im;
  
  int i,j;
  const int end_i = _totalSize;
  const int shift_i = end_i / 2;
/*
 * Tr(r)[i][j] = -1^(p+1) Tr(-r)[i+max_i/2][j]
 * see doc/2WayM2L.txt
 */
  if(_p%2 == 0) { 
    j=shift_i;
    for(i=0; i<shift_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * -Tr_Im[i] - S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * -Tr_Re[i] + S2_Im[j] * Tr_Im[i];
      ++j;
    }
    j=0;
    for(;i<end_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * -Tr_Im[i] - S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * -Tr_Re[i] + S2_Im[j] * Tr_Im[i];
      ++j;
    }
  } else {

    j=shift_i;
    for(i=0; i<shift_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * Tr_Im[i] + S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * Tr_Re[i] - S2_Im[j] * Tr_Im[i];
      ++j;
    }
    j=0;
    for(;i<end_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * Tr_Im[i] + S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * Tr_Re[i] - S2_Im[j] * Tr_Im[i];
      ++j;
    }
  }
}

void FFTAcceleration_2wayM2L::FFT_M2L_2way_vec(FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2, 
      FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction)
{   
  FFTDataContainer_arrays* Source1Data = getFFTData(Source1);
  FFT_precision* & S1_Re = Source1Data->Re;
  FFT_precision* & S1_Im = Source1Data->Im;
  
  FFTDataContainer_arrays* Source2Data = getFFTData(Source2);
  FFT_precision* & S2_Re = Source2Data->Re;
  FFT_precision* & S2_Im = Source2Data->Im;
  
  FFTDataContainer_arrays* Target1Data = getFFTData(Target1);
  FFT_precision* & T1_Re = Target1Data->Re;
  FFT_precision* & T1_Im = Target1Data->Im;
  
  FFTDataContainer_arrays* Target2Data = getFFTData(Target2);
  FFT_precision* & T2_Re = Target2Data->Re;
  FFT_precision* & T2_Im = Target2Data->Im;

  FFTDataContainer_arrays* TFData = static_cast<FFTDataContainer_arrays*>(TransferFunction);
  FFT_precision* & Tr_Re = TFData->Re;
  FFT_precision* & Tr_Im = TFData->Im;
  
  int i,j;
  const int end_i = _totalSize;
  const int shift_i = end_i / 2;
/*
 * Tr(r)[i][j] = -1^(p+1) Tr(-r)[i+max_i/2][j]
 * see doc/2WayM2L.txt
 */
  if(_p%2 == 0) { 
    j=shift_i;
	#pragma omp simd aligned (T1_Re, T1_Im, S1_Re, S1_Im, \
							  T2_Re, T2_Im, S2_Re, S2_Im, \
							  Tr_Re, Tr_Im: __FFT_MATRIX_ALIGNMENT__)
    for(i=0; i<shift_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * -Tr_Im[i] - S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * -Tr_Re[i] + S2_Im[j] * Tr_Im[i];
      ++j;
    }
    j=0;
	#pragma omp simd aligned (T1_Re, T1_Im, S1_Re, S1_Im, \
							  T2_Re, T2_Im, S2_Re, S2_Im, \
							  Tr_Re, Tr_Im: __FFT_MATRIX_ALIGNMENT__)
    for(i=shift_i; i<end_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * -Tr_Im[i] - S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * -Tr_Re[i] + S2_Im[j] * Tr_Im[i];
      ++j;
    }
  } else {

    j=shift_i;
	#pragma omp simd aligned (T1_Re, T1_Im, S1_Re, S1_Im, \
							  T2_Re, T2_Im, S2_Re, S2_Im, \
							  Tr_Re, Tr_Im: __FFT_MATRIX_ALIGNMENT__)
    for(i=0; i<shift_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * Tr_Im[i] + S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * Tr_Re[i] - S2_Im[j] * Tr_Im[i];
      ++j;
    }
    j=0;
	#pragma omp simd aligned (T1_Re, T1_Im, S1_Re, S1_Im, \
							  T2_Re, T2_Im, S2_Re, S2_Im, \
							  Tr_Re, Tr_Im: __FFT_MATRIX_ALIGNMENT__)
    for(i=shift_i;i<end_i; ++i) {
      T2_Re[i] += S1_Re[i] * Tr_Re[i] - S1_Im[i] * Tr_Im[i];
      T2_Im[i] += S1_Re[i] * Tr_Im[i] + S1_Im[i] * Tr_Re[i];
      T1_Im[j] += S2_Re[j] * Tr_Im[i] + S2_Im[j] * Tr_Re[i];
      T1_Re[j] += S2_Re[j] * Tr_Re[i] - S2_Im[j] * Tr_Im[i];
      ++j;
    }
  }
}
