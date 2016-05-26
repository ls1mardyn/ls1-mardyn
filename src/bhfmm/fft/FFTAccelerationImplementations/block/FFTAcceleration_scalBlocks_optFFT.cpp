/*
 * FFTAcceleration_scalBlocks_optFFT.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_scalBlocks_optFFT.h"

FFTAcceleration_scalBlocks_optFFT::FFTAcceleration_scalBlocks_optFFT(int order) 
{
  //printf("FFTAcceleration_scalBlocks_optFFT\n");
  _p = order+1;
  _nbLinePerBlock = 4;
  if( _p % _nbLinePerBlock != 0)
  {
    throw std::invalid_argument("FFTAcceleration: _p % _nbLinePerBlock != 0");
  }
  
  _nbBlocks = _p / _nbLinePerBlock;
  _fft_nx = 2*_nbLinePerBlock; // litterature uses 2*_nbLinePerBlock (no -1)
  _fft_ny = _p;
  
  _Re_tmp = alloc_matrix(_fft_nx, _fft_ny);
  _Im_tmp = alloc_matrix(_fft_nx, _fft_ny);
  
  _optFFT_API = optFFT_API_Factory::getOptFFT_API(order, false);
  _blockSize = new int[_nbBlocks];
  
  //initializing _blockSize
  int repetition, i;
  
  for(i=0; i<_nbBlocks/2; i++) { //first half should have smaller blocks
    repetition = _fft_ny / (_nbLinePerBlock * (i+1)); //will be at least 2
    while(_fft_ny % repetition != 0)
      repetition--;
    _blockSize[i] = _fft_ny / repetition;
  }
  for(;i<_nbBlocks; i++)
    _blockSize[i] = _fft_ny;
}

FFTDataContainer_scalBlocks* FFTAcceleration_scalBlocks_optFFT::getFFTData(FFTAccelerableExpansion & Expansion)
{
  FFTDataContainer_scalBlocks* FFTData = NULL;
  if(! Expansion.issetFFTData()) {
    FFTData = new FFTDataContainer_scalBlocks(_nbBlocks, _fft_nx, _fft_ny, false, _blockSize);
    FFTData->Re = alloc_blocks_arr(_nbBlocks, _fft_nx, _fft_ny);
    FFTData->Im = alloc_blocks_arr(_nbBlocks, _fft_nx, _fft_ny);
    Expansion._FFTData = FFTData;
  }
  else
  {
    FFTData = static_cast<FFTDataContainer_scalBlocks*>(Expansion._FFTData);
  }
  
  return FFTData;
}

FFTDataContainer_scalBlocks* FFTAcceleration_scalBlocks_optFFT::getFFTData_scal(FFTAccelerableExpansion & Expansion)
{
  FFTDataContainer_scalBlocks* FFTData = NULL;
  if(! Expansion.issetFFTData()) {
    FFTData = new FFTDataContainer_scalBlocks(_nbBlocks, _fft_nx, _fft_ny, true, _blockSize);
    FFTData->Re = alloc_scalBlocks_arr(_nbBlocks, _fft_nx, _blockSize);
    FFTData->Im = alloc_scalBlocks_arr(_nbBlocks, _fft_nx, _blockSize);
    Expansion._FFTData = FFTData;
  }
  else
  {
    FFTData = static_cast<FFTDataContainer_scalBlocks*>(Expansion._FFTData);
  }
  
  return FFTData;
}

void FFTAcceleration_scalBlocks_optFFT::FFT_initialize_Target(FFTAccelerableExpansion & Expansion)
{
  FFTDataContainer_scalBlocks* FFTData = getFFTData(Expansion);
  
  clear_blocks_arr(FFTData->Re, _nbBlocks, _fft_nx, _fft_ny);
  clear_blocks_arr(FFTData->Im, _nbBlocks, _fft_nx, _fft_ny);
}

void FFTAcceleration_scalBlocks_optFFT::FFT_initialize_Source(FFTAccelerableExpansion & Expansion, double radius)
{
  FFTDataContainer_scalBlocks* FFTData = getFFTData(Expansion);
  FFT_precision** & Re_arr = FFTData->Re;
  FFT_precision** & Im_arr = FFTData->Im;
  
  int b,i,j,trueI;
  FFT_precision scaling = 1.0/(FFT_precision)radius; //include radius scaling and -1^i
  
  for(b=0; b<_nbBlocks; b++)
  {
    for(i=0; i<_nbLinePerBlock; i++) {
      trueI = i + _nbLinePerBlock*b;
      for(j=0; j<=trueI; j++) {
        _Re_tmp[i][j] =  (FFT_precision)Expansion.get_C(trueI,j) * scaling;
        _Im_tmp[i][j] = -(FFT_precision)Expansion.get_S(trueI,j) * scaling; //we want to use use conj(Source) 
      } 
      for(;j<_fft_ny;j++) {
        _Re_tmp[i][j] = 0.0;
        _Im_tmp[i][j] = 0.0;
      }
      scaling /= -(FFT_precision)radius;
    }
    
    for(;i<_fft_nx;i++)
      for(j=0; j<_fft_ny; j++) {
        _Re_tmp[i][j] = 0.0;
        _Im_tmp[i][j] = 0.0;
      }
    
    _optFFT_API->optimizedFFT(_Re_tmp, _Im_tmp, _fft_nx, _fft_ny);
    
    for(i=0;i<_fft_nx;i++)
      for(j=0;j<_fft_ny;j++) {
        Re_arr[b][j*_fft_nx+i] = _Re_tmp[i][j];
        Im_arr[b][j*_fft_nx+i] = _Im_tmp[i][j];
      }
  }
}

void FFTAcceleration_scalBlocks_optFFT::FFT_finalize_Target(FFTAccelerableExpansion & Expansion, double radius)
{
  FFTDataContainer_scalBlocks* FFTData = getFFTData(Expansion);
  FFT_precision** & Re_arr = FFTData->Re;
  FFT_precision** & Im_arr = FFTData->Im;  
  
  int b,i,j,trueI;
  double scaling  = 1.0;// / (FFT_precision)(_fft_nx * _fft_ny);
  double scaling2 = 1.0;// / (FFT_precision)(_fft_nx * _fft_ny);
  double minus_one_power_j;

  for(b=0;b<_nbBlocks;b++) {
    
    for(i=0;i<_fft_nx;i++)
      for(j=0;j<_fft_ny;j++) {
       _Re_tmp[i][j] =  Re_arr[b][j*_fft_nx+i];
       _Im_tmp[i][j] =  Im_arr[b][j*_fft_nx+i];
      }
    
    _optFFT_API->optimizedIFFT(_Re_tmp, _Im_tmp, _fft_nx, _fft_ny);
    
    for (i=0; i<_nbLinePerBlock; i++) 
    {
      trueI = i + _nbLinePerBlock*b;
      minus_one_power_j = 1.0;
      for (j=0; j<=trueI; j++) 
      {
        Expansion.get_C(trueI,j) +=  minus_one_power_j * (double)(_Re_tmp[_nbLinePerBlock-i-1][j]) * scaling;
        Expansion.get_S(trueI,j) += -minus_one_power_j * (double)(_Im_tmp[_nbLinePerBlock-i-1][j]) * scaling;
        minus_one_power_j *= -1.0;
      }
      scaling /= radius;
    }
    
    if(b!=0) {
      for (i=1; i<_nbLinePerBlock; i++) 
      {
        trueI = i + _nbLinePerBlock*(b-1);
        minus_one_power_j = 1.0;
        for (j=0; j<=trueI; j++) 
        {
          Expansion.get_C(trueI,j) +=  minus_one_power_j * (double)(_Re_tmp[_nbLinePerBlock*2-i-1][j]) * scaling2;
          Expansion.get_S(trueI,j) += -minus_one_power_j * (double)(_Im_tmp[_nbLinePerBlock*2-i-1][j]) * scaling2;
          minus_one_power_j *= -1.0;
        }
        scaling2 /= radius;
      }
    }
    scaling2 /= (FFT_precision)radius;
  }
}


template <bool Vect, bool OrderRed>
void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_template(
  FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction, int order)
{
  FFTDataContainer_scalBlocks* S_FFTData = getFFTData(Source);
  FFT_precision** & S_Re = S_FFTData->Re;
  FFT_precision** & S_Im = S_FFTData->Im;
  
  FFTDataContainer_scalBlocks* T_FFTData = getFFTData(Target);
  FFT_precision** & T_Re = T_FFTData->Re;
  FFT_precision** & T_Im = T_FFTData->Im;
  
  FFTDataContainer_scalBlocks* Tf_FFTData = static_cast<FFTDataContainer_scalBlocks*>(TransferFunction);
  FFT_precision** & Tf_Re = Tf_FFTData->Re;
  FFT_precision** & Tf_Im = Tf_FFTData->Im;
  
  int i,b;
  int currentBlock, tfBlock;
  int lastBlock;
  
  if(OrderRed) {
    lastBlock = (order/_nbLinePerBlock); //compute up to order
    if(order % _nbLinePerBlock == 0) lastBlock--;
    if(lastBlock < 0) lastBlock = 0;
    else if(lastBlock > _nbBlocks-1) lastBlock = _nbBlocks-1;
  } else {
    lastBlock = _nbBlocks-1; // could be reduced for far way interaction
  }
  
  int repetition;
  int j, j_rep, j_dec;
  const int end_i = _fft_nx * _fft_ny; //size of a block
  FFT_precision * t_re;
  FFT_precision * t_im;
  FFT_precision * s_re;
  FFT_precision * s_im;
  FFT_precision * tf_re;
  FFT_precision * tf_im;
  
  
  for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
    t_re = &(T_Re[currentBlock][0]);
    t_im = &(T_Im[currentBlock][0]);
    for(b=0; b<=lastBlock-currentBlock; b++) {
      tfBlock = b+currentBlock;
      s_re = &(S_Re[b][0]);
      s_im = &(S_Im[b][0]);
      tf_re = &(Tf_Re[tfBlock][0]);
      tf_im = &(Tf_Im[tfBlock][0]);
      if(tfBlock < _nbBlocks/2) { //scaling blocks scheme
        repetition = _fft_ny / _blockSize[tfBlock];
        for(j=0; j<_blockSize[tfBlock];j++) {
          j_dec = j*_fft_nx;
          j_rep = repetition * j_dec;
          if(Vect) {
			#pragma omp simd aligned (t_re, t_im, s_re, s_im, \
									  tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
            for(i=0; i<_fft_nx; i++) {
              t_re[i+j_rep] += s_re[i+j_rep] * tf_re[i+j_dec] - s_im[i+j_rep] * tf_im[i+j_dec];
              t_im[i+j_rep] += s_re[i+j_rep] * tf_im[i+j_dec] + s_im[i+j_rep] * tf_re[i+j_dec];
            }
          } else {
            for(i=0; i<_fft_nx; i++) {
              t_re[i+j_rep] += s_re[i+j_rep] * tf_re[i+j_dec] - s_im[i+j_rep] * tf_im[i+j_dec];
              t_im[i+j_rep] += s_re[i+j_rep] * tf_im[i+j_dec] + s_im[i+j_rep] * tf_re[i+j_dec];
            }
          }
        }
          
      }
      else { //same size block, default scheme
        if(Vect) {
		  #pragma omp simd aligned (t_re, t_im, s_re, s_im, \
		                            tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
          for(i=0; i<end_i; ++i) {
            t_re[i] += s_re[i] * tf_re[i] - s_im[i] * tf_im[i];
            t_im[i] += s_re[i] * tf_im[i] + s_im[i] * tf_re[i];
          }
        } else {
          for(i=0; i<end_i; ++i) {
            t_re[i] += s_re[i] * tf_re[i] - s_im[i] * tf_im[i];
            t_im[i] += s_re[i] * tf_im[i] + s_im[i] * tf_re[i];
          }
        }
      }
    }
  }
}

void FFTAcceleration_scalBlocks_optFFT::FFT_M2L(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction) 
{
  return FFT_M2L_template<false, false>(Source, Target, TransferFunction, 0);
}

void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_vec(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction) 
{
  return FFT_M2L_template<true, false>(Source, Target, TransferFunction, 0);
}

void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_OrderReduction(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction, int order)
{
  return FFT_M2L_template<false, true>(Source, Target, TransferFunction, order);
}

void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_OrderReduction_vec(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction, int order)
{
  return FFT_M2L_template<true, true>(Source, Target, TransferFunction, order);
}

template <bool Vect, bool OrderRed>
void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_2way_template(FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2, 
        FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction, int order)
{
  FFTDataContainer_scalBlocks* Source1Data = getFFTData(Source1);
  FFT_precision** & S1_Re = Source1Data->Re;
  FFT_precision** & S1_Im = Source1Data->Im;
  
  FFTDataContainer_scalBlocks* Source2Data = getFFTData(Source2);
  FFT_precision** & S2_Re = Source2Data->Re;
  FFT_precision** & S2_Im = Source2Data->Im;
  
  FFTDataContainer_scalBlocks* Target1Data = getFFTData(Target1);
  FFT_precision** & T1_Re = Target1Data->Re;
  FFT_precision** & T1_Im = Target1Data->Im;
  
  FFTDataContainer_scalBlocks* Target2Data = getFFTData(Target2);
  FFT_precision** & T2_Re = Target2Data->Re;
  FFT_precision** & T2_Im = Target2Data->Im;

  FFTDataContainer_scalBlocks* TFData = static_cast<FFTDataContainer_scalBlocks*>(TransferFunction);
  FFT_precision** & Tf_Re = TFData->Re;
  FFT_precision** & Tf_Im = TFData->Im;
  
  FFT_precision* t1_re;
  FFT_precision* t1_im;
  FFT_precision* s1_re;
  FFT_precision* s1_im;
  FFT_precision* t2_re;
  FFT_precision* t2_im;
  FFT_precision* s2_re;
  FFT_precision* s2_im;
  FFT_precision* tf_re;
  FFT_precision* tf_im;
  
  int i,b,n;
  int currentBlock, tfBlock;
  
  int lastBlock;
  
  if(OrderRed) {
    lastBlock = (order/_nbLinePerBlock); //compute up to order
    if(order % _nbLinePerBlock == 0) lastBlock--;
    if(lastBlock < 0) lastBlock = 0;
    else if(lastBlock > _nbBlocks-1) lastBlock = _nbBlocks-1;
  } else {
    lastBlock = _nbBlocks-1;
  }
  
  int repetition;
  int j, j_rep, j_dec;
  const int shift_2way = _nbLinePerBlock;
  const int end_shift_2way = 2*_nbLinePerBlock;
/*
 * Tr(r)[i][j] = -1^(p+1) Tr(-r)[i+max_i/2][j]
 * see doc/2WayM2L.txt
 */
  for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
    t1_re = &(T1_Re[currentBlock][0]);
    t1_im = &(T1_Im[currentBlock][0]);
    t2_re = &(T2_Re[currentBlock][0]);
    t2_im = &(T2_Im[currentBlock][0]);
    for(b=0; b<=lastBlock-currentBlock; b++) {
      tfBlock = b+currentBlock;
      s1_re = &(S1_Re[b][0]);
      s1_im = &(S1_Im[b][0]);
      s2_re = &(S2_Re[b][0]);
      s2_im = &(S2_Im[b][0]);
      tf_re = &(Tf_Re[tfBlock][0]);
      tf_im = &(Tf_Im[tfBlock][0]);
      
      if(tfBlock < _nbBlocks/2) { //scaling blocks scheme
        repetition = _fft_ny / _blockSize[tfBlock];
        for(j=0; j<_blockSize[tfBlock];j++) {
          j_dec = j*end_shift_2way;
          j_rep = repetition * j_dec;
          if(Vect) {
            #pragma omp simd aligned (t1_re, t1_im, s1_re, s1_im, \
                                      t2_re, t2_im, s2_re, s2_im, \
                                      tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
            for(i=0; i<shift_2way; i++) {
              t2_re[i+j_rep] += s1_re[i+j_rep] * tf_re[i+j_dec] - s1_im[i+j_rep] * tf_im[i+j_dec];
              t2_im[i+j_rep] += s1_re[i+j_rep] * tf_im[i+j_dec] + s1_im[i+j_rep] * tf_re[i+j_dec];
              t1_im[i+j_rep+shift_2way] -= s2_re[i+j_rep+shift_2way] * tf_im[i+j_dec] + s2_im[i+j_rep+shift_2way] * tf_re[i+j_dec];
              t1_re[i+j_rep+shift_2way] -= s2_re[i+j_rep+shift_2way] * tf_re[i+j_dec] - s2_im[i+j_rep+shift_2way] * tf_im[i+j_dec];
            }
            #pragma omp simd aligned (t1_re, t1_im, s1_re, s1_im, \
                                      t2_re, t2_im, s2_re, s2_im, \
                                      tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
            for(i=shift_2way; i<end_shift_2way; i++) {
              t2_re[i+j_rep] += s1_re[i+j_rep] * tf_re[i+j_dec] - s1_im[i+j_rep] * tf_im[i+j_dec];
              t2_im[i+j_rep] += s1_re[i+j_rep] * tf_im[i+j_dec] + s1_im[i+j_rep] * tf_re[i+j_dec];
              t1_im[i+j_rep-shift_2way] -= s2_re[i+j_rep-shift_2way] * tf_im[i+j_dec] + s2_im[i+j_rep-shift_2way] * tf_re[i+j_dec];
              t1_re[i+j_rep-shift_2way] -= s2_re[i+j_rep-shift_2way] * tf_re[i+j_dec] - s2_im[i+j_rep-shift_2way] * tf_im[i+j_dec];
            }
          } else {
            for(i=0; i<shift_2way; i++) {
              t2_re[i+j_rep] += s1_re[i+j_rep] * tf_re[i+j_dec] - s1_im[i+j_rep] * tf_im[i+j_dec];
              t2_im[i+j_rep] += s1_re[i+j_rep] * tf_im[i+j_dec] + s1_im[i+j_rep] * tf_re[i+j_dec];
              t1_im[i+j_rep+shift_2way] -= s2_re[i+j_rep+shift_2way] * tf_im[i+j_dec] + s2_im[i+j_rep+shift_2way] * tf_re[i+j_dec];
              t1_re[i+j_rep+shift_2way] -= s2_re[i+j_rep+shift_2way] * tf_re[i+j_dec] - s2_im[i+j_rep+shift_2way] * tf_im[i+j_dec];
            }
            for(;i<end_shift_2way; i++) {
              t2_re[i+j_rep] += s1_re[i+j_rep] * tf_re[i+j_dec] - s1_im[i+j_rep] * tf_im[i+j_dec];
              t2_im[i+j_rep] += s1_re[i+j_rep] * tf_im[i+j_dec] + s1_im[i+j_rep] * tf_re[i+j_dec];
              t1_im[i+j_rep-shift_2way] -= s2_re[i+j_rep-shift_2way] * tf_im[i+j_dec] + s2_im[i+j_rep-shift_2way] * tf_re[i+j_dec];
              t1_re[i+j_rep-shift_2way] -= s2_re[i+j_rep-shift_2way] * tf_re[i+j_dec] - s2_im[i+j_rep-shift_2way] * tf_im[i+j_dec];
            }
          }
        }
          
      }
      else { //same size block, default scheme
        if(Vect) {
          for(int j=0;j<_fft_ny;j++) {
            i=j*_fft_nx;
            #pragma omp simd aligned (t1_re, t1_im, s1_re, s1_im, \
                                      t2_re, t2_im, s2_re, s2_im, \
                                      tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
            for(n=0; n<shift_2way; n++) {
              t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
              t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
              t1_im[i+shift_2way] -= s2_re[i+shift_2way] * tf_im[i] + s2_im[i+shift_2way] * tf_re[i];
              t1_re[i+shift_2way] -= s2_re[i+shift_2way] * tf_re[i] - s2_im[i+shift_2way] * tf_im[i];
              ++i;
            }
            #pragma omp simd aligned (t1_re, t1_im, s1_re, s1_im, \
                                      t2_re, t2_im, s2_re, s2_im, \
                                      tf_re, tf_im: __FFT_MATRIX_ALIGNMENT__)
            for(n=shift_2way; n<end_shift_2way; n++) {
              t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
              t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
              t1_im[i-shift_2way] -= s2_re[i-shift_2way] * tf_im[i] + s2_im[i-shift_2way] * tf_re[i];
              t1_re[i-shift_2way] -= s2_re[i-shift_2way] * tf_re[i] - s2_im[i-shift_2way] * tf_im[i];
              ++i;
            }
          }
        } else {
          for(int j=0;j<_fft_ny;j++) {
            i=j*_fft_nx;
            for(n=0; n<shift_2way; n++) {
              t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
              t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
              t1_im[i+shift_2way] -= s2_re[i+shift_2way] * tf_im[i] + s2_im[i+shift_2way] * tf_re[i];
              t1_re[i+shift_2way] -= s2_re[i+shift_2way] * tf_re[i] - s2_im[i+shift_2way] * tf_im[i];
              ++i;
            }
            for(;n<end_shift_2way; n++) {
              t2_re[i] += s1_re[i] * tf_re[i] - s1_im[i] * tf_im[i];
              t2_im[i] += s1_re[i] * tf_im[i] + s1_im[i] * tf_re[i];
              t1_im[i-shift_2way] -= s2_re[i-shift_2way] * tf_im[i] + s2_im[i-shift_2way] * tf_re[i];
              t1_re[i-shift_2way] -= s2_re[i-shift_2way] * tf_re[i] - s2_im[i-shift_2way] * tf_im[i];
              ++i;
            }
          }
        }
      }
    }
  }
}

void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_2way(FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2, 
      FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction)
{
  return FFT_M2L_2way_template<false, false>(Source1, Source2, Target1, Target2, TransferFunction, 0);
}
      
void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_2way_vec(FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2, 
      FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction)
{
  return FFT_M2L_2way_template<true, false>(Source1, Source2, Target1, Target2, TransferFunction, 0);
}
    
void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_2way_ORed(FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2, 
      FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction, int order)
{
  return FFT_M2L_2way_template<false, true>(Source1, Source2, Target1, Target2, TransferFunction, order);
}
      
      
void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_2way_ORed_vec(FFTAccelerableExpansion & Source1, FFTAccelerableExpansion & Source2, 
      FFTAccelerableExpansion & Target1, FFTAccelerableExpansion & Target2, FFTDataContainer* TransferFunction, int order)
{
  return FFT_M2L_2way_template<true, true>(Source1, Source2, Target1, Target2, TransferFunction, order);
}

/*
void FFTAcceleration_scalBlocks_optFFT::FFT_M2L(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction)
{
  FFTDataContainer_scalBlocks* S_FFTData = getFFTData(Source);
  FFT_precision** & S_Re = S_FFTData->Re;
  FFT_precision** & S_Im = S_FFTData->Im;
  
  FFTDataContainer_scalBlocks* T_FFTData = getFFTData(Target);
  FFT_precision** & T_Re = T_FFTData->Re;
  FFT_precision** & T_Im = T_FFTData->Im;
  
  FFTDataContainer_scalBlocks* Tf_FFTData = static_cast<FFTDataContainer_scalBlocks*>(TransferFunction);
  FFT_precision** & Tf_Re = Tf_FFTData->Re;
  FFT_precision** & Tf_Im = Tf_FFTData->Im;
  
  int i,b;
  int currentBlock, tfBlock;
  int lastBlock = _nbBlocks-1; // could be reduced for far way interaction
  int repetition;
  int j, j_rep, j_dec;
  const int end_i = _fft_nx * _fft_ny; //size of a block
  FFT_precision * t_re;
  FFT_precision * t_im;
  FFT_precision * s_re;
  FFT_precision * s_im;
  FFT_precision * tr_re;
  FFT_precision * tr_im;
  
  
  for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
    t_re = &(T_Re[currentBlock][0]);
    t_im = &(T_Im[currentBlock][0]);
    for(b=0; b<=lastBlock-currentBlock; b++) {
      tfBlock = b+currentBlock;
      s_re = &(S_Re[b][0]);
      s_im = &(S_Im[b][0]);
      tr_re = &(Tf_Re[tfBlock][0]);
      tr_im = &(Tf_Im[tfBlock][0]);
      if(tfBlock < _nbBlocks/2) { //scaling blocks scheme
        repetition = _fft_ny / _blockSize[tfBlock];
        for(j=0; j<_blockSize[tfBlock];j++) {
          j_dec = j*_fft_nx;
          j_rep = repetition * j_dec;
          for(i=0; i<_fft_nx; i++) {
            t_re[i+j_rep] += s_re[i+j_rep] * tr_re[i+j_dec] - s_im[i+j_rep] * tr_im[i+j_dec];
            t_im[i+j_rep] += s_re[i+j_rep] * tr_im[i+j_dec] + s_im[i+j_rep] * tr_re[i+j_dec];
          }
        }
          
      }
      else { //same size block, default scheme
        for(i=0; i<end_i; ++i) {
          t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
          t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
        }
      }
    }
  }
}

void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_vec(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction)
{
  FFTDataContainer_scalBlocks* S_FFTData = getFFTData(Source);
  FFT_precision** & S_Re = S_FFTData->Re;
  FFT_precision** & S_Im = S_FFTData->Im;
  
  FFTDataContainer_scalBlocks* T_FFTData = getFFTData(Target);
  FFT_precision** & T_Re = T_FFTData->Re;
  FFT_precision** & T_Im = T_FFTData->Im;
  
  FFTDataContainer_scalBlocks* Tf_FFTData = static_cast<FFTDataContainer_scalBlocks*>(TransferFunction);
  FFT_precision** & Tf_Re = Tf_FFTData->Re;
  FFT_precision** & Tf_Im = Tf_FFTData->Im;
  
  int i,b;
  int currentBlock, tfBlock;
  int lastBlock = _nbBlocks-1; // could be reduced for far way interaction
  int repetition;
  int j, j_rep, j_dec;
  const int end_i = _fft_nx * _fft_ny; //size of a block
  FFT_precision * t_re;
  FFT_precision * t_im;
  FFT_precision * s_re;
  FFT_precision * s_im;
  FFT_precision * tr_re;
  FFT_precision * tr_im;
  
  
  for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
    t_re = &(T_Re[currentBlock][0]);
    t_im = &(T_Im[currentBlock][0]);
    for(b=0; b<=lastBlock-currentBlock; b++) {
      tfBlock = b+currentBlock;
      s_re = &(S_Re[b][0]);
      s_im = &(S_Im[b][0]);
      tr_re = &(Tf_Re[tfBlock][0]);
      tr_im = &(Tf_Im[tfBlock][0]);
      if(tfBlock < _nbBlocks/2) { //scaling blocks scheme
        repetition = _fft_ny / _blockSize[tfBlock];
        for(j=0; j<_blockSize[tfBlock];j++) {
          j_dec = j*_fft_nx;
          j_rep = repetition * j_dec;
#pragma ivdep
#pragma vector aligned
          for(i=0; i<_fft_nx; i++) {
            t_re[i+j_rep] += s_re[i+j_rep] * tr_re[i+j_dec] - s_im[i+j_rep] * tr_im[i+j_dec];
            t_im[i+j_rep] += s_re[i+j_rep] * tr_im[i+j_dec] + s_im[i+j_rep] * tr_re[i+j_dec];
          }
        }
          
      }
      else { //same size block, default scheme
#pragma ivdep
#pragma vector aligned
        for(i=0; i<end_i; ++i) {
          t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
          t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
        }
      }
    }
  }

  
}

void FFTAcceleration_scalBlocks_optFFT::FFT_M2L_OrderReduction(FFTAccelerableExpansion & Source, FFTAccelerableExpansion & Target, FFTDataContainer* TransferFunction, int order)
{
  FFTDataContainer_scalBlocks* S_FFTData = getFFTData(Source);
  FFT_precision** & S_Re = S_FFTData->Re;
  FFT_precision** & S_Im = S_FFTData->Im;
  
  FFTDataContainer_scalBlocks* T_FFTData = getFFTData(Target);
  FFT_precision** & T_Re = T_FFTData->Re;
  FFT_precision** & T_Im = T_FFTData->Im;
  
  FFTDataContainer_scalBlocks* Tf_FFTData = static_cast<FFTDataContainer_scalBlocks*>(TransferFunction);
  FFT_precision** & Tf_Re = Tf_FFTData->Re;
  FFT_precision** & Tf_Im = Tf_FFTData->Im;
  
  int i,b;
  int currentBlock, tfBlock;
  
  int lastBlock = (order/_nbLinePerBlock); //compute up to order
  if(order % _nbLinePerBlock == 0) lastBlock--;
  if(lastBlock < 0) lastBlock = 0;
  else if(lastBlock > _nbBlocks-1) lastBlock = _nbBlocks-1;
  
  int repetition;
  int j, j_rep, j_dec;
  const int end_i = _fft_nx * _fft_ny; //size of a block
  FFT_precision * t_re;
  FFT_precision * t_im;
  FFT_precision * s_re;
  FFT_precision * s_im;
  FFT_precision * tr_re;
  FFT_precision * tr_im;
  
  
  for(currentBlock=0;currentBlock<=lastBlock;currentBlock++) {
    t_re = &(T_Re[currentBlock][0]);
    t_im = &(T_Im[currentBlock][0]);
    for(b=0; b<=lastBlock-currentBlock; b++) {
      tfBlock = b+currentBlock;
      s_re = &(S_Re[b][0]);
      s_im = &(S_Im[b][0]);
      tr_re = &(Tf_Re[tfBlock][0]);
      tr_im = &(Tf_Im[tfBlock][0]);
      if(tfBlock < _nbBlocks/2) { //scaling blocks scheme
        repetition = _fft_ny / _blockSize[tfBlock];
        for(j=0; j<_blockSize[tfBlock];j++) {
          j_dec = j*_fft_nx;
          j_rep = repetition * j_dec;
#pragma ivdep
#pragma vector aligned
          for(i=0; i<_fft_nx; i++) {
            t_re[i+j_rep] += s_re[i+j_rep] * tr_re[i+j_dec] - s_im[i+j_rep] * tr_im[i+j_dec];
            t_im[i+j_rep] += s_re[i+j_rep] * tr_im[i+j_dec] + s_im[i+j_rep] * tr_re[i+j_dec];
          }
        }
          
      }
      else { //same size block, default scheme
#pragma ivdep
#pragma vector aligned
        for(i=0; i<end_i; ++i) {
          t_re[i] += s_re[i] * tr_re[i] - s_im[i] * tr_im[i];
          t_im[i] += s_re[i] * tr_im[i] + s_im[i] * tr_re[i];
        }
      }
    }
  }
}
*/


void FFTAcceleration_scalBlocks_optFFT::FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion)
{
  FFTDataContainer_scalBlocks* FFTData = getFFTData_scal(Expansion);
  FFT_precision** & Re_arr = FFTData->Re;
  FFT_precision** & Im_arr = FFTData->Im;
  
  int b,i,j,trueI;
  FFT_precision minus_one_power_j = 1.0;
  FFT_precision repetition;

  for(b=0; b<_nbBlocks; b++)
  {
    for(i=0; i<_nbLinePerBlock; i++) {
      trueI = i + _nbLinePerBlock*b;
       minus_one_power_j = 1.0;
      for(j=0; j<=trueI; j++) {
        _Re_tmp[_nbLinePerBlock-i-1][j] = (FFT_precision)Expansion.get_C(trueI,j) * minus_one_power_j;
        _Im_tmp[_nbLinePerBlock-i-1][j] = (FFT_precision)Expansion.get_S(trueI,j) * -minus_one_power_j;
        minus_one_power_j *= -1.0;
      }
      
      for(;j<_blockSize[b];j++) {
        _Re_tmp[_nbLinePerBlock-i-1][j] = 0.0;
        _Im_tmp[_nbLinePerBlock-i-1][j] = 0.0;
      }
    }
    
    for(;i<_fft_nx;i++) 
      for(j=0; j<_blockSize[b]; j++) {
        _Re_tmp[i][j] = 0.0;
        _Im_tmp[i][j] = 0.0;
      }
 
    _optFFT_API->optimizedFFT(_Re_tmp, _Im_tmp, _fft_nx, _blockSize[b]);
    
    if(b<_nbBlocks/2)
    {
      repetition = (FFT_precision)(_fft_ny/_blockSize[b]);
      for(i=0;i<_fft_nx;i++)
        for(j=0;j<_blockSize[b];j++){
          _Re_tmp[i][j] *= repetition;
          _Im_tmp[i][j] *= repetition;
        }
    }
    
    for(i=0;i<_fft_nx;i++)
      for(j=0;j<_blockSize[b];j++) {
        Re_arr[b][j*_fft_nx+i] = _Re_tmp[i][j];
        Im_arr[b][j*_fft_nx+i] = _Im_tmp[i][j];
      }
    
    
  }
}

