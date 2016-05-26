/*
 * FFTAcceleration_blocks_optFFT.cpp
 *
 *  Created on: Mar 08, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_blocks_optFFT.h"


FFTAcceleration_blocks_optFFT::FFTAcceleration_blocks_optFFT(int order)
{
  //printf("FFTAcceleration_blocks_optFFT \n");
  
  _p = order+1;
  _nbLinePerBlock = 4;
  if( _p % _nbLinePerBlock != 0)
  {
    throw std::invalid_argument("FFTAcceleration: _p % _nbLinePerBlock != 0");
  }
  
  _nbBlocks = _p / _nbLinePerBlock;
  _fft_nx = 2*_nbLinePerBlock; 
  if(!FFTSettings::USE_2WAY_M2L) //requires an odd number of line in 2way
    _fft_nx--; // litterature uses 2*_nbLinePerBlock (no -1)
  _fft_ny = _p;
  
  _optFFT_API = optFFT_API_Factory::getOptFFT_API(order, false);
}


void FFTAcceleration_blocks_optFFT::FFT_initialize_Source(FFTAccelerableExpansion & Expansion, double radius) 
{
  FFTDataContainer_blocks* FFTData = getFFTData(Expansion);
  FFT_precision*** & Re = FFTData->Re;
  FFT_precision*** & Im = FFTData->Im;
  
  int b,i,j,trueI;
  FFT_precision scaling = 1.0/(FFT_precision)radius; //include radius scaling and -1^i
  
  for(b=0; b<_nbBlocks; b++)
  {
    for(i=0; i<_nbLinePerBlock; i++) {
      trueI = i + _nbLinePerBlock*b;
      for(j=0; j<=trueI; j++) {
        Re[b][i][j] =  (FFT_precision)Expansion.get_C(trueI,j) * scaling;
        Im[b][i][j] = -(FFT_precision)Expansion.get_S(trueI,j) * scaling; //we want to use use conj(Source) 
      } 
      for(;j<_fft_ny;j++) {
        Re[b][i][j] = 0.0;
        Im[b][i][j] = 0.0;
      }
      scaling /= -(FFT_precision)radius;
    }
    
    for(;i<_fft_nx;i++)
      for(j=0; j<_fft_ny; j++) {
        Re[b][i][j] = 0.0;
        Im[b][i][j] = 0.0;
      }
    
    _optFFT_API->optimizedFFT(Re[b], Im[b], _fft_nx, _fft_ny);
  }

}


void FFTAcceleration_blocks_optFFT::FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion)
{
  FFTDataContainer_blocks* FFTData = getFFTData(Expansion);
  FFT_precision*** & Re = FFTData->Re;
  FFT_precision*** & Im = FFTData->Im;
  
  int b,i,j,trueI;
  FFT_precision minus_one_power_j = 1.0;

  for(b=0; b<_nbBlocks; b++)
  {
    for(i=0; i<_nbLinePerBlock; i++) {
      trueI = i + _nbLinePerBlock*b;
       minus_one_power_j = 1.0;
      for(j=0; j<=trueI; j++) {
        Re[b][_nbLinePerBlock-i-1][j] = (FFT_precision)Expansion.get_C(trueI,j) * minus_one_power_j;
        Im[b][_nbLinePerBlock-i-1][j] = (FFT_precision)Expansion.get_S(trueI,j) * -minus_one_power_j;
        minus_one_power_j *= -1.0;
      }
      
      for(;j<_fft_ny;j++) {
        Re[b][_nbLinePerBlock-i-1][j] = 0.0;
        Im[b][_nbLinePerBlock-i-1][j] = 0.0;
      }
    }
    
    for(;i<_fft_nx;i++) 
      for(j=0; j<_fft_ny; j++) {
        Re[b][i][j] = 0.0;
        Im[b][i][j] = 0.0;
      }
 
    _optFFT_API->optimizedFFT(Re[b], Im[b], _fft_nx, _fft_ny);
  }
}


void FFTAcceleration_blocks_optFFT::FFT_finalize_Target(FFTAccelerableExpansion & Expansion, double radius)
{
  FFTDataContainer_blocks* FFTData = getFFTData(Expansion);
  FFT_precision*** & Re = FFTData->Re;
  FFT_precision*** & Im = FFTData->Im;  
  
  int b,i,j,trueI;
  double scaling  = 1.0;// / (FFT_precision)(_fft_nx * _fft_ny);
  double scaling2 = 1.0;// / (FFT_precision)(_fft_nx * _fft_ny);
  double minus_one_power_j;

  for(b=0;b<_nbBlocks;b++) {
    
    _optFFT_API->optimizedIFFT(Re[b], Im[b], _fft_nx, _fft_ny);
    
    for (i=0; i<_nbLinePerBlock; i++) 
    {
      trueI = i + _nbLinePerBlock*b;
      minus_one_power_j = 1.0;
      for (j=0; j<=trueI; j++) 
      {
        Expansion.get_C(trueI,j) +=  minus_one_power_j * (double)(Re[b][_nbLinePerBlock-i-1][j]) * scaling;
        Expansion.get_S(trueI,j) += -minus_one_power_j * (double)(Im[b][_nbLinePerBlock-i-1][j]) * scaling;
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
          Expansion.get_C(trueI,j) +=  minus_one_power_j * (double)(Re[b][_nbLinePerBlock*2-i-1][j]) * scaling2;
          Expansion.get_S(trueI,j) += -minus_one_power_j * (double)(Im[b][_nbLinePerBlock*2-i-1][j]) * scaling2;
          minus_one_power_j *= -1.0;
        }
        scaling2 /= radius;
      }
    }
    scaling2 /= (FFT_precision)radius;
  }
}
