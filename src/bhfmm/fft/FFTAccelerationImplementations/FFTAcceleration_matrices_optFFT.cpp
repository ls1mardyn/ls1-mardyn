/*
 * FFTAcceleration_matrices_optFFT.cpp
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#include "FFTAcceleration_matrices_optFFT.h"

FFTAcceleration_matrices_optFFT::FFTAcceleration_matrices_optFFT(int order)
{
  _p = order+1;
  _fft_nx = 2*_p-1; //Litterature uses 2*_p but -1 works, however will force fakedOptFFT instead of Kurzak's FFT (see optFFT_API_Factory)
  _fft_ny = _p;
  
  _optFFT_API = optFFT_API_Factory::getOptFFT_API(order, (_fft_nx == 2*_fft_ny));
}


void FFTAcceleration_matrices_optFFT::FFT_initialize_Source(FFTAccelerableExpansion & Expansion, double radius) 
{
  FFTDataContainer_matrices* FFTData = getFFTData(Expansion);
  FFT_precision** & Re = FFTData->Re;
  FFT_precision** & Im = FFTData->Im;
  
  int n, m;
  FFT_precision minus_one_power_n = 1.0;
  FFT_precision r = (FFT_precision)radius;

  for (n = 0; n < _p; n++) {
    for (m = 0; m <= n; m++)
    {
      Re[n][m] =  minus_one_power_n  * (FFT_precision)Expansion.get_C(n,m) / r;
      Im[n][m] =  -minus_one_power_n * (FFT_precision)Expansion.get_S(n,m) / r; //we want to use use conj(Source)
    }
    for (; m<_fft_ny; m++)
    {
      Re[n][m] = 0.0;
      Im[n][m] = 0.0;
    }
    minus_one_power_n *= -1.0;
    r *= (FFT_precision)radius;
  }
  
  for(; n<_fft_nx; n++)
   for (m = 0; m<_fft_ny; m++)
    {
      Re[n][m] = 0.0;
      Im[n][m] = 0.0;
    }

  _optFFT_API->optimizedFFT(Re, Im, _fft_nx, _fft_ny);  
}


void FFTAcceleration_matrices_optFFT::FFT_initialize_TransferFunction(FFTAccelerableExpansion & Expansion)
{
  FFTDataContainer_matrices* FFTData = getFFTData(Expansion);
  FFT_precision** & Re = FFTData->Re;
  FFT_precision** & Im = FFTData->Im;
  int n, m;
  FFT_precision minus_one_power_m = 1.0;

  for (n = 0; n < _p; n++) {
    minus_one_power_m = 1.0;
    for (m = 0; m <= n; m++)
    {
      //SH_Expansion.signed_acc_const_CS(n,-m,Re[_p-n-1][m],Im[_p-n-1][m]); //flipmatrix included
      Re[_p-n-1][m] =  (FFT_precision)Expansion.get_C(n,m) * minus_one_power_m;
      Im[_p-n-1][m] =  (FFT_precision)Expansion.get_S(n,m) * -minus_one_power_m;
      minus_one_power_m *= -1.0;
    }
    for (; m < _fft_ny; m++)
    {
      Re[_p-n-1][m] = 0;
      Im[_p-n-1][m] = 0;
    }
  }
  
  for (; n < _fft_nx; n++)
    for (m = 0 ; m < _fft_ny; m++)
    {
      Re[n][m] = 0;
      Im[n][m] = 0;
    }
      
  _optFFT_API->optimizedFFT(Re, Im, _fft_nx, _fft_ny);
}

    
void FFTAcceleration_matrices_optFFT::FFT_finalize_Target(FFTAccelerableExpansion & Expansion, double radius)
{
  FFTDataContainer_matrices* FFTData = getFFTData(Expansion);
  FFT_precision** & Re = FFTData->Re;
  FFT_precision** & Im = FFTData->Im;
  int n, m;

  _optFFT_API->optimizedIFFT(Re, Im, _fft_nx, _fft_ny);
  
  FFT_precision scaling = 1.0;
  double minus_1_power_m;
  for (n=0; n<_p; n++) 
  {
    minus_1_power_m = 1.0;
    for (m=0; m<=n; m++) 
    {
      Expansion.get_C(n,m) += (double)(Re[(_p-n-1)][m] * scaling) * minus_1_power_m;
      Expansion.get_S(n,m) += (double)(Im[(_p-n-1)][m] * scaling) * -minus_1_power_m;
      minus_1_power_m *= -1.0;
    }
    scaling /= (FFT_precision)radius;
  }

}

