/*
 * BasicOptFFT.h
 *
 *  Created on: Mar 06, 2015
 *      Author: gallardjm
 */

#ifndef BASICOPTFFT_H_
#define BASICOPTFFT_H_


#include "bhfmm/fft/tools/optimizedFFT/optFFT_API.h"
#include "bhfmm/fft/tools/optimizedFFT/basic/optimizedFFT.h"

/**
 * API using Kurzak's code, require size_x = 2*size_y and 5<size_y<17
 * 
 * Use the code in the basic folder
 */
class BasicOptFFT : public optFFT_API 
{
  public:
  
    BasicOptFFT() {}
    ~BasicOptFFT() {}
  
    void optimizedFFT (FFT_precision** & Real, FFT_precision** & Imag, const int size_x, const int size_y) 
    { ::optimizedFFT(Real, Imag, size_y); }
    
    void optimizedIFFT(FFT_precision** & Real, FFT_precision** & Imag, const int size_x, const int size_y) 
    { ::optimizedIFFT(Real, Imag, size_y); }
};

#endif
