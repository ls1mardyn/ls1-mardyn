/*
 * optimizedFFT.h
 *
 *  Created on: Sep 14, 2015
 *      Author: gallardjm
 */



#ifndef OPT_FFT_H_
#define OPT_FFT_H_

#include <stdexcept>

#include "bhfmm/fft/tools/optimizedFFT/basic/uhfft_r.h"
#include "bhfmm/fft/tools/optimizedFFT/basic/uhfft_c.h"

/**
 * Perform in place forward FFT on the matrices using Kurzak's optimisation
 * Throws Invalid_Argument if p not in [6,16]
 * @param Real, real part
 * @param Imag, imag part
 * @param p = SH_Expansion's order+1
 */
void optimizedFFT(FFT_precision** & Real, FFT_precision** & Imag, const int p);

/**
 * Perform in place backward FFT on the matrices using Kurzak's optimisation
 * Throws Invalid_Argument if p not in [6,16]
 * @param Real, Fourier space real part
 * @param Imag, Fourier space imag part
 * @param p = SH_Expansion's order+1
 */
void optimizedIFFT(FFT_precision** & Real, FFT_precision** & Imag, const int p);

#endif
