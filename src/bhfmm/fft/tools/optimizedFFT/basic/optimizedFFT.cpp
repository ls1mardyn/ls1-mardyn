/*
 * optimizedFFT.cpp
 *
 *  Created on: Sep 14, 2015
 *      Author: gallardjm
 */

#include "optimizedFFT.h"


/*
 * Compute forward FFT
 * p is the size of the expansion (order+1)
 */
void optimizedFFT(FFT_precision** & Real, FFT_precision** & Imag, const int p)
{
  int n,m;
  
  switch (p) {
    case 6:
      for (n = 0; n < p; n++)
        MFFTR12(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC12(&Real[0][m],&Imag[0][m]);
    break;
    case 7:
      for (n = 0; n < p; n++)
        MFFTR14(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC14(&Real[0][m],&Imag[0][m]);
    break;
    case 8:
      for (n = 0; n < p; n++)
        MFFTR16(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC16(&Real[0][m],&Imag[0][m]);
    break;
    case 9:
      for (n = 0; n < p; n++)
        MFFTR18(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC18(&Real[0][m],&Imag[0][m]);
    break;
    case 10:
      for (n = 0; n < p; n++)
        MFFTR20(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC20(&Real[0][m],&Imag[0][m]);
    break;
    case 11:
      for (n = 0; n < p; n++)
        MFFTR22(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC22(&Real[0][m],&Imag[0][m]);
    break;
    case 12:
      for (n = 0; n < p; n++)
        MFFTR24(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC24(&Real[0][m],&Imag[0][m]);
    break;
    case 13:
      for (n = 0; n < p; n++)
        MFFTR26(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC26(&Real[0][m],&Imag[0][m]);
    break;
    case 14:
      for (n = 0; n < p; n++)
        MFFTR28(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC28(&Real[0][m],&Imag[0][m]);
    break;
    case 15:
      for (n = 0; n < p; n++)
        MFFTR30(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC30(&Real[0][m],&Imag[0][m]);
    break;
    case 16:
      for (n = 0; n < p; n++)
        MFFTR32(&Real[n][0],&Imag[n][0]);
      for (m = 0; m < p; m++)
        MFFTC32(&Real[0][m],&Imag[0][m]);
    break;
    default:
      throw std::invalid_argument("order must be between 5 and 15" );
    break;
  }
}

/*
 * Compute backward FFT
 * p is the size of the expansion (order+1)
 */
void optimizedIFFT(FFT_precision** & Real, FFT_precision** & Imag, const int p)
{
  int n,m;
  
  switch (p) {
    case 6:
      for (m = 0; m < p; m++)
        MIFFTC12(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR12(&Real[n][0], &Imag[n][0]);
    break;
    case 7:
      for (m = 0; m < p; m++)
        MIFFTC14(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR14(&Real[n][0], &Imag[n][0]);
    break;
    case 8:
      for (m = 0; m < p; m++)
        MIFFTC16(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR16(&Real[n][0], &Imag[n][0]);
    break;
    case 9:
      for (m = 0; m < p; m++)
        MIFFTC18(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR18(&Real[n][0], &Imag[n][0]);
    break;
    case 10:
      for (m = 0; m < p; m++)
        MIFFTC20(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR20(&Real[n][0], &Imag[n][0]);
    break;
    case 11:
      for (m = 0; m < p; m++)
        MIFFTC22(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR22(&Real[n][0], &Imag[n][0]);
    break;
    case 12:
      for (m = 0; m < p; m++)
        MIFFTC24(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR24(&Real[n][0], &Imag[n][0]);
    break;
    case 13:
      for (m = 0; m < p; m++)
        MIFFTC26(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR26(&Real[n][0], &Imag[n][0]);
    break;
    case 14:
      for (m = 0; m < p; m++)
        MIFFTC28(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR28(&Real[n][0], &Imag[n][0]);
    break;
    case 15:
      for (m = 0; m < p; m++)
        MIFFTC30(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR30(&Real[n][0], &Imag[n][0]);
    break;
    case 16:
      for (m = 0; m < p; m++)
        MIFFTC32(&Real[0][m], &Imag[0][m]);
      for (n = 0; n < 2*p; n++)
        MIFFTR32(&Real[n][0], &Imag[n][0]);
    break;
    default:
      throw std::invalid_argument("Order must be between 5 and 15" );
    break;
  }

  for (n = 0; n < 2*p; n++)
    for (m = 0; m < p; m++)
    {
      Real[n][m] /= (double)(p*p*4);
      Imag[n][m] /= (double)(p*p*4);
    }
    
}
