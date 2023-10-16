/*
 * FakeOptFFT.cpp
 *
 *  Created on: Mar 06, 2015
 *      Author: gallardjm
 */

#include "FakedOptFFT.h"

FFTW_API* FakedOptFFT::getFFTW_API(const int size_x, const int size_y) {
	pos p;
	p.x = size_x;
	p.y = size_y;
	if (_fftw_api_map.count(p) != 0) {
		return _fftw_api_map[p][mardyn_get_thread_num()];
	} else {
		FFTW_API** apis = new FFTW_API*[mardyn_get_max_threads()];
		for (int i = 0; i < mardyn_get_max_threads(); ++i) {
			apis[i] = new FFTW_API(size_x, size_y);
		}
		_fftw_api_map.insert(pair<pos, FFTW_API**>(p, apis));
		return apis[mardyn_get_thread_num()];
	}
}

/*
//Old version of the code, no prunning
void FakedOptFFT::optimizedFFT (FFT_precision** & Real, FFT_precision** & Imag, const int size_x, const int size_y)
{
  const int size_y_fft = size_y*2;
  FFTW_API* fftw_api = getFFTW_API(size_x, size_y_fft);
#if defined(__SINGLE_PRECISION_FFT__)
  fftwf_complex* in = fftw_api->getIn_Forward();
  fftwf_complex* out;
#else
  fftw_complex*  in = fftw_api->getIn_Forward();
  fftw_complex*  out;
#endif

  int i,j;
  FFT_precision minus_1_power_j;

  for(i=0;i<size_x;++i) {
    //each line looks like [x(0), .. , x(n), 0, x(-n), .. , x(-1)]
    //x(-i) = -1^i conj(x(i))
    for(j=0;j<size_y;++j) { //copy the first half
      in[i*size_y_fft+j][0] = Real[i][j];
      in[i*size_y_fft+j][1] = Imag[i][j];
    }
    in[i*size_y_fft+size_y][0] = 0.0; //add on zero pad
    in[i*size_y_fft+size_y][1] = 0.0;
    minus_1_power_j = -1.0;
    for(j=1;j<size_y;++j) { //get the neg part from the symmetry
      in[i*size_y_fft+size_y_fft-j][0] =  Real[i][j] * minus_1_power_j;
      in[i*size_y_fft+size_y_fft-j][1] = -Imag[i][j] * minus_1_power_j;
      minus_1_power_j *= -1.0;
    }
  }

  out = fftw_api->FFTAndGetOutput_Forward();

  for(i=0;i<size_x;++i) {
    for(j=0;j<size_y;++j) { //only need the first half
      Real[i][j] = out[i*size_y_fft+j][0];
      Imag[i][j] = out[i*size_y_fft+j][1];
    }
  }
}
*/

//New version with prunning
void FakedOptFFT::optimizedFFT(FFT_precision** & Real, FFT_precision** & Imag,
		const int size_x, const int size_y) {
	const int size_y_fft = size_y * 2;
	FFTW_API* fftw_api;
#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* in, *out;
#else
	fftw_complex* in, *out;
#endif

	const int nonZeroLines = (size_x + 1) / 2; //expect a padding of at least floor(size_x/2) at the bottom half

	int i, j;
	FFT_precision minus_1_power_j;

	//lines fft
	fftw_api = getFFTW_API(size_y_fft, 1);
	in = fftw_api->getIn_Forward();
	for (i = 0; i < nonZeroLines; ++i) {
		//each line looks like [x(0), .. , x(n), 0, x(-n), .. , x(-1)]
		//x(-i) = -1^i conj(x(i))
		for (j = 0; j < size_y; ++j) { //copy the first half
			in[j][0] = Real[i][j];
			in[j][1] = Imag[i][j];
		}
		in[size_y][0] = 0.0; //add on zero pad
		in[size_y][1] = 0.0;
		minus_1_power_j = -1.0;
		for (j = 1; j < size_y; ++j) { //get the neg part from the symmetry
			in[size_y_fft - j][0] = Real[i][j] * minus_1_power_j;
			in[size_y_fft - j][1] = -Imag[i][j] * minus_1_power_j;
			minus_1_power_j *= -1.0;
		}

		out = fftw_api->FFTAndGetOutput_Forward();

		//second half is useless
		for (j = 0; j < size_y; ++j) {
			Real[i][j] = out[j][0];
			Imag[i][j] = out[j][1];
		}

	}

	//column
	fftw_api = getFFTW_API(size_x, 1);
	in = fftw_api->getIn_Forward();
	for (j = 0; j < size_y; ++j) {
		for (i = 0; i < size_x; ++i) {
			in[i][0] = Real[i][j];
			in[i][1] = Imag[i][j];
		}

		out = fftw_api->FFTAndGetOutput_Forward();

		for (i = 0; i < size_x; ++i) {
			Real[i][j] = out[i][0];
			Imag[i][j] = out[i][1];
		}
	}

}

void FakedOptFFT::optimizedIFFT(FFT_precision** & Real, FFT_precision** & Imag,
		const int size_x, const int size_y) {
	FFTW_API* fftw_api;
#if defined(__SINGLE_PRECISION_FFT__)
	fftwf_complex* in, *out;
#else
	fftw_complex* in, *out;
#endif

	int i, j;
	const int size_y_fft = size_y * 2;

	//start by processing the columns
	fftw_api = getFFTW_API(size_x, 1);
	in = fftw_api->getIn_Backward();

	for (j = 0; j < size_y; ++j) {

		for (i = 0; i < size_x; ++i) {
			in[i][0] = Real[i][j];
			in[i][1] = Imag[i][j];
		}

		out = fftw_api->FFTAndGetOutput_Backward();

		for (i = 0; i < size_x; ++i) {
			Real[i][j] = out[i][0];
			Imag[i][j] = out[i][1];
		}
	}

	//process the rows
	//each row looks like [x(0), .. , x(n), conj(x(0)), .. , conj(x(n))]

	fftw_api = getFFTW_API(size_y * 2, 1);
	in = fftw_api->getIn_Backward();

	const FFT_precision scaling = 1.0 / (FFT_precision) (size_x * size_y_fft);

	for (i = 0; i < size_x; ++i) {
		for (j = 0; j < size_y; ++j) {
			in[j][0] = Real[i][j];
			in[j][1] = Imag[i][j];
			in[j + size_y][0] = Real[i][j];
			in[j + size_y][1] = -Imag[i][j];
		}

		out = fftw_api->FFTAndGetOutput_Backward();

		for (j = 0; j < size_y; ++j) { //only need the first half and rescale
			Real[i][j] = out[j][0] * scaling;
			Imag[i][j] = out[j][1] * scaling;
		}
	}
}
