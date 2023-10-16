#ifndef FFTSETTINGS_PP_H_
#define FFTSETTINGS_PP_H_

//Preprocessor instruction and option to choose at compile time

//!single or double precision
//#define __SINGLE_PRECISION_FFT__ //use float instead of double, if using fftw link -lfftw3f in makefile
#if defined(__SINGLE_PRECISION_FFT__)
typedef float FFT_precision;
#else
typedef double FFT_precision;
#endif

//#define __TEST_FAKE_VECTORIZATION__ //To fake vectorization on a computer not supporting it (replace _mm_malloc with malloc and _mm_free with free)
//!Vectorization settings are defined at compile time
//Intel SSE = 16, AVX = 32, AVX-512 = 64 (https://software.intel.com/en-us/articles/fdiag15126)
#define __FFT_MATRIX_ALIGNMENT__  32

#endif
