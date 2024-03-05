/*
 * FFTW_Helper.cpp
 *
 *  Created on: Sep 14, 2015
 *      Author: gallardjm
 */

#include "FFTW_Helper.h"

#if defined(__SINGLE_PRECISION_FFT__)

FFTW_Helper::FFTW_Helper(const int p, const int nx, const int ny) : _p(p),_nx(nx), _ny(ny) {
	S2FFT_in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
	S2FFT_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
	T2FFT_in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
	T2FFT_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
	FFT2L_in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
	FFT2L_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*ny);

	//TODO
	//Fill all in array with dummy value

	S2FFT = fftwf_plan_dft_2d(nx, ny, S2FFT_in, S2FFT_out, FFTW_FORWARD, FFTW_ESTIMATE );// could use FFTW_MEASURE
	T2FFT = fftwf_plan_dft_2d(nx, ny, T2FFT_in, T2FFT_out, FFTW_FORWARD, FFTW_ESTIMATE );
	FFT2L = fftwf_plan_dft_2d(nx, ny, FFT2L_in, FFT2L_out, FFTW_BACKWARD, FFTW_ESTIMATE );
}

FFTW_Helper::~FFTW_Helper() {
	fftwf_destroy_plan(S2FFT);
	fftwf_destroy_plan(T2FFT);
	fftwf_destroy_plan(FFT2L);
	fftwf_free(S2FFT_in);
	fftwf_free(S2FFT_out);
	fftwf_free(T2FFT_in);
	fftwf_free(T2FFT_out);
	fftwf_free(FFT2L_in);
	fftwf_free(FFT2L_out);

	fftw_cleanup();
}

#else

FFTW_Helper::FFTW_Helper(const int p, const int nx, const int ny) :
		_p(p), _nx(nx), _ny(ny) {
	S2FFT_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	S2FFT_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	T2FFT_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	T2FFT_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	FFT2L_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	FFT2L_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);

	//TODO
	//Fill all in array with dummy value

	S2FFT = fftw_plan_dft_2d(nx, ny, S2FFT_in, S2FFT_out, FFTW_FORWARD, FFTW_ESTIMATE); // could use FFTW_MEASURE
	T2FFT = fftw_plan_dft_2d(nx, ny, T2FFT_in, T2FFT_out, FFTW_FORWARD, FFTW_ESTIMATE);
	FFT2L = fftw_plan_dft_2d(nx, ny, FFT2L_in, FFT2L_out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

FFTW_Helper::~FFTW_Helper() {
	fftw_destroy_plan(S2FFT);
	fftw_destroy_plan(T2FFT);
	fftw_destroy_plan(FFT2L);
	fftw_free(S2FFT_in);
	fftw_free(S2FFT_out);
	fftw_free(T2FFT_in);
	fftw_free(T2FFT_out);
	fftw_free(FFT2L_in);
	fftw_free(FFT2L_out);

	fftw_cleanup();
}

#endif

#if defined(__SINGLE_PRECISION_FFT__)
fftwf_complex* FFTW_Helper::Source2FFT(FFTAccelerableExpansion & Expansion, double radius)
#else
fftw_complex* FFTW_Helper::Source2FFT(FFTAccelerableExpansion & Expansion,
		double radius)
#endif
		{
	int n, m, fm;
	FFT_precision minus_one_power_n = 1.0;
	FFT_precision minus_one_neg_m = 1.0;
	FFT_precision r = (FFT_precision) radius;

	for (n = 0; n < _nx; n++) {
		for (m = 0; m < _ny; m++) {
			fm = (m + 1 + (_ny / 2)) % _ny;
			//fm = m;
			if (n < _p && _p - 1 <= m + n && _p - 1 >= m - n) {
				if (m - _p + 1 < 0) {
					minus_one_neg_m = 1.0 - 2.0 * (FFT_precision) ((_p - 1 - m) % 2);
					S2FFT_in[n * _ny + fm][0] = (FFT_precision) Expansion.get_C(n, _p - 1 - m) * minus_one_neg_m;
					S2FFT_in[n * _ny + fm][1] = (FFT_precision) Expansion.get_S(n, _p - 1 - m) * -minus_one_neg_m;
				} else {
					S2FFT_in[n * _ny + fm][0] = (FFT_precision) Expansion.get_C(n, m - _p + 1);
					S2FFT_in[n * _ny + fm][1] = (FFT_precision) Expansion.get_S(n, m - _p + 1);
				}

				S2FFT_in[n * _ny + fm][0] = minus_one_power_n * S2FFT_in[n * _ny + fm][0] / r;
				S2FFT_in[n * _ny + fm][1] = -minus_one_power_n * S2FFT_in[n * _ny + fm][1] / r; //we want to use use conj(Source)
			} else {
				S2FFT_in[n * _ny + fm][0] = 0.0;
				S2FFT_in[n * _ny + fm][1] = 0.0;
			}
		}
		minus_one_power_n *= -1.0;
		r *= (FFT_precision) radius;
	}

	execute_S2FFT();

	return S2FFT_out;
}

#if defined(__SINGLE_PRECISION_FFT__)
fftwf_complex* FFTW_Helper::TransferFunction2FFT(FFTAccelerableExpansion & Expansion)
#else
fftw_complex* FFTW_Helper::TransferFunction2FFT(
		FFTAccelerableExpansion & Expansion)
#endif
		{
	int n, m, fm;
	FFT_precision minus_one_neg_m = 1.0;

	for (n = 0; n < _nx; n++)
		for (m = 0; m < _ny; m++) {
			fm = (m + (_ny / 2)) % _ny;
			//fm = m;
			if (n < _p && _p <= m + n && _p >= m - n) {
				if (_p - 1 - m < 0) {
					minus_one_neg_m = 1.0 - 2.0 * (FFT_precision) ((m - _p) % 2);
					T2FFT_in[(_p - n - 1) * _ny + fm][0] = (FFT_precision) Expansion.get_C(n, m - _p) * minus_one_neg_m;
					T2FFT_in[(_p - n - 1) * _ny + fm][1] = (FFT_precision) Expansion.get_S(n, m - _p) * -minus_one_neg_m;
				} else {
					T2FFT_in[(_p - n - 1) * _ny + fm][0] = (FFT_precision) Expansion.get_C(n, _p - m);
					T2FFT_in[(_p - n - 1) * _ny + fm][1] = (FFT_precision) Expansion.get_S(n, _p - m);
				}

			} else {
				T2FFT_in[n * _ny + fm][0] = 0.0;
				T2FFT_in[n * _ny + fm][1] = 0.0;
			}
		}

	execute_T2FFT();

	return T2FFT_out;
}

void FFTW_Helper::FFT2Local(FFTAccelerableExpansion & Expansion,
		double radius) {
	execute_FFT2L();

	int n, m;
	double minus_1_power_m;
	FFT_precision scaling = 1.0 / (FFT_precision) (_nx * _ny);

	for (n = 0; n < _p; n++) {
		minus_1_power_m = 1.0;
		for (m = 0; m <= n; m++) {
			//Expansion.get_C(n,m) += (double)(FFT2L_out[(_p-n-1)*_ny+(_ny-m-1)][0] * scaling); //if fm=m
			//Expansion.get_S(n,m) += (double)(FFT2L_out[(_p-n-1)*_ny+(_ny-m-1)][1] * scaling);
			Expansion.get_C(n, m) += (double) (FFT2L_out[(_p - n - 1) * _ny + m][0] * scaling) * minus_1_power_m;
			Expansion.get_S(n, m) += (double) (FFT2L_out[(_p - n - 1) * _ny + m][1] * scaling) * -minus_1_power_m;
			minus_1_power_m *= -1.0;
		}
		scaling /= (FFT_precision) radius;
	}
}

