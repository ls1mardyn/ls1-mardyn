/*
 * FakeOptFFT.h
 *
 *  Created on: Mar 06, 2015
 *      Author: gallardjm
 */

#ifndef FAKEOPTFFT_H_
#define FAKEOPTFFT_H_

#include "WrapOpenMP.h"

#include "bhfmm/fft/tools/optimizedFFT/optFFT_API.h"
#include "bhfmm/fft/tools/FFTW_API.h"
#include "bhfmm/fft/tools/fft_utils.h"

#include <map>
#include <iostream>

//struct storing a int[2] for the map
struct pos {
	int x, y;
};

//comparator of pos (=int[2])
struct pos_comp {
	bool operator()(const pos &l, const pos &r) const {
		return (l.x < r.x || (l.x == r.x && l.y < r.y));
	}
};

/**
 * Implementation of the scheme of the optimized FFT that reconstruct the whole matrix
 * and use FFTW to perform the FFT/IFFT
 * Use a map<pos,FFTW_API*> to store the various FFTW_API needed to perform the optimized FFT
 */
class FakedOptFFT: public optFFT_API {

public:

	FakedOptFFT() {
	} //cout << "Using faked opt FFT (see bhfmm/fft/tools/optimizedFFT)" <<  endl;}
	~FakedOptFFT() { //free all entry of the map and the map
		auto itr = _fftw_api_map.begin();
		while (itr != _fftw_api_map.end()) {
			for (int i = 0; i < mardyn_get_max_threads(); ++i) {
				delete ((*itr).second[i]);
			}
			delete ((*itr).second);
			_fftw_api_map.erase(itr++);
		}
		_fftw_api_map.clear();
	}

	void optimizedFFT(FFT_precision** & Real, FFT_precision** & Imag,
			const int size_x, const int size_y);
	void optimizedIFFT(FFT_precision** & Real, FFT_precision** & Imag,
			const int size_x, const int size_y);

private:
	std::map<pos, FFTW_API**, pos_comp> _fftw_api_map; //storage of the various FFTW_API required

	FFTW_API* getFFTW_API(const int size_x, const int size_y); //memoized function using the map
};

#endif
