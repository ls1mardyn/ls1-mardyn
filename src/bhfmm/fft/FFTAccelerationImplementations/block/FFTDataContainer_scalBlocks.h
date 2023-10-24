/*
 * FFTDataContainer_scalBlocks.h
 *
 *  Created on: Mar 16, 2016
 *      Author: gallardjm
 */
#ifndef FFTDATA_SCALBLOCK_H_
#define FFTDATA_SCALBLOCK_H_

#include "bhfmm/fft/FFTSettings_preprocessor.h" //tmp include for the typedef FFT_precision
#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/tools/fft_utils.h"
#include <stdlib.h>

/*
 * Storage for blocks, can be used with scaling Blocks (scalBlock)
 *
 * blocks[i] is an array reprensenting the ith block
 * Convention: block 0 contain the terms of lowest order
 *  => source matrix in ascending block (block 0 is the top of the matrix),
 *      tf and target matrices in descending blocks (block 0 is the lower part of the matrix)
 */
class FFTDataContainer_scalBlocks: public FFTDataContainer {

public:
	FFTDataContainer_scalBlocks(int nbBlocks, int nx, int ny, bool scal,
			int* blockSize) :
			_nbBlocks(nbBlocks), _nx(nx), _ny(ny), _scal(scal), _blockSize(
					blockSize) {
	}

	FFT_precision** Re;
	FFT_precision** Im;
	int _nbBlocks;
	int _nx;
	int _ny;
	bool _scal;
	int* _blockSize;

	~FFTDataContainer_scalBlocks() {
		delete_blocks_arr(Re, _nbBlocks);
		delete_blocks_arr(Im, _nbBlocks);
	}

	FFTDataContainer* copyContainer() {
		FFTDataContainer_scalBlocks* copy = new FFTDataContainer_scalBlocks(
				_nbBlocks, _nx, _ny, _scal, _blockSize);
		if (_scal) {
			copy->Re = copy_scalBlocks_arr(Re, _nbBlocks, _nx, _blockSize);
			copy->Im = copy_scalBlocks_arr(Im, _nbBlocks, _nx, _blockSize);
		} else {
			copy->Re = copy_blocks_arr(Re, _nbBlocks, _nx, _ny);
			copy->Im = copy_blocks_arr(Im, _nbBlocks, _nx, _ny);
		}

		return copy;
	}
};

#endif
