/*
 * FFTDataContainer_blocks.h
 *
 *  Created on: Mar 15, 2016
 *      Author: gallardjm
 */
#ifndef FFTDATA_SCALBLOCK_V0_H_
#define FFTDATA_SCALBLOCK_V0_H_

#include "bhfmm/fft/FFTSettings_preprocessor.h" //tmp include for the typedef FFT_precision
#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/tools/fft_utils.h"
#include <stdlib.h>

/*
 * Storage for blocks, can be used with scaling Blocks (scalBlock)
 *
 * blocks[i] is a matrix reprensenting the ith block
 * Convention: block 0 contain the terms of lowest order
 *  => source matrix in ascending block (block 0 is the top of the matrix),
 *      tf and target matrices in descending blocks (block 0 is the lower part of the matrix)
 */
class FFTDataContainer_scalBlocks_v0: public FFTDataContainer {

public:
	FFTDataContainer_scalBlocks_v0(int nbBlocks, int nx, int ny, bool scal,
			int* blockSize) :
			_nbBlocks(nbBlocks), _nx(nx), _ny(ny), _scal(scal), _blockSize(
					blockSize) {
	}

	FFT_precision*** Re;
	FFT_precision*** Im;
	int _nbBlocks;
	int _nx;
	int _ny;
	bool _scal;
	int* _blockSize; //direct copy of the pointer, do not modify or free

	~FFTDataContainer_scalBlocks_v0() {
		delete_blocks(Re, _nbBlocks);
		delete_blocks(Im, _nbBlocks);
	}

	FFTDataContainer* copyContainer() {
		FFTDataContainer_scalBlocks_v0* copy =
				new FFTDataContainer_scalBlocks_v0(_nbBlocks, _nx, _ny, _scal,
						_blockSize);
		if (_scal) {
			copy->Re = copy_scalBlocks(Re, _nbBlocks, _nx, _blockSize);
			copy->Im = copy_scalBlocks(Im, _nbBlocks, _nx, _blockSize);
		} else {
			copy->Re = copy_blocks(Re, _nbBlocks, _nx, _ny);
			copy->Im = copy_blocks(Im, _nbBlocks, _nx, _ny);
		}

		return copy;
	}
};

#endif
