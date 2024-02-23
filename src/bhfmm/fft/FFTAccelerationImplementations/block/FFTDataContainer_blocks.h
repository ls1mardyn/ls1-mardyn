/*
 * FFTDataContainer_blocks.h
 *
 *  Created on: Feb 12, 2016
 *      Author: gallardjm
 */
#ifndef FFTDATA_BLOCK_H_
#define FFTDATA_BLOCK_H_

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
class FFTDataContainer_blocks: public FFTDataContainer {

public:
	FFTDataContainer_blocks(int nbBlocks, int nx, int ny) :
			_nbBlocks(nbBlocks), _nx(nx), _ny(ny) {
	}

	FFT_precision*** Re;
	FFT_precision*** Im;
	int _nbBlocks;
	int _nx;
	int _ny;

	~FFTDataContainer_blocks() {
		delete_blocks(Re, _nbBlocks);
		delete_blocks(Im, _nbBlocks);
	}

	FFTDataContainer* copyContainer() {
		FFTDataContainer_blocks* copy = new FFTDataContainer_blocks(_nbBlocks,
				_nx, _ny);
		copy->Re = copy_blocks(Re, _nbBlocks, _nx, _ny);
		copy->Im = copy_blocks(Im, _nbBlocks, _nx, _ny);

		return copy;
	}
};

#endif
