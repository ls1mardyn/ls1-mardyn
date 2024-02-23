/*
 * fft_utils.cpp
 *
 *  Created on: Sep 14, 2015
 *      Author: gallardjm
 */
#include "fft_utils.h"

/*
 * Allocate a contiguous matrix
 * This way if a = &Mat[x][y];
 * a[z] <=> Mat[x+(y+z)/x][(y+z) mod x]
 * In particular a[x] = Mat[x+1][y] (easy column manipulation for UHFFT)
 */
FFT_precision** alloc_matrix(int n, int m) {
	FFT_precision **mat = new FFT_precision*[n];
	FFT_precision *ptrpool;
	if (FFTSettings::USE_VECTORIZATION)
		ptrpool = (FFT_precision*) _mm_malloc(n * m * sizeof(FFT_precision),
				__FFT_MATRIX_ALIGNMENT__);
	else
		ptrpool = new FFT_precision[n * m];
	for (int i = 0; i < n; i++, ptrpool += m)
		mat[i] = ptrpool;

	return mat;
}

/*
 * Set a matrix to 0
 */
void clear_matrix(FFT_precision** mat, int n, int m) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			mat[i][j] = 0.0;
}

/*
 * Delete a contiguous matrix alocated with alloc_matrix
 */
void delete_matrix(FFT_precision** mat) {
	if (mat == NULL)
		return;
	if (FFTSettings::USE_VECTORIZATION)
		_mm_free(mat[0]);
	else
		delete[] mat[0];
	delete[] mat;
}

/*
 * Print a matrix for debug
 */
void print_matrix(FFT_precision** mat, int n, int m) {
	std::cout << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			std::cout << mat[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

FFT_precision** copy_matrix(FFT_precision** mat, int n, int m) {
	FFT_precision** nMat = alloc_matrix(n, m);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			nMat[i][j] = mat[i][j];

	return nMat;
}

FFT_precision* copy_aligned_array(FFT_precision* a, int size) {
	FFT_precision* copy = alloc_aligned_array(size);
	for (int i = 0; i < size; ++i) {
		copy[i] = a[i];
	}

	return copy;
}

FFT_precision*** alloc_blocks(int nbBlocks, int nx, int ny) {
	FFT_precision*** blocks = new FFT_precision**[nbBlocks];
	for (int i = 0; i < nbBlocks; i++)
		blocks[i] = alloc_matrix(nx, ny);

	return blocks;
}

void clear_blocks(FFT_precision*** blocks, int nbBlocks, int nx, int ny) {
	for (int i = 0; i < nbBlocks; i++)
		clear_matrix(blocks[i], nx, ny);
}

FFT_precision*** copy_blocks(FFT_precision*** blocks, int nbBlocks, int nx,
		int ny) {
	FFT_precision*** copy = alloc_blocks(nbBlocks, nx, ny);
	for (int i = 0; i < nbBlocks; i++)
		for (int j = 0; j < nx; j++)
			for (int k = 0; k < ny; k++)
				copy[i][j][k] = blocks[i][j][k];

	return copy;
}

void delete_blocks(FFT_precision*** blocks, int nbBlocks) {
	for (int i = 0; i < nbBlocks; i++)
		delete_matrix(blocks[i]);

	delete[] blocks;
}

void print_block(FFT_precision*** blocks, int nbBlocks, int nx, int ny,
		bool ascending) {
	int b;
	if (ascending) {
		for (b = 0; b < nbBlocks; b++)
			print_matrix(blocks[b], nx, ny);
	} else {
		for (b = nbBlocks - 1; b >= 0; b--)
			print_matrix(blocks[b], nx, ny);
	}
}

FFT_precision*** alloc_scalBlocks(int nbBlocks, int nx, int* ny_sizes) {
	FFT_precision*** blocks = new FFT_precision**[nbBlocks];
	for (int i = 0; i < nbBlocks; i++)
		blocks[i] = alloc_matrix(nx, ny_sizes[i]);

	return blocks;
}

FFT_precision*** copy_scalBlocks(FFT_precision*** blocks, int nbBlocks, int nx,
		int* ny_sizes) {
	FFT_precision*** copy = alloc_scalBlocks(nbBlocks, nx, ny_sizes);
	for (int i = 0; i < nbBlocks; i++)
		for (int j = 0; j < nx; j++)
			for (int k = 0; k < ny_sizes[i]; k++)
				copy[i][j][k] = blocks[i][j][k];

	return copy;
}

//Block using arrays instead of matrices

FFT_precision** alloc_blocks_arr(int nbBlocks, int nx, int ny) {
	FFT_precision** blocks = new FFT_precision*[nbBlocks];
	for (int i = 0; i < nbBlocks; i++)
		blocks[i] = alloc_aligned_array(nx * ny);

	return blocks;
}

void clear_blocks_arr(FFT_precision** blocks, int nbBlocks, int nx, int ny) {
	for (int i = 0; i < nbBlocks; i++)
		clear_array(blocks[i], nx * ny);
}

FFT_precision** copy_blocks_arr(FFT_precision** blocks, int nbBlocks, int nx,
		int ny) {
	FFT_precision** copy = alloc_blocks_arr(nbBlocks, nx, ny);
	for (int i = 0; i < nbBlocks; i++)
		for (int j = 0; j < nx * ny; j++)
			copy[i][j] = blocks[i][j];

	return copy;
}

void delete_blocks_arr(FFT_precision** blocks, int nbBlocks) {
	for (int i = 0; i < nbBlocks; i++)
		delete_aligned_array(blocks[i]);

	delete[] blocks;
}

FFT_precision** alloc_scalBlocks_arr(int nbBlocks, int nx, int* ny_sizes) {
	FFT_precision** blocks = new FFT_precision*[nbBlocks];
	for (int i = 0; i < nbBlocks; i++)
		blocks[i] = alloc_aligned_array(nx * ny_sizes[i]);

	return blocks;
}

FFT_precision** copy_scalBlocks_arr(FFT_precision** blocks, int nbBlocks,
		int nx, int* ny_sizes) {
	FFT_precision** copy = alloc_scalBlocks_arr(nbBlocks, nx, ny_sizes);
	for (int i = 0; i < nbBlocks; i++)
		for (int j = 0; j < nx * ny_sizes[i]; j++)
			copy[i][j] = blocks[i][j];

	return copy;
}
