/*
 * fft_utils.h
 *
 *  Created on: Sep 14, 2015
 *      Author: gallardjm
 */

#ifndef FFT_UTILS_
#define FFT_UTILS_

#include <stdio.h>
#include <string.h> //memset
#include <iostream>
#include "bhfmm/fft/FFTSettings_preprocessor.h"
#include "bhfmm/fft/FFTSettings.h"


#if defined(__TEST_FAKE_VECTORIZATION__)
  #include <stdlib.h>
  #define _mm_malloc(a,b) malloc(a)
  #define _mm_free free
#else
  #include <mm_malloc.h>
#endif

/**
 * Allocate a contiguous matrix of double of size n*m
 * Contiguous <=> (a = &Mat[x][y]; => (a[z] <=> Mat[x+(y+z)/x][(y+z) mod x]))
 * Will create an aligned Matrix if specified in FFT_Settings.h (see vectorization)
 * @param n
 * @param m
 * @return double** matrix
 */
FFT_precision** alloc_matrix(int n, int m);


/**
 * clear a matrix (mat) of double of size n*m
 * @param mat
 * @param n
 * @param m
 */
void clear_matrix(FFT_precision** mat, int n, int m);


/**
 * Delete a matrix allocated using alloc_matrix
 * @param mat
 */
void delete_matrix(FFT_precision** mat);


/**
 * Deep copy a matrix mat of size n*m
 * @param mat
 * @param n
 * @param m
 * @return double** mat_copy
 */
FFT_precision** copy_matrix(FFT_precision** mat, int n, int m);


/**
 * Print a matrix, used for debug
 */
void print_matrix(FFT_precision** mat, int n, int m);

/**
 * Compute the number of zero necessary to pad a 2way array and preserve
 * data alignement (see FFTAcceretionImplementation/_2way_)
 */
inline int nb_zeroes_2wayFFT(int n, int m)
{
  int nb_double_aligned = __FFT_MATRIX_ALIGNMENT__ / sizeof(FFT_precision);
  
  return (nb_double_aligned - ((n*m/2) % nb_double_aligned)) % nb_double_aligned;
}

/**
 * Allocate an aligned array using _mm_malloc
 */
inline FFT_precision* alloc_aligned_array(int size)
{
  return (FFT_precision*)_mm_malloc(size * sizeof(FFT_precision), __FFT_MATRIX_ALIGNMENT__ );
}

/**
 * Set an array to 0
 */
inline void clear_array(FFT_precision* a, int size)
{
  for(int i=0;i<size;i++)
    a[i] = 0.0;
}

/**
 * deep copy an array
 */
FFT_precision* copy_aligned_array(FFT_precision* a, int size);

/**
 * free an aligned array with _mm_free
 */
inline void delete_aligned_array(FFT_precision* a)
{
  _mm_free(a);
}


//blocks scheme, same as matrix one
FFT_precision*** alloc_blocks(int nbBlocks, int nx, int ny);

void clear_blocks(FFT_precision*** blocks, int nbBlocks, int nx, int ny);

FFT_precision*** copy_blocks(FFT_precision*** blocks, int nbBlocks, int nx, int ny);

void delete_blocks(FFT_precision*** blocks, int nbBlocks);

void print_block(FFT_precision*** blocks, int nbBlocks, int nx, int ny, bool ascending);

//used by scalBlock_v0
FFT_precision*** alloc_scalBlocks(int nbBlocks, int nx, int* ny_sizes);

FFT_precision*** copy_scalBlocks(FFT_precision*** blocks, int nbBlocks, int nx, int* ny_sizes);



//array blocks for scalBlock
FFT_precision** alloc_blocks_arr(int nbBlocks, int nx, int ny);

void clear_blocks_arr(FFT_precision** blocks, int nbBlocks, int nx, int ny);

FFT_precision** copy_blocks_arr(FFT_precision** blocks, int nbBlocks, int nx, int ny);

void delete_blocks_arr(FFT_precision** blocks, int nbBlocks);

FFT_precision** alloc_scalBlocks_arr(int nbBlocks, int nx, int* ny_sizes);

FFT_precision** copy_scalBlocks_arr(FFT_precision** blocks, int nbBlocks, int nx, int* ny_sizes);



#endif /* FFT_UTILS_ */
