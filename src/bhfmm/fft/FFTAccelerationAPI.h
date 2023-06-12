/*
 * FFTAcceleration.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */
#ifndef FFTACC_H_
#define FFTACC_H_

#include <stdexcept>

#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/FFTAccelerableExpansion.h"

/**
 * Abstract class defining the interface of a class handling the M2L's
 * FFT acceleration
 *
 * Various implementations will offer various schemes and can be choosed
 * at runtime using a factory and strategy pattern (see FFTSetting and FFTFactory)
 */
class FFTAccelerationAPI {

public:

	//! destructor, child class will be upcasted, virtual destructor required to call the right child class destructor
	virtual ~FFTAccelerationAPI() {
	}

	//Initialization, to be done before the M2L

	/**
	 * Initialize a source expansion's FFTDatacontainer,
	 * including rescaling using the radius parameter (see doc/Rescalling)
	 *
	 * @param FFTAccelerableExpansion & Expansion
	 * @param double radius
	 */
	virtual void FFT_initialize_Source(FFTAccelerableExpansion & Expansion,
			double radius) =0;

	/**
	 * Initialize a Target expansion's FFTDatacontainer
	 * (schould be to a full 0 FFTDatacontainer)
	 *
	 * @param FFTAccelerableExpansion & Expansion
	 */
	virtual void FFT_initialize_Target(FFTAccelerableExpansion & Expansion) =0;

	/**
	 * Initialize a source expansion's FFTDatacontainer,
	 * all subsequent transformations (flip of the matrix or rescaling of
	 * scaling blocks for example) should be done here.
	 *
	 * @param FFTAccelerableExpansion & Expansion
	 */
	virtual void FFT_initialize_TransferFunction(
			FFTAccelerableExpansion & Expansion) =0;

	//M2L operator, with and without vectorization
	/**
	 * M2L operator, performs the M2L in fourier space, usually an entrywise
	 * product (or more complexe scheme if needed, see blocks implementations)
	 *
	 * All parameters needs to be initialized using the previous methods
	 * The TransferFunction should come from a TransferFunctionManager
	 *
	 * @param FFTAccelerableExpansion & Source
	 * @param FFTAccelerableExpansion & Target
	 * @param FFTDataContainer* TransferFunction
	 */
	virtual void FFT_M2L(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction) =0;

	/**
	 * M2L operator, performs the M2L in fourier space, usually an entrywise
	 * product (or more complexe scheme if needed, see blocks implementations)
	 *
	 * Uses vectorization of the loops
	 *
	 * All parameters needs to be initialized using the previous methods
	 * The TransferFunction should come from a TransferFunctionManager
	 *
	 * @param FFTAccelerableExpansion & Source
	 * @param FFTAccelerableExpansion & Target
	 * @param FFTDataContainer* TransferFunction
	 */
	virtual void FFT_M2L_vec(FFTAccelerableExpansion & Source,
			FFTAccelerableExpansion & Target,
			FFTDataContainer* TransferFunction) =0;

	//Finalize method, to use after the M2L operations on each target
	/**
	 * Finalize a Target expansion by adding the result of the M2L stored
	 * in Fourier space to the expansion values.
	 *
	 * Include the rescaling using the radius parameter, should be the same
	 * as the one used by FFT_initialize_Source
	 *
	 * @param FFTAccelerableExpansion & Expansion
	 * @param double radius
	 */
	virtual void FFT_finalize_Target(FFTAccelerableExpansion & Expansion,
			double radius) =0;

protected:
	int _p;      //! order of the expansions (start at 0 so usually order+1)
	int _fft_nx; //! number of line of the FFT matrices (or similar)
	int _fft_ny; //! number of column of the FFT matrices (or similar)
};

#endif
