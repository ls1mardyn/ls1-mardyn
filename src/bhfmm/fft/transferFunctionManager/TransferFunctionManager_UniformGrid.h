/*
 * TransferFunctionManager_UniformGrid.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#ifndef TRANSFERFUNCTIONMANAGER_UNI_H_
#define TRANSFERFUNCTIONMANAGER_UNI_H_

#include <math.h>
#include <iostream>
#include <stdlib.h>     /* abs */
#include "bhfmm/fft/transferFunctionManager/TransferFunctionManager.h"
#include "bhfmm/fft/TransferFunctionManagerAPI.h"
#include "bhfmm/fft/transferFunctionManager/DummyExpansion.h"
#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/FFTAccelerationAPI.h"

/**
 * Uniform grid variant of the generalist TransferFunctionManager class
 * Use precomputation since they are only ~300 different tf to be used
 * Rescaling is implied and required
 */
class TransferFunctionManager_UniformGrid: public TransferFunctionManagerAPI {

public:
	//! Constructor, create the storage and call buildTransferFunction() to set it
	TransferFunctionManager_UniformGrid(int ord, FFTAccelerationAPI* FFTA,
			bool verbose);
	//! destructor, clean the storage (free all memory used)
	~TransferFunctionManager_UniformGrid();

	/**
	 * get the transfer function, x,y,z are cell coordinate differences
	 * in a 3d cell storage, -3 <= x,y,z <= 3
	 * assumes cell_size_x, cell_size_y, cell_size_z = 2/sqrt(3)
	 */
	FFTDataContainer* getTransferFunction(int x, int y, int z,
			double cell_size_x, double cell_size_y, double cell_size_z);

private:
	/**
	 * Precompute the transfer functions and store them in the _storage
	 */
	void buildTransferFunction();
	FFTDataContainer**** _storage; //! 3d array storage for the precomputed tf

	int _ord; //order of the expansions
	bool _verbose; //verbose setting
	int _asked; //stat
	int _builded; //stat
	FFTAccelerationAPI* _FFTAcceleration; //FFTAcceleration to delegate FFT operations
};

#endif //TRANSFERFUNCTIONMANAGER_UNI_H_
