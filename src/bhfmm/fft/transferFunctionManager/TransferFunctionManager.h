/*
 * TransferFunctionManager.h
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#ifndef TRANSFERFUNCTIONMANAGER_H_
#define TRANSFERFUNCTIONMANAGER_H_

#include <math.h>
#include <iostream>
#include "bhfmm/fft/transferFunctionManager/DummyExpansion.h"
#include "bhfmm/fft/FFTDataContainer.h"
#include "bhfmm/fft/FFTAccelerationAPI.h"
#include "bhfmm/fft/TransferFunctionManagerAPI.h"

/**
 * Class managing the TransferFunction
 * To use in the main M2L loop to get the FFTDataContainer* transferfunction
 * required by the FFTAcceleration's M2L
 *
 * Delegate FFT conversions and handling of the true implementation of
 * the abstract FFTDataContainer to a FFTAcceleration
 */
class TransferFunctionManager: public TransferFunctionManagerAPI {

public:
	/**
	 * Constructor
	 *
	 * @param int ord, order of the expansions
	 * @param FFTAcceleration* FFTA, FFTAcceleration implementation used for the FMM
	 * @param bool verbose, output some stats during destruction
	 * @param bool v_child, as a verbose child that will be upcasted, shouldn't be verbose at this level
	 */
	TransferFunctionManager(int ord, FFTAccelerationAPI* FFTA, bool verbose);

	//! virtual destructor since inheriting classes will be upcasted
	virtual ~TransferFunctionManager();

	/**
	 * get a transfer function corresponding the the input vector
	 *
	 * @param utils::Vector<double,3> r, vector between the source and the target
	 * @return FFTDataContainer*, the transfer function
	 */
	FFTDataContainer* getTransferFunction(int x, int y, int z,
			double cell_size_x, double cell_size_y, double cell_size_z);

protected:
	int _ord; //order of the expansions
	bool _verbose; //verbose setting
	int _asked; //stat
	int _builded; //stat
	FFTAccelerationAPI* _FFTAcceleration; //FFTAcceleration to delegate FFT operations

};

#endif //TRANSFERFUNCTIONMANAGER_H_
