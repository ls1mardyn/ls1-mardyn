/*
 * TransferFunctionManager_UniformGrid.cpp
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#include "TransferFunctionManager_UniformGrid.h"

//! Constructor, set the storage
TransferFunctionManager_UniformGrid::TransferFunctionManager_UniformGrid(
		int ord, FFTAccelerationAPI* FFTA, bool verbose) :
		_ord(ord), _verbose(verbose), _asked(0), _builded(0), _FFTAcceleration(
				FFTA) {

	_storage = new FFTDataContainer***[7];
	for (int i = 0; i < 7; ++i) {
		_storage[i] = new FFTDataContainer**[7];
		for (int j = 0; j < 7; ++j) {
			_storage[i][j] = new FFTDataContainer*[7];
			for (int k = 0; k < 7; ++k)
				_storage[i][j][k] = NULL;
		}
	}

	buildTransferFunction();
}

//! destructor, clean the storage (free all memory used)
TransferFunctionManager_UniformGrid::~TransferFunctionManager_UniformGrid() {
	if (_verbose)
		std::cout << "TransferFunctionManager_UniformGrid's stats: TF asked="
				<< _asked << " | TF builded=" << _builded << std::endl;

	for (int i = 0; i < 7; ++i) {
		for (int j = 0; j < 7; ++j) {
			for (int k = 0; k < 7; ++k) {
				if (_storage[i][j][k] != NULL)
					delete _storage[i][j][k];
			}
			delete[] _storage[i][j];
		}
		delete[] _storage[i];
	}
	delete[] _storage;
}

void TransferFunctionManager_UniformGrid::buildTransferFunction() {

	double base_unit = 2.0 / sqrt(3);
	double x, y, z;
	int i, j, k;
	TransferFunctionManager tfman(_ord, _FFTAcceleration, false);

	for (i = -3; i < 4; i++) {
		x = ((double) i) * base_unit;
		for (j = -3; j < 4; j++) {
			y = ((double) j) * base_unit;
			for (k = -3; k < 4; k++) {
				if (std::max(abs(i), std::max(abs(j), abs(k))) < 2)
					continue; //skip useless tf (cells too close)

				z = ((double) k) * base_unit;
				_storage[i + 3][j + 3][k + 3] = tfman.getTransferFunction(x, y,
						z, base_unit, base_unit, base_unit); //Use general class
				_builded++;
			}
		}
	}
}

FFTDataContainer* TransferFunctionManager_UniformGrid::getTransferFunction(
		int x, int y, int z, double cell_size_x, double cell_size_y,
		double cell_size_z) {
	_asked++;
	return _storage[x + 3][y + 3][z + 3];
}

