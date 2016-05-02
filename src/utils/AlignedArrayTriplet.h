/*
 * AlignedArrayTriplet.h
 *
 *  Created on: Mar 23, 2016
 *      Author: tchipevn
 */

#ifndef ALIGNEDARRAYTRIPLET_H
#define ALIGNEDARRAYTRIPLET_H

#include "AlignedArray.h"
#include <cassert>

class AlignedArrayTriplet {
public:
	AlignedArrayTriplet(size_t initialSize) : _data(), _numEntriesPerArray(0) {
		resize(initialSize);
	}

	double& x(size_t i = 0) {
		assert(i < _numEntriesPerArray);
		return _data[i];
	}

	double& x(size_t i = 0) const {
		assert(i < _numEntriesPerArray);
		return _data[i];
	}

	double& y(size_t i = 0) {
		assert(i < _numEntriesPerArray);
		return _data[i + _numEntriesPerArray];
	}

	double& y(size_t i = 0) const {
		assert(i < _numEntriesPerArray);
		return _data[i + _numEntriesPerArray];
	}

	double& z(size_t i = 0) {
		assert(i < _numEntriesPerArray);
		return _data[i + 2 * _numEntriesPerArray];
	}

	double& z(size_t i = 0) const {
		assert(i < _numEntriesPerArray);
		return _data[i + 2 * _numEntriesPerArray];
	}

	/**
	 * \brief Reallocate the array. All content may be lost.
	 */
	void resize(size_t nEntriesPerArray) {
		if (nEntriesPerArray == 0)
			return;

		// make sure that this is divisible by eight, otherwise we break alignment
		// get the rounded-up number of entries
		_numEntriesPerArray = _data._round_up(nEntriesPerArray);
		_data.resize(3 * _numEntriesPerArray);
		// set elements from nEntriesPerArray to _numEntriesPerArray to zero:
		size_t elements = _numEntriesPerArray - nEntriesPerArray;
		std::memset(&(x()), 0, elements * sizeof(double));
		std::memset(&(y()), 0, elements * sizeof(double));
		std::memset(&(z()), 0, elements * sizeof(double));
	}

	size_t get_dynamic_memory() const {
		return _data.get_dynamic_memory();
	}

private:
	AlignedArray<double> _data;

	//! at any point, _data contains precisely three times this storage
	size_t _numEntriesPerArray;
};

#endif /* ALIGNEDARRAYTRIPLET_H */
