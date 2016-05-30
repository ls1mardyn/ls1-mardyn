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

	double* xBegin() {
		return _numEntriesPerArray > 0 ? &(_data[0 * _numEntriesPerArray]) : 0;
	}

	double* xBegin() const {
		return _numEntriesPerArray > 0 ? &(_data[0 * _numEntriesPerArray]) : 0;
	}

	double* yBegin() {
		return _numEntriesPerArray > 0 ? &(_data[1 * _numEntriesPerArray]) : 0;
	}

	double* yBegin() const {
		return _numEntriesPerArray > 0 ? &(_data[1 * _numEntriesPerArray]) : 0;
	}

	double* zBegin() {
		return _numEntriesPerArray > 0 ? &(_data[2 * _numEntriesPerArray]) : 0;
	}

	double* zBegin() const {
		return _numEntriesPerArray > 0 ? &(_data[2 * _numEntriesPerArray]) : 0;
	}

	double& x(size_t i) {
		assert(i < _numEntriesPerArray);
		return _data[i];
	}

	double& x(size_t i) const {
		assert(i < _numEntriesPerArray);
		return _data[i];
	}

	double& y(size_t i) {
		assert(i < _numEntriesPerArray);
		return _data[i + _numEntriesPerArray];
	}

	double& y(size_t i) const {
		assert(i < _numEntriesPerArray);
		return _data[i + _numEntriesPerArray];
	}

	double& z(size_t i) {
		assert(i < _numEntriesPerArray);
		return _data[i + 2 * _numEntriesPerArray];
	}

	double& z(size_t i) const {
		assert(i < _numEntriesPerArray);
		return _data[i + 2 * _numEntriesPerArray];
	}

	/**
	 * \brief Reallocate the array. All content may be lost.
	 */
	void resize(size_t nEntriesPerArray) {
		if (nEntriesPerArray == 0 and _numEntriesPerArray == 0)
			return;

		// make sure that this is divisible by eight, otherwise we break alignment
		// get the rounded-up number of entries
		_numEntriesPerArray = _data._round_up(nEntriesPerArray);
		_data.resize(3 * _numEntriesPerArray);

		if (_numEntriesPerArray > nEntriesPerArray) {
			// set elements from nEntriesPerArray to _numEntriesPerArray to zero:
			size_t elements = _numEntriesPerArray - nEntriesPerArray;
			std::memset(&(x(nEntriesPerArray)), 0, elements * sizeof(double));
			std::memset(&(y(nEntriesPerArray)), 0, elements * sizeof(double));
			std::memset(&(z(nEntriesPerArray)), 0, elements * sizeof(double));
		}
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
