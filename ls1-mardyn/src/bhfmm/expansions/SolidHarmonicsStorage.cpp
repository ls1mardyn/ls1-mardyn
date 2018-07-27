/*
 * SolidHarmonicsStorage.cpp
 *
 *  Created on: Nov 20, 2014
 *      Author: tchipevn, rajats
 */

#include "SolidHarmonicsStorage.h"
#include <algorithm>

namespace bhfmm {

// CONSTRUCTORS //
SolidHarmonicsStorage::SolidHarmonicsStorage(int numRows, bool initializeToZero) :
		_numRows(numRows), _totalNumValues(numRows * (numRows + 1) / 2), _values(0) {
	mardyn_assert(numRows > 0);

	_values = new double[_totalNumValues];

	if (initializeToZero) {
		setToZero();
	}
}

SolidHarmonicsStorage::SolidHarmonicsStorage(const SolidHarmonicsStorage& s) :
		_numRows(s._numRows), _totalNumValues(s._totalNumValues), _values(0) {
	_values = new double[_totalNumValues];

	std::copy(s._values, s._values + s._totalNumValues, _values);
}

// DESTRUCTOR //
SolidHarmonicsStorage::~SolidHarmonicsStorage() {
	delete[] _values;
	_values = 0;
	_numRows = 0;
	_totalNumValues = 0;
}

void swap(SolidHarmonicsStorage& s1, SolidHarmonicsStorage& s2) {
	mardyn_assert(s1.getTotalNumValues() == s2.getTotalNumValues());
	std::swap(s1._values, s2._values);
}

// OPERATORS //
SolidHarmonicsStorage& SolidHarmonicsStorage::operator=(SolidHarmonicsStorage rhs) {
	mardyn_assert(this->_totalNumValues == rhs._totalNumValues);
	swap(*this, rhs);
	return *this;
}

SolidHarmonicsStorage& SolidHarmonicsStorage::operator+=(const SolidHarmonicsStorage& rhs) {
	mardyn_assert(this->_totalNumValues == rhs._totalNumValues);
	for (int i = 0; i < _totalNumValues; ++i) {
		this->_values[i] += rhs._values[i];
	}
	return *this;
}

inline SolidHarmonicsStorage operator+(SolidHarmonicsStorage lhs, const SolidHarmonicsStorage& rhs) {
	lhs += rhs;
	return lhs;
}

SolidHarmonicsStorage& SolidHarmonicsStorage::operator*=(double scalar) {
	for (int i = 0; i < _totalNumValues; ++i) {
		this->_values[i] *= scalar;
	}
	return *this;
}

inline SolidHarmonicsStorage operator*(double scalar, SolidHarmonicsStorage rhs) {
	rhs *= scalar;
	return rhs;
}

// METHODS //
inline int SolidHarmonicsStorage::getNumRows() const {
	return _numRows;
}

int SolidHarmonicsStorage::getTotalNumValues() const {
	return _totalNumValues;
}

void SolidHarmonicsStorage::setToZero() {
	std::fill(_values, _values + _totalNumValues, 0.0);
}

} /* namespace bhfmm */
