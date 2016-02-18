/*
 * SolidHarmonicsStorage.h
 *
 *  Created on: Nov 20, 2014
 *      Author: tchipevn, rajats
 */

#ifndef SOLIDHARMONICSSTORAGE_H_
#define SOLIDHARMONICSSTORAGE_H_

#include <cassert>

namespace bhfmm {
class SolidHarmonicsStorage;

/**
 * swap function (for the copy-and-swap idiom)
 * @param s1
 * @param s2
 */
void swap(SolidHarmonicsStorage & s1, SolidHarmonicsStorage & s2);

/**
 * Class to store the solid-harmonics expansions in the form of a triangular matrix.
 * Provides, storage, means to access it and trivial (entrywise) math operations.
 */
class SolidHarmonicsStorage {
public:
	/**
	 * constructor
	 * @param numRows - the number of rows that the (lower) triangular matrix should have
	 * @param initializeToZero - if true, values are initialized to zero, otherwise left uninitialized
	 */
	SolidHarmonicsStorage(int numRows = 0, bool initializeToZero = true);

	/**
	 * copy constructor
	 * @param s
	 */
	SolidHarmonicsStorage(const SolidHarmonicsStorage& s);

	/**
	 * Destructor - frees allocated memory
	 */
	~SolidHarmonicsStorage();

	/**
	 * swap function (for the copy-and-swap idiom)
	 * @param s1
	 * @param s2
	 */
	friend void bhfmm::swap(SolidHarmonicsStorage& s1, SolidHarmonicsStorage& s2);

	/**
	 * assignment operator: values are copied entrywise
	 * @param s
	 * @return
	 */
	SolidHarmonicsStorage& operator =(SolidHarmonicsStorage s);

	/**
	 * entrywise addition
	 * @param s
	 * @return
	 */
	SolidHarmonicsStorage& operator+=(const SolidHarmonicsStorage& s);

	/**
	 * scalar multiplication
	 * @param scalar
	 * @return
	 */
	SolidHarmonicsStorage& operator*=(double scalar);

	/**
	 * @return the number of rows
	 */
	int getNumRows() const;

	/**
	 * @return the total number of entries
	 */
	int getTotalNumValues() const;

	/**
	 * set all entries to zero
	 */
	void setToZero();

	/**
	 * calculate the sequential index, corresponding to entry (l,m)
	 * @param l row
	 * @param m column
	 * @return corresponding sequential index
	 */
	int index(int l, int m) const {
		assert(l >= 0);
		assert(m >= 0);
		assert(l <= _numRows);
		assert(m <= l);
		return l * (l + 1) / 2 + m;
	}

	/**
	 * access to value in matrix form
	 * @param l row
	 * @param m column
	 * @return entry at (l,m) by reference (for modification)
	 */
	double& getValue(int l, int m) {
		return _values[index(l, m)];
	}

	/**
	 * access to value sequentially
	 * @param i index
	 * @return entry at (i) by reference (for modification)
	 */
	double& getValueSequential(int i) {
		return _values[i];
	}

	/**
	 * const-access to value in matrix form
	 * @param l row
	 * @param m column
	 * @return copy of entry at (l,m) (for read-only access)
	 */
	double getValueConst(int l, int m) const {
		return _values[index(l, m)];
	}

	/**
	 * const-access to value sequentially
	 * @param i entry
	 * @return copy of entry at (i) (for read-only access)
	 */
	double getValueConstSequential(int i) const {
		return _values[i];
	}

private:
	/**
	 * the number of rows
	 */
	int _numRows;

	/**
	 * total number of entries
	 */
	int _totalNumValues;

	/**
	 * values
	 */
	double* _values;
};

/**
 * operator+ entrywise addition
 * @param lhs
 * @param rhs
 * @return
 */
SolidHarmonicsStorage operator+(SolidHarmonicsStorage lhs, const SolidHarmonicsStorage & rhs);

/**
 * operator* scalar multiplication
 * @param scalar
 * @param s
 * @return
 */
SolidHarmonicsStorage operator*(double scalar, SolidHarmonicsStorage s);

} /* namespace bhfmm */

#endif /* SOLIDHARMONICSSTORAGE_H_ */
