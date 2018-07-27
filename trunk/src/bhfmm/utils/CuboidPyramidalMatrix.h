/*
 * CuboidPyramidalMatrix.h
 *
 *  Created on: Dec 11, 2014
 *      Author: uwe
 */

#ifndef CUBOIDPYRAMIDALMATRIX_H_
#define CUBOIDPYRAMIDALMATRIX_H_

#include "utils/mardyn_assert.h"
#include <cstdlib>

namespace bhfmm {

class CuboidPyramidalMatrix {

public:
	CuboidPyramidalMatrix() : _totalNumEntries(0), _numSlices(0), _entries(NULL) {}

	CuboidPyramidalMatrix(unsigned N, bool clean = false) :
			_totalNumEntries(N*(N+1)*(4*N-1)/6), _numSlices(N) {
		if (N != 0) {
			_entries = new double[_totalNumEntries];
		} else {
			_entries = NULL;
		}
		if (clean) clear();
	}

	CuboidPyramidalMatrix(const CuboidPyramidalMatrix& c) :
			_totalNumEntries(c._totalNumEntries), _numSlices(c._numSlices) {

		_entries = new double[_totalNumEntries];
		std::copy(c._entries, c._entries + c._totalNumEntries, _entries);
	}

	/* Assignment using copy-and-swap idiom */
	CuboidPyramidalMatrix& operator=(CuboidPyramidalMatrix rhs) {
		std::swap(_totalNumEntries, rhs._totalNumEntries);
		std::swap(_numSlices, rhs._numSlices);
		std::swap(_entries, rhs._entries);
		return *this;
	}

	~CuboidPyramidalMatrix() {if (_entries != NULL) delete[] _entries;}

	void clear() { for (unsigned i = 0; i < _totalNumEntries; ++i) _entries[i] = 0.0; }

	unsigned get_num_entries() const {return _totalNumEntries;}

    inline double & access(unsigned i, unsigned j, int k) {return _entries[index(i,j,k)];}
    inline const double& access_const(unsigned i, unsigned j, int k) const {return _entries[index(i,j,k)];}
    inline double & access_seq(unsigned i) {return _entries[i];}
    inline double access_seq_const(unsigned i) const {return _entries[i];}

    inline unsigned index(unsigned i, unsigned j, int k) const {
      mardyn_assert(i <= _numSlices);
      mardyn_assert(j <= i);
      mardyn_assert(static_cast<unsigned>(abs(k))<= i);
      return i*(i+1)*(4*i-1)/6 + (2*i + 1)*j + i + k;  // valid for i < 1024

      /* i*(i+1)*(4*i-1)/6	: select slice.
       * + (2*i + 1)*j		: select row
       * + i + k  			: select column. k is signed -> start at the middle column i */
    }

private:

	unsigned _totalNumEntries;
	unsigned _numSlices;

	double* _entries;
};

} /* namespace bhfmm */




#endif /* CUBOIDPYRAMIDALMATRIX_H_ */
