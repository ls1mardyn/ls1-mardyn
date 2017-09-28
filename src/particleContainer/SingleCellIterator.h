/*
 * SingleCellIterator.h
 *
 *  Created on: 27 Sep 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_SINGLECELLITERATOR_H_
#define SRC_PARTICLECONTAINER_SINGLECELLITERATOR_H_

#include "molecules/Molecule.h"

class ParticleCellBase;

class SingleCellIterator {
public:
	SingleCellIterator();
	SingleCellIterator(ParticleCellBase * cell_arg, size_t index_arg = 0);
	SingleCellIterator& operator=(const SingleCellIterator& other);
	~SingleCellIterator(){}

	void operator ++();

	bool operator == (const SingleCellIterator& other) const;
	bool operator != (const SingleCellIterator& other) const;

	Molecule& operator *  () const;
	Molecule* operator -> () const;

	static SingleCellIterator invalid();

	void deleteCurrentParticle();

	size_t getIndex() const {
		return _mol_index;
	}

	ParticleCellBase * getCell() const { return _cell; }

	void make_invalid ();

private:

	ParticleCellBase * _cell;
	size_t _mol_index;
	bool _currentParticleDeleted;

	// note: possible to add offset for a threaded entry within cells

#ifdef ENABLE_REDUCED_MEMORY_MODE
	Molecule _AoSMoleculeReservoir;
#endif

};

inline SingleCellIterator::SingleCellIterator() : _cell(nullptr), _mol_index(0), _currentParticleDeleted(false) {
	make_invalid();
}

inline SingleCellIterator& SingleCellIterator::operator=(const SingleCellIterator& other) {
	_cell = other._cell;
	_mol_index = other._mol_index;
	_currentParticleDeleted = other._currentParticleDeleted;
	// no operator= for Molecule and it should not be needed.
	return *this;
}

inline bool SingleCellIterator :: operator == (const SingleCellIterator& other) const {
	// _currentParticleDeleted not be needed?
	return (_cell == other._cell) and (_mol_index == other._mol_index);
}

inline bool SingleCellIterator :: operator != (const SingleCellIterator& other) const {
	return not (*this == other);
}

// no clue why this returns a pointer
inline Molecule* SingleCellIterator:: operator -> () const {
	return &(this->operator*());
}

inline SingleCellIterator SingleCellIterator :: invalid () {
	return SingleCellIterator();
}

inline void SingleCellIterator :: make_invalid () {
	_cell = nullptr;
	_mol_index = size_t(-1);
	_currentParticleDeleted = false;
}

#endif /* SRC_PARTICLECONTAINER_SINGLECELLITERATOR_H_ */
