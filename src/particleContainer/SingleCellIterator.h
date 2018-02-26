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

	void next() {
		operator++();
	}

	Molecule& operator *  () const;
	Molecule* operator -> () const;

	void deleteCurrentParticle();

	size_t getIndex() const {
		return _mol_index;
	}

	ParticleCellBase * getCell() const { return _cell; }

	bool hasNext() const {
		return isValid();
	}


private:
	void operator ++();
	bool isValid() const;

	ParticleCellBase * _cell;
	size_t _mol_index;
	bool _currentParticleDeleted;

	// note: possible to add offset for a threaded entry within cells

#ifdef ENABLE_REDUCED_MEMORY_MODE
	Molecule _AoSMoleculeReservoir;
#endif

};

inline SingleCellIterator::SingleCellIterator() : _cell(nullptr), _mol_index(0), _currentParticleDeleted(false) {
}

inline SingleCellIterator& SingleCellIterator::operator=(const SingleCellIterator& other) {
	_cell = other._cell;
	_mol_index = other._mol_index;
	_currentParticleDeleted = other._currentParticleDeleted;
	// no operator= for Molecule and it should not be needed.
	return *this;
}

// no clue why this returns a pointer
inline Molecule* SingleCellIterator:: operator -> () const {
	return &(this->operator*());
}

#endif /* SRC_PARTICLECONTAINER_SINGLECELLITERATOR_H_ */
