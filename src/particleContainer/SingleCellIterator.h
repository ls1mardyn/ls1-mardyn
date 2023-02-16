/*
 * SingleCellIterator.h
 *
 *  Created on: 27 Sep 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_SINGLECELLITERATOR_H_
#define SRC_PARTICLECONTAINER_SINGLECELLITERATOR_H_

#include "molecules/Molecule.h"

template<class ParticleCell>
class SingleCellIterator {
public:
	SingleCellIterator();
	SingleCellIterator(ParticleCell * cell_arg,
		size_t index_arg = 0) : _cell(cell_arg), _mol_index(index_arg), _currentParticleDeleted(false) {
	}
	SingleCellIterator& operator=(const SingleCellIterator& other);
	~SingleCellIterator(){}

	constexpr SingleCellIterator(const SingleCellIterator&) = default;

	Molecule& operator *  () const {
		// .at method performs automatically an out-of-bounds check
		Molecule *moleculePtr = nullptr;

	#ifdef ENABLE_REDUCED_MEMORY_MODE
		moleculePtr = const_cast<Molecule *>(& _AoSMoleculeReservoir);
	#endif

		_cell->moleculesAtNew(_mol_index, moleculePtr);

		return *moleculePtr;
	}
	Molecule* operator -> () const;

	void deleteCurrentParticle() {
		_cell->deleteMoleculeByIndex(_mol_index);
		_currentParticleDeleted = true;
	}

	size_t getIndex() const {
		return _mol_index;
	}

	ParticleCell * getCell() const { return _cell; }

	bool isValid() const {
		return _cell != nullptr and _mol_index < static_cast<size_t>(_cell->getMoleculeCount());
	}

	void operator ++() {
		// if the "current" particle was deleted, then there is a new particle at _mol_index
		// and the _mol_index value should not be incremented.
		_mol_index += _currentParticleDeleted ? 0 : 1;

		_currentParticleDeleted = false;
	}

private:

	ParticleCell * _cell;
	size_t _mol_index;
	bool _currentParticleDeleted;

	// note: possible to add offset for a threaded entry within cells

#ifdef ENABLE_REDUCED_MEMORY_MODE
	Molecule _AoSMoleculeReservoir;
#endif

};

template<class ParticleCell>
inline SingleCellIterator<ParticleCell>::SingleCellIterator() : _cell(nullptr), _mol_index(0), _currentParticleDeleted(false) {
}

template<class ParticleCell>
inline SingleCellIterator<ParticleCell>& SingleCellIterator<ParticleCell>::operator=(const SingleCellIterator& other) {
	_cell = other._cell;
	_mol_index = other._mol_index;
	_currentParticleDeleted = other._currentParticleDeleted;
	// no operator= for Molecule and it should not be needed.
	return *this;
}

// no clue why this returns a pointer
template<class ParticleCell>
inline Molecule* SingleCellIterator<ParticleCell>:: operator -> () const {
	return &(this->operator*());
}

#endif /* SRC_PARTICLECONTAINER_SINGLECELLITERATOR_H_ */
