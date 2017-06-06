/*
 * ParticleIterator.cpp
 *
 *  Created on: 23 May 2017
 *      Author: tchipevn
 */

#include "ParticleIterator.h"
#include "ParticleContainer.h"

ParticleIterator :: ParticleIterator () : _type(ALL_CELLS), _cells (nullptr), _cell_index (0), _mol_index (0), _currentParticleDeleted(false), _stride (1) {
	make_invalid();
}

ParticleIterator :: ParticleIterator (Type t_arg, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize) :
		_type(t_arg), _cells (cells_arg), _cell_index (offset_arg), _mol_index (0), _currentParticleDeleted(false), _stride (stride_arg) {
	if ( !initialize ) {
		return;
	}

	mardyn_assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;

	if(_cell_index >= cells.getNumCells()) {
		make_invalid();
	}
	else {
		if(cells.getCell(_cell_index)->isEmpty()) {
			next_non_empty_cell();
		}
		/*
		else {
			// leave object as is
		}
		*/
	}
}

ParticleIterator& ParticleIterator::operator=(const ParticleIterator& other) {
	mardyn_assert(_stride == other._stride);
	_type = other._type;
	_cells = other._cells;
	_cell_index = other._cell_index;
	_mol_index = other._mol_index;
	_currentParticleDeleted = other._currentParticleDeleted;
	return *this;
}

void ParticleIterator :: next_non_empty_cell() {
	mardyn_assert(*this != ParticleIterator :: invalid());
	mardyn_assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;
	const CellIndex_T numCells = cells.getNumCells();

	// find the next non-empty cell
	bool validCellFound = false;

	for (_cell_index += _stride; _cell_index < numCells; _cell_index += _stride) {

		const ParticleCellBase * c = cells.getCell(_cell_index);

		if(c->isNotEmpty() and (_type == ALL_CELLS or not c->isHaloCell())) {
			validCellFound = true;
			_mol_index = 0;
			break;
		}
	}

	if (not validCellFound) {
		make_invalid();
	}
}

void ParticleIterator :: operator ++ () {

	// if the "current" particle was deleted, then there is a new particle at _mol_index
	// and the _mol_index value should not be incremented.
	_mol_index += _currentParticleDeleted ? 0 : 1;

	const CellContainer_T& cells = *_cells;

	if (_mol_index >= static_cast<MolIndex_T>(cells.getCell(_cell_index)->getMoleculeCount())) {
		next_non_empty_cell();
	}

	// at this stage, we are pointing to a new particle, or are invalid, in any case we can set this:
	_currentParticleDeleted = false;
}

bool ParticleIterator :: operator == (const ParticleIterator& other) const {
	return (_cell_index == other._cell_index) and (_mol_index == other._mol_index);
}

bool ParticleIterator :: operator != (const ParticleIterator& other) const {
	return not (*this == other);
}

Molecule& ParticleIterator :: operator * () const {
	// .at method performs automatically an out-of-bounds check
	return _cells->getCell(_cell_index)->moleculesAt(_mol_index);
}

// no clue why this returns a pointer
Molecule* ParticleIterator:: operator -> () const {
	return &(this->operator*());
}

ParticleIterator ParticleIterator :: invalid () {
	return ParticleIterator();
}

void ParticleIterator :: make_invalid () {
	_cell_index = CellIndex_T(-1);
	_mol_index = MolIndex_T(-1);
	_currentParticleDeleted = false;
}

void ParticleIterator :: deleteCurrentParticle () {
	_cells->getCell(_cell_index)->deleteMoleculeByIndex(_mol_index);
	_currentParticleDeleted = true;
}


