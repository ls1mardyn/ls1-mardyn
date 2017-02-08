/***********************************************************************************//**
 *
 * \file ParticleIterator.h
 *
 * \brief Iterator class for cell based particle containers
 *
 * \author Wolfgang Hölzl (hoelzlw), hoelzlw AT in.tum.de
 *
 * \details This class implements a forward iterator for the cell based particle containers.
 * It does a coarse iteration over the cells and, within each cell, a fine iteration over the particles.
 *
 **************************************************************************************/

#ifndef  ParticleIterator_INC
#define  ParticleIterator_INC

#include <vector>
#include <stdexcept>
#include <cassert>
#include "ParticleCell.h"


class ParticleIterator {
public:
	typedef std::vector<ParticleCell> CellContainer_T;
	typedef CellContainer_T* CellContainer_T_ptr;
	typedef CellContainer_T::size_type CellIndex_T;
	typedef std::vector<Molecule>::size_type MolIndex_T;

	ParticleIterator ();
	ParticleIterator (CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize=true);
	ParticleIterator& operator=(const ParticleIterator& other);

	virtual ~ParticleIterator(){}

	void operator ++ ();

	bool operator == (const ParticleIterator& other) const;
	bool operator != (const ParticleIterator& other) const;

	Molecule& operator *  () const;
	Molecule* operator -> () const;

	static ParticleIterator invalid ();

	void deleteCurrentParticle();

protected:
	CellContainer_T_ptr _cells;

	CellIndex_T _cell_index;
	MolIndex_T _mol_index;

	bool _currentParticleDeleted;

	const CellIndex_T _stride;

	virtual void make_invalid ();
	virtual void next_non_empty_cell();
};

inline ParticleIterator :: ParticleIterator () : _cells (nullptr), _cell_index (0), _mol_index (0), _currentParticleDeleted(false), _stride (1) {
	make_invalid();
}

inline ParticleIterator :: ParticleIterator (CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize) :
		_cells (cells_arg), _cell_index (offset_arg), _mol_index (0), _currentParticleDeleted(false), _stride (stride_arg) {
	if ( !initialize ) {
		return;
	}

	assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;

	if(_cell_index >= cells.size()) {
		make_invalid();
	}
	else {
		if(cells[_cell_index].isEmpty()) {
			next_non_empty_cell();
		}
		/*
		else {
			// leave object as is
		}
		*/
	}
}

inline ParticleIterator& ParticleIterator::operator=(const ParticleIterator& other) {
	assert(_stride == other._stride);
	_cells = other._cells;
	_cell_index = other._cell_index;
	_mol_index = other._mol_index;
	return *this;
}

inline void ParticleIterator :: next_non_empty_cell() {
	assert(*this != ParticleIterator :: invalid());
	assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;
	const CellIndex_T numCells = cells.size();

	// find the next non-empty cell
	do {
		_cell_index += _stride;
	} while (_cell_index < numCells and cells[_cell_index].isEmpty());

	if (_cell_index < numCells) {
		// if we found a non-empty cell..
		// valid
		_mol_index = 0;
	}
	else {
		// else there is no next non-empty cell..
		// invalid
		make_invalid();
	}
}

inline void ParticleIterator :: operator ++ () {

	// if the "current" particle was deleted, then there is a new particle at _mol_index
	// and the _mol_index value should not be incremented.
	_mol_index += _currentParticleDeleted ? 0 : 1;

	const CellContainer_T& cells = *_cells;

	if (_mol_index >= static_cast<MolIndex_T>(cells[_cell_index].getMoleculeCount())) {
		next_non_empty_cell();
	}

	// at this stage, we are pointing to a new particle, or are invalid, in any case we can set this:
	_currentParticleDeleted = false;
}

inline bool ParticleIterator :: operator == (const ParticleIterator& other) const {
	return (_cell_index == other._cell_index) and (_mol_index == other._mol_index);
}

inline bool ParticleIterator :: operator != (const ParticleIterator& other) const {
	return not (*this == other);
}

inline Molecule& ParticleIterator :: operator * () const {
	// .at method performs automatically an out-of-bounds check
	return _cells->at(_cell_index).moleculesAt(_mol_index);
}

// no clue why this returns a pointer
inline Molecule* ParticleIterator:: operator -> () const {
	return &(this->operator*());
}

inline ParticleIterator ParticleIterator :: invalid () {
	return ParticleIterator();
}

inline void ParticleIterator :: make_invalid () {
	_cell_index = CellIndex_T(-1);
	_mol_index = MolIndex_T(-1);
}

inline void ParticleIterator :: deleteCurrentParticle () {
	_cells->at(_cell_index).deleteMoleculeByIndex(_mol_index);
	_currentParticleDeleted = true;
}
#endif /* #ifndef ParticleIterator_INC */
