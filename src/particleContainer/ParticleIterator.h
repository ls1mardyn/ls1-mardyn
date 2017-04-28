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
#include "utils/mardyn_assert.h"
#include "ParticleCell.h"

class ParticleIterator {
public:
	enum Type {
		ALL_CELLS=0, /* iterates every cell */
		ONLY_INNER_AND_BOUNDARY=1, /* iterates every cell except halo cells */
	};

	typedef std::vector<ParticleCell> CellContainer_T;
	typedef CellContainer_T* CellContainer_T_ptr;
	typedef CellContainer_T::size_type CellIndex_T;
	typedef std::vector<Molecule>::size_type MolIndex_T;

	ParticleIterator ();
	ParticleIterator (Type t_arg, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize=true);
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
	Type _type;
	CellContainer_T_ptr _cells;

	CellIndex_T _cell_index;
	MolIndex_T _mol_index;

	bool _currentParticleDeleted;

	const CellIndex_T _stride;

	virtual void make_invalid ();
	virtual void next_non_empty_cell();
};

inline ParticleIterator :: ParticleIterator () : _type(ALL_CELLS), _cells (nullptr), _cell_index (0), _mol_index (0), _currentParticleDeleted(false), _stride (1) {
	make_invalid();
}

inline ParticleIterator :: ParticleIterator (Type t_arg, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize) :
		_type(t_arg), _cells (cells_arg), _cell_index (offset_arg), _mol_index (0), _currentParticleDeleted(false), _stride (stride_arg) {
	if ( !initialize ) {
		return;
	}

	mardyn_assert(_cells != nullptr);

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
	mardyn_assert(_stride == other._stride);
	_type = other._type;
	_cells = other._cells;
	_cell_index = other._cell_index;
	_mol_index = other._mol_index;
	_currentParticleDeleted = other._currentParticleDeleted;
	return *this;
}

inline void ParticleIterator :: next_non_empty_cell() {
	mardyn_assert(*this != ParticleIterator :: invalid());
	mardyn_assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;
	const CellIndex_T numCells = cells.size();

	// find the next non-empty cell
	bool validCellFound = false;

	for (_cell_index += _stride; _cell_index < numCells; _cell_index += _stride) {

		const ParticleCell & c = cells[_cell_index];

		if(c.isNotEmpty() and (_type == ALL_CELLS or not c.isHaloCell())) {
			validCellFound = true;
			_mol_index = 0;
			break;
		}
	}

	if (not validCellFound) {
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
	_currentParticleDeleted = false;
}

inline void ParticleIterator :: deleteCurrentParticle () {
	_cells->at(_cell_index).deleteMoleculeByIndex(_mol_index);
	_currentParticleDeleted = true;
}

#endif /* #ifndef ParticleIterator_INC */
