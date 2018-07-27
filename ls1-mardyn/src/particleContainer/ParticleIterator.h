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
#include "SingleCellIterator.h"
#include "ParticleCell.h"
#include "WrapOpenMP.h"

class ParticleIterator {
public:
	enum Type {
		ALL_CELLS=0, /* iterates every cell */
		ONLY_INNER_AND_BOUNDARY=1, /* iterates only inner and boundary cells, i.e. no halo cells */
	};

	typedef std::vector<ParticleCell> CellContainer_T;
	typedef CellContainer_T* CellContainer_T_ptr;
	typedef size_t CellIndex_T;
	typedef size_t MolIndex_T;

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

	CellIndex_T getCellIndex(){return _cell_index;}

protected:
	SingleCellIterator _cell_iterator;
	Type _type;
	CellContainer_T_ptr _cells;

	CellIndex_T _cell_index;

	const CellIndex_T _stride;

	virtual void make_invalid ();
	virtual void next_non_empty_cell();
	virtual void updateCellIteratorCell();
};


inline ParticleIterator :: ParticleIterator () : _cell_iterator(), _type(ALL_CELLS), _cells (nullptr), _cell_index (0), _stride (1) {
	make_invalid();
}

inline ParticleIterator :: ParticleIterator (Type t_arg, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize) :
		_cell_iterator(&(cells_arg->front())), _type(t_arg), _cells (cells_arg), _cell_index (offset_arg), _stride (stride_arg) {
	if ( !initialize ) {
		return;
	}

#ifdef ENABLE_REDUCED_MEMORY_MODE
	const unsigned long my_start = _cells->size() * mardyn_get_thread_num() / mardyn_get_num_threads();
	_cell_index = static_cast<CellIndex_T>(my_start);
#endif
	updateCellIteratorCell();

	mardyn_assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;

	if(_cell_index >= cells.size()) {
		make_invalid();
	}
	else {
		if(cells.at(_cell_index).isEmpty()) {
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
	_cell_iterator = other._cell_iterator;
	_type = other._type;
	_cells = other._cells;
	_cell_index = other._cell_index;
	return *this;
}

inline void ParticleIterator :: next_non_empty_cell() {
	mardyn_assert(*this != ParticleIterator :: invalid());
	mardyn_assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;
	const CellIndex_T numCells = cells.size();

	// find the next non-empty cell
	bool validCellFound = false;

#ifndef ENABLE_REDUCED_MEMORY_MODE
	for (_cell_index += _stride; _cell_index < numCells; _cell_index += _stride) {
#else
	const unsigned long my_end = _cells->size() * (mardyn_get_thread_num() + 1) / mardyn_get_num_threads();
	for (_cell_index++; _cell_index < my_end; ++_cell_index) {
#endif

		const ParticleCellBase & c = cells.at(_cell_index);

		if(c.isNotEmpty() and (_type == ALL_CELLS or not c.isHaloCell())) {
			validCellFound = true;
			updateCellIteratorCell();
			break;
		}
	}

	if (not validCellFound) {
		make_invalid();
	}
}

inline void ParticleIterator :: operator ++ () {

	if (not (_cell_iterator == SingleCellIterator::invalid())) {
		++_cell_iterator;
	}

	// don't merge into if-else, _cell_iterator may becoeme invalid after ++

	if (_cell_iterator == SingleCellIterator::invalid()) {
		next_non_empty_cell();
	}
}

inline bool ParticleIterator :: operator == (const ParticleIterator& other) const {
	return (_cell_index == other._cell_index) and (_cell_iterator == other._cell_iterator) and (_cells == other._cells);
}

inline bool ParticleIterator :: operator != (const ParticleIterator& other) const {
	return not (*this == other);
}

inline Molecule& ParticleIterator :: operator * () const {
	// .at method performs automatically an out-of-bounds check
	mardyn_assert(&_cells->at(_cell_index) == dynamic_cast<const ParticleCell*>(_cell_iterator.getCell()));
	return _cell_iterator.operator *();
}

// no clue why this returns a pointer
inline Molecule* ParticleIterator:: operator -> () const {
	return &(this->operator*());
}

inline ParticleIterator ParticleIterator :: invalid () {
	return ParticleIterator();
}

inline void ParticleIterator :: make_invalid () {
	_cells = nullptr;
	_cell_index = 0;
	_cell_iterator.make_invalid();
}

inline void ParticleIterator :: deleteCurrentParticle () {
	_cell_iterator.deleteCurrentParticle();
}

inline void ParticleIterator :: updateCellIteratorCell() {
	if(_cell_index < _cells->size()) {
		_cell_iterator = SingleCellIterator(&_cells->at(_cell_index));
	} else {
		_cell_iterator.make_invalid();
	}
}


#endif /* #ifndef ParticleIterator_INC */
