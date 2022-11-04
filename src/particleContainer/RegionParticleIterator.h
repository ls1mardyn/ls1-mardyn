/**************************************************************************************
 *
 * \file RegionParticleIterator.h
 *
 * \brief Iterator class for cell based particle containers
 *
 * \author Andrei Costinescu (costines), costines AT in.tum.de
 *
 * \details This class implements a forward iterator for the cell based particle containers.
 * It does a coarse iteration over the cells in a specific region and, within each cell, a fine iteration over the particles.
 *
 **************************************************************************************/

#pragma once

#ifdef MARDYN_AUTOPAS

#include "ParticleIterator.h"
class RegionParticleIterator : public ParticleIterator {
public:
	RegionParticleIterator() = default;
	explicit RegionParticleIterator(const autopas::IteratorTraits<AutoPasSimpleMolecule>::iterator_t& parent)
		: ParticleIterator(parent){};
};

#else

#include <vector>
#include <stdexcept>
#include "utils/mardyn_assert.h"
#include "ParticleIterator.h"

class RegionParticleIterator : public ParticleIterator {
	public:
		RegionParticleIterator ();
		RegionParticleIterator (Type t, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const int startCellIndex_arg, const int regionDimensions_arg[3], const int globalDimensions_arg[3], const double startRegion_arg[3], const double endRegion_arg[3]);
		RegionParticleIterator& operator=(const RegionParticleIterator& other);
		
		~RegionParticleIterator() override = default;

		constexpr RegionParticleIterator(const RegionParticleIterator&) = default;

		void operator++() override;
	private:

		CellIndex_T getGlobalCellIndex();
		void next_non_empty_cell() override;

		CellIndex_T _baseX;
		CellIndex_T _baseY;
		CellIndex_T _baseZ;
		CellIndex_T _localCellIndex;
		CellIndex_T _regionDimensions[3];
		CellIndex_T _globalDimensions[3];

		double _startRegion[3];
		double _endRegion[3];
};

inline RegionParticleIterator :: RegionParticleIterator () : ParticleIterator(), _localCellIndex(0) {
}

inline RegionParticleIterator :: RegionParticleIterator (Type t, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const int startCellIndex_arg, const int regionDimensions_arg[3], const int globalDimensions_arg[3], const double startRegion_arg[3], const double endRegion_arg[3]) :
		ParticleIterator(t, cells_arg, offset_arg, stride_arg, false), _localCellIndex (offset_arg) {

#ifdef ENABLE_REDUCED_MEMORY_MODE
	_cell_index = offset_arg;
#endif

	for (int d = 0; d < 3; d++) {
		_regionDimensions[d] = regionDimensions_arg[d];
		_globalDimensions[d] = globalDimensions_arg[d];
		_startRegion[d] = startRegion_arg[d];
		_endRegion[d] = endRegion_arg[d];
	}

	_baseZ = startCellIndex_arg / (_globalDimensions[0] * _globalDimensions[1]);
	_baseY = (startCellIndex_arg % (_globalDimensions[0] * _globalDimensions[1])) / _globalDimensions[0];
	_baseX = (startCellIndex_arg % (_globalDimensions[0] * _globalDimensions[1])) % _globalDimensions[0];

	_cell_index = getGlobalCellIndex();
	updateCellIteratorCell();

	mardyn_assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;
	const CellIndex_T numCellsInRegion = _regionDimensions[2] * _regionDimensions[1] * _regionDimensions[0]; //cells.size();

	// if _cell_index is out of the region => invalid <=> (_localCellIndex >= numCellsInRegion)
	if (_localCellIndex < numCellsInRegion) {
		/*
		if (cells[_cell_index].isEmpty()) {
			next_non_empty_cell();
		}
		*/
//		_currentParticleDeleted = true;

		if (cells[_cell_index].isNotEmpty() and (this->operator*()).inBox(_startRegion, _endRegion)) {
			// if the particle is good, leave it
		} else {
			// else, find a particle in the box
			this->operator++();
		}
	} else {
		// set to invalid
		_cells = nullptr;
	}
}

inline RegionParticleIterator& RegionParticleIterator::operator=(const RegionParticleIterator& other) {
	mardyn_assert(_stride == other._stride);
	_cells = other._cells;
	_cell_index = other._cell_index;
	_cell_iterator = other._cell_iterator;
	_baseX = other._baseX;
	_baseY = other._baseY;
	_baseZ = other._baseZ;
	_localCellIndex = other._localCellIndex;
	for(int d = 0; d < 3; d++){
		_regionDimensions[d] = other._regionDimensions[d];
		_globalDimensions[d] = other._globalDimensions[d];
		_startRegion[d] = other._startRegion[d];
		_endRegion[d] = other._endRegion[d];
	}
	return *this;
}

inline void RegionParticleIterator :: operator ++() {
	do{
		ParticleIterator :: operator++();
	} while (isValid() and !(this->operator*()).inBox(_startRegion, _endRegion));
}

inline void RegionParticleIterator :: next_non_empty_cell() {
	//cellIndex should always be the index in the cell array (_cells member variable in LinkedCells)
	mardyn_assert(_cells != nullptr);

	const CellContainer_T& cells = *_cells;
	const CellIndex_T numCellsInRegion = _regionDimensions[2] * _regionDimensions[1] * _regionDimensions[0];

	// find the next non-empty cell
	for (_localCellIndex += _stride; _localCellIndex < numCellsInRegion; _localCellIndex += _stride) {

		const ParticleCellBase & c = cells.at(getGlobalCellIndex());

		if (c.isNotEmpty() and (_type == ALL_CELLS or not c.isHaloCell())) {
			_cell_index = getGlobalCellIndex();
			updateCellIteratorCell();
			break;
		}
	}
}

inline ParticleIterator::CellIndex_T RegionParticleIterator :: getGlobalCellIndex() {
	CellIndex_T dz = _localCellIndex / (_regionDimensions[0] * _regionDimensions[1]);
	CellIndex_T dy = (_localCellIndex % (_regionDimensions[0] * _regionDimensions[1])) / _regionDimensions[0];
	CellIndex_T dx = (_localCellIndex % (_regionDimensions[0] * _regionDimensions[1])) % _regionDimensions[0];
	return (_baseX + dx) + (_baseY + dy) * _globalDimensions[0] + (_baseZ + dz) * _globalDimensions[0] * _globalDimensions[1];
}

#endif  // MARDYN_AUTOPAS
