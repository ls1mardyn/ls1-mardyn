/*
 * CellBorderManager.h
 *
 *  Created on: 10 Oct 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_CELLBORDERANDFLAGMANAGER_H_
#define SRC_PARTICLECONTAINER_CELLBORDERANDFLAGMANAGER_H_

#include <array>
#include <vector>

#include "utils/mardyn_assert.h"
#include "utils/threeDimensionalMapping.h"

// TODO: use in LinkedCells::initializeCells

class CellBorderAndFlagManager {
public:
	typedef unsigned long GlobalLinearizedIndex_t;
	typedef unsigned long Global1DIndex_t;
	typedef std::array<Global1DIndex_t, 3> Global3DIndex_t;
	enum class IsCell_t { HALO = 0, BOUNDARY = 1, INNER = 2, INNERMOST = 3};

	CellBorderAndFlagManager() : _initCalled(false) {}

	void init(int cellsPerDim[3],
			double haloBoxMin[3], double haloBoxMax[3],
			double boxMin[3], double boxMax[3],
			double cellLength[3]) {
		for (int d = 0; d < 3; ++d) {
			_cellsPerDimension[d] = cellsPerDim[d];
			_haloBoundingBoxMin[d] = haloBoxMin[d];
			_haloBoundingBoxMax[d] = haloBoxMax[d];
			_boundingBoxMin[d] = boxMin[d];
			_boundingBoxMax[d] = boxMax[d];
			_cellLength[d] = cellLength[d];
		}
		_initCalled = true;
	}

	bool isHaloCell(GlobalLinearizedIndex_t cellIndex) const { return isCell(IsCell_t::HALO, cellIndex); }
	bool isBoundaryCell(GlobalLinearizedIndex_t cellIndex) const { return isCell(IsCell_t::BOUNDARY, cellIndex); }
	bool isInnerCell(GlobalLinearizedIndex_t cellIndex) const { return isCell(IsCell_t::INNER, cellIndex); }
	bool isInnerMostCell(GlobalLinearizedIndex_t cellIndex) const { return isCell(IsCell_t::INNERMOST, cellIndex); }

	double getBoundingBoxMin(GlobalLinearizedIndex_t cellIndex, int dimension) const {
		Global3DIndex_t r = rIndex(cellIndex);
		Global1DIndex_t indexInDimension = r[dimension];
		return getBorder(indexInDimension, dimension);
	}

	double getBoundingBoxMax(GlobalLinearizedIndex_t cellIndex, int dimension) const {
		Global3DIndex_t r = rIndex(cellIndex);
		Global1DIndex_t indexInDimension = r[dimension];
		return getBorder(indexInDimension + 1, dimension);
	}

private:
	bool isCell(IsCell_t type, GlobalLinearizedIndex_t cellIndex) const {
		bool ret;
		Global3DIndex_t r = rIndex(cellIndex);

		enum {HALOORBOUNDARY = 0, INNERORINNERMOST=1} t;
		if (type == IsCell_t::HALO or type == IsCell_t::BOUNDARY) {
			t = HALOORBOUNDARY;
		} else {
			t = INNERORINNERMOST;
		}
		switch(t) {
		case HALOORBOUNDARY:
			// if halo or boundary in any dimension, then true
			ret = false;
			for (int d = 0; d < 3; ++d){
				ret = ret or isCell1D(type, r[d], d);
			}
			break;
		case INNERORINNERMOST:
			// if inner or innermost in all dimensions, then true
			ret = true;
			for (int d = 0; d < 3; ++d){
				ret = ret and isCell1D(type, r[d], d);
			}
			break;
		}
		return ret;
	}
	bool isCell1D(IsCell_t type, Global1DIndex_t index, int dimension) const {
		mardyn_assert(_initCalled);
		bool ret;
		switch(type) {
		case IsCell_t::HALO:
			ret = (index == 0) or (index == _cellsPerDimension[dimension]-1);
			break;
		case IsCell_t::BOUNDARY:
			ret = (index == 1) or (index == _cellsPerDimension[dimension]-2);
			break;
		case IsCell_t::INNER:
			ret = (index >= 2) and (index <= _cellsPerDimension[dimension]-3);
			break;
		case IsCell_t::INNERMOST:
			ret = (index >= 3) and (index <= _cellsPerDimension[dimension]-4);
			break;
		}
		return ret;
	}
	double getBorder(Global1DIndex_t index, int dimension) const {
		mardyn_assert(_initCalled);
		mardyn_assert(index <= _cellsPerDimension[dimension]);

		double ret;
		// why oh why don't we have switch statements on variables
		if (index == 0) {
			ret = _haloBoundingBoxMin[dimension];
		} else if (index == 1) {
			ret = _boundingBoxMin[dimension];
		} else if (index == _cellsPerDimension[dimension] - 1) {
			ret = _boundingBoxMax[dimension];
		} else if (index == _cellsPerDimension[dimension]) {
			ret = _haloBoundingBoxMax[dimension];
		} else {
			mardyn_assert(index > 1 and index < _cellsPerDimension[dimension] - 1);
			ret = index * _cellLength[dimension] + _haloBoundingBoxMin[dimension];
		}
		return ret;
	}

	Global3DIndex_t rIndex(GlobalLinearizedIndex_t index) const {
		return threeDimensionalMapping::oneToThreeD(index, _cellsPerDimension);
	}

	//member fields
	bool _initCalled;
	Global3DIndex_t _cellsPerDimension;
	std::array<double, 3> _boundingBoxMin, _boundingBoxMax;
	std::array<double, 3> _haloBoundingBoxMin, _haloBoundingBoxMax;
	std::array<double, 3> _cellLength;

};

#endif /* SRC_PARTICLECONTAINER_CELLBORDERANDFLAGMANAGER_H_ */
