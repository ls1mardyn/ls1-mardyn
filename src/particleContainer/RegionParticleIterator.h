/***********************************************************************************//**
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

#ifndef  RegionParticleIterator_INC
#define  RegionParticleIterator_INC

#include <vector>
#include <stdexcept>
#include "utils/mardyn_assert.h"
#include "ParticleIterator.h"

class RegionParticleIterator : public ParticleIterator {
	public:
		RegionParticleIterator ();
		RegionParticleIterator (Type t, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const int startCellIndex_arg, const int regionDimensions_arg[3], const int globalDimensions_arg[3], const double startRegion_arg[3], const double endRegion_arg[3]);
		RegionParticleIterator& operator=(const RegionParticleIterator& other);

		void operator ++ ();

		static RegionParticleIterator invalid();

	private:
		CellIndex_T getGlobalCellIndex();
		void make_invalid();
		void next_non_empty_cell();

		CellIndex_T _baseX;
		CellIndex_T _baseY;
		CellIndex_T _baseZ;
		CellIndex_T _localCellIndex;
		CellIndex_T _regionDimensions[3];
		CellIndex_T _globalDimensions[3];

		double _startRegion[3];
		double _endRegion[3];
};

#endif /* #ifndef RegionParticleIterator_INC */
