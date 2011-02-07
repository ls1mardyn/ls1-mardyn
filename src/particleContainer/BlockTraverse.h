/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef BLOCKTRAVERSE_H_
#define BLOCKTRAVERSE_H_

#include <vector>

class ParticleContainer;
class ParticlePairsHandler;
class Cell;

//! @brief BlockTraverse datastructure
//! @author Johannes Weißl
//!
//! This datastructure is there to prevent code duplication in LinkedCells and AdaptiveSubCells classes.
//! The traversePairs() method in these classes are very similar, except for
//! - naming (cells <-> subCells), etc.
//! - offsets (NeighbourOffsets equal for every cell <-> different for every cell)
//! This class contains a generic traversePairs() method and all needed member variables. There are two
//! constructors to deal with the different datastructures for the offsets in LinkedCells and AdaptiveSubCells.
class BlockTraverse {
public:
	//! @brief initialize BlockTraverse structure
	//!
	//! Use this initializer if there already is vector of vectors for the neighbourOffsets.
	//! A BlockTraverse structure created this way doesn't need to be updated using assignOffsets().
	BlockTraverse(
            ParticleContainer* moleculeContainer, std::vector<Cell>& cells,
            std::vector<unsigned long>& innerCellIndices,
            std::vector<unsigned long>& boundaryCellIndices,
            std::vector<unsigned long>& haloCellIndices,
            std::vector<std::vector<unsigned long> >& forwardNeighbourOffsets,
            std::vector<std::vector<unsigned long> >& backwardNeighbourOffsets
	);

	//! @brief initialize BlockTraverse structure
	//!
	//! _(forward|backward)NeighbourOffsets members are newly allocated and must be kept updated using the assignOffsets() method.
	BlockTraverse(
	    ParticleContainer* moleculeContainer, std::vector<Cell>& cells,
	    std::vector<unsigned long>& innerCellIndices, std::vector<unsigned long>& boundaryCellIndices, std::vector<unsigned long>& haloCellIndices
	);

	//! Destructor
	~BlockTraverse();

	//! @brief calculate the forces between the molecules.
	//!
	//! Only molecules with a distance not larger than the cutoff radius are to be used. \n
	//! Only forces on the Molecules which are in the inner and boundary region have to be calculated
	//! Newton's third law should be used for faster computation:
	//! \li a loop over all inner cells calculates the forces with all forward neighbour cells
	//!     all forward cells have to be used, as none of them can be halo or outside
	//! \li a loop over the boundary cells first calculates forces with all forward cells and all
	//!     backward cells. Here it has to be checked whether the neighbour cell is halo or not.
	//!     If it is Halo, the force is calculated, if it isn't, the force is not calculated,
	//!     because the same pair of cells has already been processed in one of the other loops.
	void traversePairs();

	//! @brief assign new (forward|backward)NeighbourOffsets
	//!
	//! Assign the same forwardNeighbourOffsets and backwardNeighbourOffsets to ALL cells of the container.
	void assignOffsets(std::vector<unsigned long>& forwardNeighbourOffsets, std::vector<unsigned long>& backwardNeighbourOffsets);

private:
	//####################################
	//##### PRIVATE MEMBER VARIABLES #####
	//####################################

	//! Datastructure for finding neighbours efficiently
	ParticleContainer* _moleculeContainer;

	//! A ParticlePairsHandler is used to process pairs of particles
	ParticlePairsHandler* _particlePairsHandler;

	//! Vector containing all cells (including halo)
	std::vector<Cell>& _cells;

	//! Vector containing the indices (for the cells vector) of all inner cells (without boundary)
	std::vector<unsigned long>& _innerCellIndices;
	//! Vector containing the indices (for the cells vector) of all boundary cells
	std::vector<unsigned long>& _boundaryCellIndices;
	//! Vector containing the indices (for the cells vector) of all halo cells
	std::vector<unsigned long>& _haloCellIndices;

	//! Neighbours that come in the total ordering after a cell
	std::vector<std::vector<unsigned long> >* _forwardNeighbourOffsets;
	//! Neighbours that come in the total ordering before a cell
	std::vector<std::vector<unsigned long> >* _backwardNeighbourOffsets;

	//! If true, _(forward|backward)NeighbourOffsets are allocated and must be free'd
	bool _allocatedOffsets;
};

#endif /*BLOCKTRAVERSE_H_*/
