/*
 * ParticleCellPointers.h
 *
 *  Created on: 20 Sep 2016
 *      Author: tchipevn
 */

#ifndef SRC_BHFMM_CONTAINERS_PARTICLECELLPOINTERS_H_
#define SRC_BHFMM_CONTAINERS_PARTICLECELLPOINTERS_H_

#include <vector>

#include "particleContainer/Cell.h"
#include "particleContainer/adapter/CellDataSoA.h"

#include "molecules/MoleculeForwardDeclaration.h"

#ifdef QUICKSCHED
#include "quicksched.h"
#endif
/**
 * simpler copy of the class ParticleCellPointers.
 * Uses pointers, instead of references for main storage,
 * in order to save storage.
 *
 * Important! this class is not responsible for the storage
 * which is pointed by the pointers! I.e. no deallocation
 * of the storage should be performed!
 *
 * @author tchipevn
 */

namespace bhfmm {

#ifdef QUICKSCHED
struct taskData{
	qsched_res_t  _resourceId;
	qsched_task_t _preprocessId;
	qsched_task_t _postprocessId;
	qsched_task_t _P2PId;
};
#endif

class ParticleCellPointers: public Cell {
public:
	/**
	 * \brief Initialize data pointers to 0.
	 * \author Johannes Heckl
	 */
	ParticleCellPointers();
	/**
	 * \brief Destructor.
	 * \author Johannes Heckl
	 */
	~ParticleCellPointers() ;

	//! removes all elements from the list molecules without deallocating them
	void removeAllParticles();

	//! insert a single molecule into this cell
	bool addParticle(Molecule* particle_ptr);

	Molecule& moleculesAt(size_t i) {
		return *_molecules.at(i);
	}

	bool isEmpty() const;

	//! return the number of molecules contained in this cell
	int getMoleculeCount() const;

	/**
	 * \brief Get the structure of arrays for VectorizedCellProcessor.
	 * \author Johannes Heckl
	 */
	CellDataSoA& getCellDataSoA() {
		return _cellDataSoA;
	}

#ifdef QUICKSCHED
	taskData& getTaskData(){
		return _taskData;
	}
    void setResourceId(qsched_res_t id){
        _taskData._resourceId = id;
    }
    void setP2PId(qsched_res_t id){
        _taskData._P2PId = id;
    }
    void setPreprocessId(qsched_res_t id){
        _taskData._preprocessId = id;
    }
    void setPostprocessId(qsched_res_t id){
        _taskData._postprocessId = id;
    }
#endif // QUICKSCHED

	void assignCellToHaloRegion() { haloCell = true; }
	void assignCellToBoundaryRegion() { boundaryCell = true; }
	void assignCellToInnerRegion() { innerCell = true; }
	void assignCellToInnerMostAndInnerRegion() { innerCell = true; innerMostCell = true; }

	void skipCellFromHaloRegion() { haloCell = false; }
	void skipCellFromBoundaryRegion() { boundaryCell = false; }
	void skipCellFromInnerRegion() { innerCell = false; }
	void skipCellFromInnerMostRegion() { innerMostCell = false; }

	bool isHaloCell() const { return haloCell; }
	bool isBoundaryCell() const { return boundaryCell; }
	bool isInnerCell() const { return innerCell; }
	bool isInnerMostCell() const { return innerMostCell; }

	double getBoxMin(int d) const { return _boxMin[d]; }
	double getBoxMax(int d) const { return _boxMax[d]; }
	void setBoxMin(const double b[3]) { for (int d = 0; d < 3; ++d) { _boxMin[d] = b[d]; } }
	void setBoxMax(const double b[3]) { for (int d = 0; d < 3; ++d) { _boxMax[d] = b[d]; } }

private:
	//! true when the cell is in the halo region
	bool haloCell;
	//! true when the cell is in the boundary region
	bool boundaryCell;
	//! true when the cell is in the inner region. Innermost cells are always also innerCells.
	bool innerCell;
	//! true when the cell is in the innermost region (does not have neighbors, that are boundary cells)
	bool innerMostCell;
	//! lower left front corner
	double _boxMin[3];
	//! upper right back corner
	double _boxMax[3];

	/**
	 * \brief A vector of pointers to the Molecules in this cell.
	 */
	std::vector<Molecule *> _molecules;

	/**
	 * \brief Structure of arrays for VectorizedCellProcessor.
	 * \author Johannes Heckl
	 */
	CellDataSoA _cellDataSoA;

#ifdef QUICKSCHED
	struct taskData _taskData;
#endif // QUICKSCHED
};

} /* namespace bhfmm */

#endif /* SRC_BHFMM_CONTAINERS_PARTICLECELLPOINTERS_H_ */
