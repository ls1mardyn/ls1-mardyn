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

#include "quicksched.h"
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

struct taskData{
	qsched_res_t  _resourceId;
	qsched_task_t _preprocessId;
	qsched_task_t _postprocessId;
	qsched_task_t _P2PId;
};

class ParticleCellPointers: public Cell {
public:
	/**
	 * \brief Initialize data pointers to 0.
	 * \author Johannes Heckl
	 */
	ParticleCellPointers() ;
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

private:
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
