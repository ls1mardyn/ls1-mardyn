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

class Molecule;

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
};

} /* namespace bhfmm */

#endif /* SRC_BHFMM_CONTAINERS_PARTICLECELLPOINTERS_H_ */
