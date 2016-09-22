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

	Molecule& moleculesAt(std::vector<Molecule *>::size_type i) {
		return *_molecules.at(i);
	}

	std::vector<Molecule *>::iterator moleculesBegin() {
		return _molecules.begin();
	}

	std::vector<Molecule *>::const_iterator moleculesCBegin() const {
		return _molecules.cbegin();
	}

	std::vector<Molecule *>::iterator moleculesEnd() {
		return _molecules.end();
	}

	std::vector<Molecule *>::const_iterator moleculesCEnd() const {
		return _molecules.cend();
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

	/**
	 * Returns the cell index
	 */
	unsigned long getCellIndex() {
		return _cellIndex;
	}

	/**
	 * Sets the cell index. On one process, this index must be unique.
	 * @param cellIndex
	 */
	void setCellIndex(unsigned long cellIndex){
		_cellIndex = cellIndex;
	}

	double getBoxMin(int d) const {
		return _boxMin[d];
	}

	void setBoxMin(const double b[3]) {
		for(int d=0; d< 3; ++d) {
			_boxMin[d] = b[d];
		}
	}

	double getBoxMax(int d) const {
		return _boxMax[d];
	}

	void setBoxMax(const double b[3]) {
		for (int d = 0; d < 3; ++d) {
			_boxMax[d] = b[d];
		}
	}

private:
	/**
	 * \brief lower left front corner
	 */
	double _boxMin[3];

	/**
	 * \brief upper right back corner
	 */
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

	/**
	 * \brief The index of the cell.
	 * On one process every index must be unique.
	 */
	unsigned long _cellIndex;
};

} /* namespace bhfmm */

#endif /* SRC_BHFMM_CONTAINERS_PARTICLECELLPOINTERS_H_ */
