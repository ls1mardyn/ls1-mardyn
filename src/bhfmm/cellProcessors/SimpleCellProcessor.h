/*
 * SimpleCellProcessor.h
 *
 *  Created on: 20 Sep 2016
 *      Author: tchipevn
 */

#ifndef SRC_BHFMM_CELLPROCESSORS_SIMPLECELLPROCESSOR_H_
#define SRC_BHFMM_CELLPROCESSORS_SIMPLECELLPROCESSOR_H_

#include <cstddef>
#include <cmath>
class Molecule;

class ParticleCell;

namespace bhfmm {

/**
 * simpler copy of the main CellProcessor.
 * @author tchipevn
 */
class SimpleCellProcessor {
public:
	SimpleCellProcessor() {}
    virtual ~SimpleCellProcessor() {}

	/**
	 * called before the traversal starts.
	 *
	 * @param numCells number of cells in window
	 */
	virtual void initTraversal() = 0;

	/**
	 * Called when this cell is the current cell.
	 *
	 * @note will not be called for empty cells.
	 */
	virtual void processCell(ParticleCell& cell) = 0;

	/**
	 * Called after the traversal finished.
	 */
	virtual void endTraversal() = 0;
};

} /* namespace bhfmm */

#endif /* SRC_BHFMM_CELLPROCESSORS_SIMPLECELLPROCESSOR_H_ */
