
#ifndef CELLPROCESSOR_H_
#define CELLPROCESSOR_H_

#include <cstddef>
#include <cmath>

#include "molecules/MoleculeForwardDeclaration.h"
#include "particleContainer/ParticleCellForwardDeclaration.h"

/**
 * Interface for traversal of cells to allow a cell-wise treatment of molecules.
 *
 * As depicted in the picture, cells are traversed in a first-in-first-out (FIFO) order,
 * e.g. cell 26 is searched for interacting particles before cell 27 is ever considered,
 * and similarily cell no. 26 will become the "current" cell before cell no. 27.
 *
 * \image CellsTraversal.jpg
 *
 * The grey cells are the ones which may e.g. be searched for interacting particles,
 * so they are called <b>active</b> cells. Before a cell will be handed to the
 * method processCell() or processCellPair() for the first time during an iteration,
 * preprocessCell will be called to allow preparatory stuff. PostprocessCell is called
 * analogously after a cell has been considered for the last time. In between,
 * processCell and processCellPair are called.
 *
 * @author eckhardw
 */
class CellProcessor {
protected:
	double _cutoffRadiusSquare;
	double _LJCutoffRadiusSquare;

public:
	CellProcessor(const double cutoffRadius, const double LJCutoffRadius) :
		_cutoffRadiusSquare(cutoffRadius * cutoffRadius),
		_LJCutoffRadiusSquare(LJCutoffRadius * LJCutoffRadius) {}
    /**
     * virtual destructor
     */
	virtual ~CellProcessor() {}

	double getCutoffRadius() const {return sqrt(_cutoffRadiusSquare);}
	double getLJCutoffRadius() const {return sqrt(_LJCutoffRadiusSquare);}
	void setCutoffRadius(const double c) {_cutoffRadiusSquare = c * c;}
	void setLJCutoffRadius(const double ljc) {_LJCutoffRadiusSquare = ljc * ljc;}


	double getCutoffRadiusSquare() const {return _cutoffRadiusSquare;}
	double getLJCutoffRadiusSquare() const {return _LJCutoffRadiusSquare;}
	void setCutoffRadiusSquare(const double c) {_cutoffRadiusSquare = c;}
	void setLJCutoffRadiusSquare(const double ljc) {_LJCutoffRadiusSquare = ljc;}


	/**
	 * called before the traversal starts.
	 *
	 * @param numCells number of cells in window
	 */
	virtual void initTraversal() = 0;

	/**
	 * Called before a cell is touched for the first time during an interation.
	 */
	virtual void preprocessCell(ParticleCell& cell) = 0;

	/**
	 * Called for each cell pair within the cutoff radius. Called exactly once per
	 * pair (i.e. pairs are not ordered).
	 *
	 * @note will not be called for empty cells.
	 * Sum up all macroscopic values (e.g. for hs) or only half of them (e.g. for fs)
         */
	virtual void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) = 0;

	/**
	 * Called when this cell is the current cell.
	 *
	 * @note will not be called for empty cells.
	 */
	virtual void processCell(ParticleCell& cell) = 0;

	virtual double processSingleMolecule(Molecule* m1, ParticleCell& cell2) = 0;

	/**
	 * Called after the cell has been considered for the last time during the traversal.
	 */
	virtual void postprocessCell(ParticleCell& cell) = 0;

	/**
	 * Called after the traversal finished.
	 */
	virtual void endTraversal() = 0;
};


#endif
