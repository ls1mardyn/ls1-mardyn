/*
 * LinkedCellsHS.h
 *
 *  Created on: 27.04.2017
 *      Author: sauermann
 */

#ifndef LINKEDCELLSHS_H_
#define LINKEDCELLSHS_H_

#include "particleContainer/LinkedCells.h"

class LinkedCellsHS: public LinkedCells {
public:
	LinkedCellsHS(): LinkedCells() {};

	LinkedCellsHS(double bBoxMin[3], double bBoxMax[3], double cutoffRadius, double LJCutoffRadius, double cellsInCutoffRadius)
		: LinkedCells(bBoxMin, bBoxMax, cutoffRadius, LJCutoffRadius,	cellsInCutoffRadius) {};

	virtual ~LinkedCellsHS() {};


	virtual void traverseCell(long int cellIndex, CellProcessor& cellProcessor) override;

	virtual void traverseCellsC08(CellProcessor& cellProcessor) override;

	virtual inline bool requiresForceExchange() override{
		return true;
	}

};

#endif /* LINKEDCELLSHS_H_ */
