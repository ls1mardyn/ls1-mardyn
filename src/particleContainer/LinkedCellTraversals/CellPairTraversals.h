/*
 * CellPairTraversalWithDependencies.h
 *
 *  Created on: 15 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_CELLPAIRTRAVERSALS_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_CELLPAIRTRAVERSALS_H_

#include <vector>
#include <array>

class CellProcessor;

template <class CellTemplate>
class CellPairTraversals {
public:
	CellPairTraversals(
		std::vector<CellTemplate>& cells,
		const std::array<unsigned long, 3>& dims): _cells(&cells), _dims(dims) {}

	virtual ~CellPairTraversals() {}

	virtual void traverseCellPairs(CellProcessor& cellProcessor) = 0;
	virtual void traverseCellPairsOuter(CellProcessor& cellProcessor) = 0;
	virtual void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount) = 0;

protected:
	virtual void processBaseCell(CellProcessor& cellProcessor, unsigned long cellIndex) const = 0;

	//TODO:
	void traverseCellPairsNoDep(CellProcessor& cellProcessor);

	std::vector<CellTemplate> * _cells;
	std::array<unsigned long, 3> _dims;
};

#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_CELLPAIRTRAVERSALS_H_ */
