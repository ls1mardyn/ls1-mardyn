/*
 * C08BasedTraversals.h
 *
 *  Created on: 16 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08BASEDTRAVERSALS_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08BASEDTRAVERSALS_H_

#include "particleContainer/LinkedCellTraversals/CellPairTraversals.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "utils/threeDimensionalMapping.h"
#include "WrapOpenMP.h"

template <class CellTemplate>
class C08BasedTraversals : public CellPairTraversals<CellTemplate> {
public:
	C08BasedTraversals(
			std::vector<CellTemplate>& cells,
			const std::array<unsigned long, 3>& dims) :
			CellPairTraversals<CellTemplate>(cells, dims) {
		computeOffsets();
	}
	virtual ~C08BasedTraversals() {
	}

	virtual void rebuild(std::vector<CellTemplate> &cells,
						 const std::array<unsigned long, 3> &dims,
						 CellPairTraversalData *data) {
		CellPairTraversals<CellTemplate>::rebuild(cells, dims, data);
		computeOffsets();
	};

protected:
	void processBaseCell(CellProcessor& cellProcessor, unsigned long cellIndex) const;
	std::array<std::pair<unsigned long, unsigned long>, 14> _cellPairOffsets;
	std::array<unsigned long, 8> _cellOffsets;

private:
	void computeOffsets();
};

template<class CellTemplate>
void C08BasedTraversals<CellTemplate>::processBaseCell(
		CellProcessor& cellProcessor, unsigned long baseIndex) const {

#ifdef MARDYN_WR
	for (int i = 0; i < 8; ++i) {
		this->_cells->at(baseIndex + _cellOffsets[i]).prefetch();
	}
#endif /* MARDYN_WR */

	using std::pair;

#ifndef NDEBUG
	// map to 3D index and check that we're not on the "right" boundary
	std::array<unsigned long, 3> threeDIndex = threeDimensionalMapping::oneToThreeD(baseIndex, this->_dims);
	for (int d = 0; d < 3; ++d) {
		mardyn_assert(threeDIndex[d] != this->_dims[d]-1);
	}
#endif

	const int num_pairs = _cellPairOffsets.size();
	for(int j = 0; j < num_pairs; ++j) {
		pair<long, long> current_pair = _cellPairOffsets[j];

		unsigned offset1 = current_pair.first;
		unsigned cellIndex1 = baseIndex + offset1;

		unsigned offset2 = current_pair.second;
		unsigned cellIndex2 = baseIndex + offset2;

		CellTemplate& cell1 = this->_cells->at(cellIndex1);
		CellTemplate& cell2 = this->_cells->at(cellIndex2);

		if(cell1.isHaloCell() and cell2.isHaloCell()) {
			continue;
		}

		if(cellIndex1 == cellIndex2) {
			cellProcessor.processCell(cell1);
		}
		else {
			if(!cell1.isHaloCell()) {
				cellProcessor.processCellPair(cell1, cell2);
			}
			else {
				cellProcessor.processCellPair(cell2, cell1);
			}
		}
	}
}

template<class CellTemplate>
void C08BasedTraversals<CellTemplate>::computeOffsets() {
	using threeDimensionalMapping::threeToOneD;
	using std::make_pair;

	std::array<long, 3> dims;
	for (int d = 0; d < 3; ++d) {
		dims[d] = static_cast<long>(this->_dims[d]);
	}

	long int o   = threeToOneD(0l, 0l, 0l, dims); // origin
	long int x   = threeToOneD(1l, 0l, 0l, dims); // displacement to the right
	long int y   = threeToOneD(0l, 1l, 0l, dims); // displacement ...
	long int z   = threeToOneD(0l, 0l, 1l, dims);
	long int xy  = threeToOneD(1l, 1l, 0l, dims);
	long int yz  = threeToOneD(0l, 1l, 1l, dims);
	long int xz  = threeToOneD(1l, 0l, 1l, dims);
	long int xyz = threeToOneD(1l, 1l, 1l, dims);

	int i = 0;
	// if incrementing along X, the following order will be more cache-efficient:
	_cellPairOffsets[i++] = make_pair(o, o  );
	_cellPairOffsets[i++] = make_pair(o, y  );
	_cellPairOffsets[i++] = make_pair(y, z  );
	_cellPairOffsets[i++] = make_pair(o, z  );
	_cellPairOffsets[i++] = make_pair(o, yz );

	_cellPairOffsets[i++] = make_pair(x, yz );
	_cellPairOffsets[i++] = make_pair(x, y  );
	_cellPairOffsets[i++] = make_pair(x, z  );
	_cellPairOffsets[i++] = make_pair(o, x  );
	_cellPairOffsets[i++] = make_pair(o, xy );
	_cellPairOffsets[i++] = make_pair(xy, z );
	_cellPairOffsets[i++] = make_pair(y, xz );
	_cellPairOffsets[i++] = make_pair(o, xz );
	_cellPairOffsets[i++] = make_pair(o, xyz);

	i = 0;
	_cellOffsets[i++] =   o;
	_cellOffsets[i++] =   y;
	_cellOffsets[i++] =   z;
	_cellOffsets[i++] =  yz;

	_cellOffsets[i++] =   x;
	_cellOffsets[i++] =  xy;
	_cellOffsets[i++] =  xz;
	_cellOffsets[i++] = xyz;

}


#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08BASEDTRAVERSALS_H_ */
