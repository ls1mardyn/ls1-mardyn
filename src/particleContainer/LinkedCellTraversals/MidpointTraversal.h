/*
 * MidpointTraversal.h
 *
 *  Created on: 06.07.2017
 *      Author: sascha
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_MIDPOINTTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_MIDPOINTTRAVERSAL_H_

#include "particleContainer/LinkedCellTraversals/CellPairTraversals.h"
#include <array>

struct MidpointTraversalData : CellPairTraversalData {
};

template <class CellTemplate>
class MidpointTraversal : public CellPairTraversals<CellTemplate>{
public:
	MidpointTraversal(std::vector<CellTemplate>& cells,	const std::array<unsigned long, 3>& dims):
	CellPairTraversals<CellTemplate>(cells, dims) {
		fillArrays(); // Array initialization could be done directly but GCC 4.9 doesn't like it
		computeOffsets3D(); // >= C++14 this could be constexpr
		computeOffsets();
	}
	virtual ~MidpointTraversal() {}

	/**
     * Reset all necessary data without reallocation.
     */
	virtual void rebuild(std::vector<CellTemplate> &cells,
		const std::array<unsigned long, 3> &dims, double cellLength[3], double cutoff, CellPairTraversalData *data) override {
		CellPairTraversals<CellTemplate>::rebuild(cells, dims, cellLength, cutoff, data);
		computeOffsets();

		_innerMostCellIndices.clear();
		_notInnerMostCellIndices.clear();

		auto maxIndex = 1;
		for (auto d : dims)
			maxIndex *= d;

		for (auto i = 0; i < maxIndex; ++i) {
			if (this->_cells->at(i).isInnerMostCell()){
				_innerMostCellIndices.push_back(i);
			} else {
				_notInnerMostCellIndices.push_back(i);
			}
		}
	}

	virtual void traverseCellPairs(CellProcessor& cellProcessor) override;
	virtual void traverseCellPairsOuter(CellProcessor& cellProcessor) override;
	virtual void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount) override;

	// Midpoint traversal requires force exchange
	virtual bool requiresForceExchange() const override {return true;}
	// Midpoint supports cell sizes down to cutoff/2
	virtual unsigned maxCellsInCutoff() const override {return 2;}

protected:
	virtual void processBaseCell(CellProcessor& cellProcessor, unsigned long baseIndex) const;

	// All pairs that have to be processed when calculating the forces (excluding self)
	std::array<std::pair<long, long>, 62> _cellPairOffsets;
	std::array<std::pair<std::array<long, 3>, std::array<long, 3>>, 62> _offsets3D;


private:

	void computeOffsets();
	void computeOffsets3D(); // This could be changed to constexpr in C++14
	void fillArrays();

	void pairOriginWithForwardNeighbors(int& index);
	void pairCellsWithPlane(std::array<long, 3>& cc,
			std::array<long, 3>& oc, int& index);
	void pairCellsWithAdjacentCorners(std::array<long, 3>& ce,
			std::array<long, 3>& oe, int& index);

	// Offsets of all faces of the surrounding cube. Opposite sides are (i+3)%6
	std::array<std::array<long, 3>, 6> _faces;

	// Offsets of all edges of the surrounding cube. Opposite sides are (i+6)%12
	std::array<std::array<long, 3>, 12> _edges;

	// Offsets of all corners of the surrounding cube. Opposite sides are (i+4)%8
	std::array<std::array<long, 3>, 8> _corners;

	std::vector<unsigned long> _innerMostCellIndices; //!< Vector containing the indices (for the cells vector) of all inner cells (without boundary)
	std::vector<unsigned long> _notInnerMostCellIndices; //!< Vector containing the indices (for the cells vector) of all outer cells
};

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::fillArrays() {
	// Faces
	_faces[0] = { 1, 0, 0};
	_faces[1] = { 0, 1, 0};
	_faces[2] = { 0, 0, 1};
	// Mirrored at origin
	_faces[3] = {-1, 0, 0};
	_faces[4] = { 0,-1, 0};
	_faces[5] = { 0, 0,-1};

	// Edges
	_edges[0] = { 1, 1, 0};
	_edges[1] = { 1, 0, 1};
	_edges[2] = { 0, 1, 1};
	_edges[3] = {-1, 1, 0};
	_edges[4] = {-1, 0, 1};
	_edges[5] = { 0,-1, 1};
	// Mirrored at origin
	_edges[6] = {-1,-1, 0};
	_edges[7] = {-1, 0,-1};
	_edges[8] = { 0,-1,-1};
	_edges[9] = { 1,-1, 0};
	_edges[10] = { 1, 0,-1};
	_edges[11] = { 0, 1,-1};

	// Corners
	_corners[0] = { 1, 1, 1};
	_corners[1] = { 1, 1,-1};
	_corners[2] = { 1,-1, 1};
	_corners[3] = { 1,-1,-1};
	// Mirrored at origin
	_corners[4] = {-1,-1,-1};
	_corners[5] = {-1,-1, 1};
	_corners[6] = {-1, 1,-1};
	_corners[7] = {-1, 1, 1};
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::computeOffsets3D() {

	// Next free index of the offsets array
	int index = 0;

	// ----------------------------------------------------
	// Opposite cells with midpoint in origin
	// Opposite cells are set up to have index ((i+(n/2)) % n)
	// ----------------------------------------------------

	// process only half of the faces to get no pair twice
	for(int i=0; i<3; ++i){ // faces
		int j = (i+3)%6;
		_offsets3D[index++] = make_pair(_faces[i], _faces[j]);
	}
	// process only half of the edges to get no pair twice
	for(int i=0; i<6; ++i){ // edges
		int j = (i+6)%12;
		_offsets3D[index++] = make_pair(_edges[i], _edges[j]);
	}
	// process only half of the corners to get no pair twice
	for(int i=0; i<4; ++i){ // corners
		int j = (i+4)%8;
		_offsets3D[index++] = make_pair(_corners[i], _corners[j]);
	}

	mardyn_assert(index == 13);

	// ----------------------------------------------------
	// Forward neighbors of origin (similar to half shell)
	// ----------------------------------------------------

	pairOriginWithForwardNeighbors(index);

	mardyn_assert(index == 13+13);


	// ----------------------------------------------------
	// 'Cone' for half the faces
	// ----------------------------------------------------
	for(int i=0; i<3; ++i){
		auto oc = _faces[(i+3)%6]; // opposite center
		auto cc = _faces[i]; // current center

		// Create pairs for each pair cc <--> oc + offset where offset is not in the direction of [CC Origin OC]
		// Therefore with every cell (except OC) in the plane with normal vector [CC Origin OC] that contains OC.
		pairCellsWithPlane(cc, oc, index);
	}

	mardyn_assert(index == 13+13+24);

	// ----------------------------------------------------
	// 'Cone' for half the edges
	// ----------------------------------------------------
	for(int i=0; i<6; ++i){
		auto oe = _edges[(i+6)%12]; // opposite edge
		auto ce = _edges[i]; // current edge

		// Create pairs for each pair ce <--> (corners adjacent to oe)
		pairCellsWithAdjacentCorners(ce, oe, index);

	}

	// We need exactly 62 cell offset pairs
	mardyn_assert(index == 13+13+24+12);
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::computeOffsets() {
	using threeDimensionalMapping::threeToOneD;

	// Dim array is int but we need it as long for some reason (copied from C08BasedTraversal)
	std::array<long, 3> dims;
	for (int d = 0; d < 3; ++d) {
		dims[d] = static_cast<long>(this->_dims[d]);
	}

	for(unsigned int i=0; i<_offsets3D.size(); ++i){

		auto a = _offsets3D[i].first;
		auto b = _offsets3D[i].second;

		auto ax = std::get<0>(a);
		auto ay = std::get<1>(a);
		auto az = std::get<2>(a);

		auto bx = std::get<0>(b);
		auto by = std::get<1>(b);
		auto bz = std::get<2>(b);

		mardyn_assert((abs(ax) <= 1) && (abs(ay) <= 1) && (abs(az) <= 1));
		mardyn_assert((abs(bx) <= 1) && (abs(by) <= 1) && (abs(bz) <= 1));

		// convert 3d index to 1d
		auto aIndex = threeToOneD(ax, ay, az, dims);
		auto bIndex = threeToOneD(bx, by, bz, dims);

		auto offsetPair = std::make_pair(aIndex, bIndex);

		// store offset pair
		_cellPairOffsets[i] = offsetPair;
	}

}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::pairOriginWithForwardNeighbors(int& index){
	using std::make_pair;

	std::array<long, 3> origin = {0l, 0l, 0l};

	for(long y=-1; y<=1; ++y){ // 3 * 4
		for(long x=-1; x<=1; ++x){ // 3
			_offsets3D[index++] = make_pair(origin, std::array<long,3>{x, y, 1l});
		}
		// 1st
		_offsets3D[index++] = make_pair(origin, std::array<long,3>{1l, y, 0l});
	}

	// 13th
	_offsets3D[index++] = make_pair(origin, std::array<long,3>{0l, 1l, 0l});

}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::pairCellsWithPlane(std::array<long, 3>& cc,
		std::array<long, 3>& oc, int& index){
	// Pairs the current cell (cc) with every cell next to oc except the origin.
	for(int i=-1; i<=1; ++i){
		for(int j=-1; j<=1; ++j){

			// Skip pair cc <--> oc
			if(i==0 && j==0) continue;

			std::array<long, 3> cell;

			if(std::get<0>(oc) == 0){
				if(std::get<1>(oc) == 0){ // x and y are 0
					// Center in +-z direction -> ring around oc in xy plane
					cell = {i, j, std::get<2>(oc)};
				} else { // x and z are 0
					// Center in +-y direction -> ring around oc in xz plane
					cell = {i, std::get<1>(oc), j};
				}
			} else { // y and z are 0
				// Center in +-x direction -> ring around oc in yz plane
				cell = {std::get<0>(oc), i, j};
			}

			_offsets3D[index++] = std::make_pair(cc, cell);
		}
	}
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::pairCellsWithAdjacentCorners(std::array<long, 3>& ce,
		std::array<long, 3>& oe, int& index){
	// Pairs the current edge cell (ce) with both adjacent corners to oe
	for(int i=-1; i<=1; i+=2){ // i=-1 and i=1

			std::array<long, 3> cell;

			if(std::get<0>(oe) == 0) { // x = 0
				// Corners are oe with x=+-1
				cell = {i, std::get<1>(oe), std::get<2>(oe)};
			} else if(std::get<1>(oe) == 0) { // y = 0
				// Corners are oe with y=+-1
				cell = {std::get<0>(oe), i, std::get<2>(oe)};
			} else { // z = 0
				// Corners are oe with z=+-1
				cell = {std::get<0>(oe), std::get<1>(oe), i};
			}

			_offsets3D[index++] = std::make_pair(ce, cell);
	}
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::traverseCellPairs(CellProcessor& cellProcessor){
	unsigned long start = 0ul;
	unsigned long end = this->_cells->size();
	for(auto i = start; i<end; ++i){
		processBaseCell(cellProcessor, i);
	}
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::traverseCellPairsOuter(CellProcessor& cellProcessor){
	//TODO ____ Test
	for(unsigned long outerCellIndex : _notInnerMostCellIndices){
		processBaseCell(cellProcessor, outerCellIndex);
	}
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount){
	//TODO ____ Test
	unsigned long start =  _innerMostCellIndices.size() * stage / stageCount;
	unsigned long end =  _innerMostCellIndices.size() * (stage+1) / stageCount;
	for (unsigned long i = start; i < end; ++i) {
		processBaseCell(cellProcessor, _innerMostCellIndices.at(i));
	}
}


template<class CellTemplate>
void MidpointTraversal<CellTemplate>::processBaseCell(CellProcessor& cellProcessor, unsigned long baseIndex) const{

	unsigned long maxIndex = this->_cells->size() - 1;

	CellTemplate& baseCell = this->_cells->at(baseIndex);

	if(!baseCell.isHaloCell()) {

		// Process all cell pairs for this cell
		for(auto& current_pair : _cellPairOffsets){

			unsigned long offset1 = current_pair.first;
			unsigned long cellIndex1 = baseIndex + offset1;

			unsigned long offset2 = current_pair.second;
			unsigned long cellIndex2 = baseIndex + offset2;

			CellTemplate& cell1 = this->_cells->at(cellIndex1);
			CellTemplate& cell2 = this->_cells->at(cellIndex2);

			const bool sumAllMacroscopic = true;
			cellProcessor.processCellPair(cell1, cell2, sumAllMacroscopic);

		}

		// Process base cell itself
		cellProcessor.processCell(baseCell);

	}
}

#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_MIDPOINTTRAVERSAL_H_ */
