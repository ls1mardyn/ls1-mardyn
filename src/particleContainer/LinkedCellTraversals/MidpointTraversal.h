/*
 * MidpointTraversal.h
 *
 *  Created on: 06.07.2017
 *      Author: sascha
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_MIDPOINTTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_MIDPOINTTRAVERSAL_H_

#include "particleContainer/LinkedCellTraversals/CellPairTraversals.h"
#include <tuple>

struct MidpointTraversalData : CellPairTraversalData {
};


template <class CellTemplate>
class MidpointTraversal : public CellPairTraversals<CellTemplate>{
public:
	MidpointTraversal(std::vector<CellTemplate>& cells,	const std::array<unsigned long, 3>& dims):
	CellPairTraversals<CellTemplate>(cells, dims) {
		computeOffsets();
	}
	virtual ~MidpointTraversal() {}

	/**
     * Reset all necessary data without reallocation.
     */
	virtual void rebuild(std::vector<CellTemplate> &cells,
		const std::array<unsigned long, 3> &dims, CellPairTraversalData *data) override;

	virtual void traverseCellPairs(CellProcessor& cellProcessor);
	virtual void traverseCellPairsOuter(CellProcessor& cellProcessor);
	virtual void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount);

	// Midpoint traversal requires force exchange
	virtual bool requiresForceExchange() const {return true;}

protected:
	virtual void processBaseCell(CellProcessor& cellProcessor, unsigned long cellIndex) const;

	// All pairs that have to be processed when calculating the forces
	std::array<std::pair<unsigned long, unsigned long>, 98> _cellPairOffsets;

private:
	void computeOffsets();

	void pairCells(std::tuple<long, long, long>& a,
			std::tuple<long, long, long>& b, int& index, std::array<long, 3>& dims);
	void pairOriginWithForewardNeighbors(int& index, std::array<long, 3>& dims);
	void pairCellsWithPlane(std::tuple<long, long, long>& cc,
			std::tuple<long, long, long>& oc, int& index, std::array<long, 3>& dims);
	void pairCellsWithAdjacentCorners(std::tuple<long, long, long>& ce,
			std::tuple<long, long, long>& oe, int& index, std::array<long, 3>& dims);

	// Offsets of all centers of the surrounding cube. Opposite sides are (i+3)%6
	std::array<std::tuple<long, long, long>, 6> _centers = {
			std::make_tuple( 1, 0, 0),
			std::make_tuple( 0, 1, 0),
			std::make_tuple( 0, 0, 1),
			// Mirrored at origin
			std::make_tuple(-1, 0, 0),
			std::make_tuple( 0,-1, 0),
			std::make_tuple( 0, 0,-1),
	};

	// Offsets of all edges of the surrounding cube. Opposite sides are (i+6)%12
	std::array<std::tuple<long, long, long>, 12> _edges = {
			std::make_tuple( 1, 1, 0),
			std::make_tuple( 1, 0, 1),
			std::make_tuple( 0, 1, 1),
			std::make_tuple(-1, 1, 0),
			std::make_tuple(-1, 0, 1),
			std::make_tuple( 0,-1, 1),
			// Mirrored at origin
			std::make_tuple(-1,-1, 0),
			std::make_tuple(-1, 0,-1),
			std::make_tuple( 0,-1,-1),
			std::make_tuple( 1,-1, 0),
			std::make_tuple( 1, 0,-1),
			std::make_tuple( 0, 1,-1),
	};

	// Offsets of all corners of the surrounding cube. Opposite sides are (i+4)%8
	std::array<std::tuple<long, long, long>, 8> _corners = {
			std::make_tuple( 1, 1, 1),
			std::make_tuple( 1, 1,-1),
			std::make_tuple( 1,-1, 1),
			std::make_tuple( 1,-1,-1),
			// Mirrored at origin
			std::make_tuple(-1,-1,-1),
			std::make_tuple(-1,-1, 1),
			std::make_tuple(-1, 1,-1),
			std::make_tuple(-1, 1, 1),
	};
};


template<class CellTemplate>
void MidpointTraversal<CellTemplate>::computeOffsets() {

	using threeDimensionalMapping::threeToOneD;
	using std::make_pair;

	// Dim array is int but we need it as long for some reason (copied from C08BasedTraversal)
	std::array<long, 3> dims;
	for (int d = 0; d < 3; ++d) {
		dims[d] = static_cast<long>(this->_dims[d]);
	}

	// Cell index of origin
	long int o = threeToOneD(0l, 0l, 0l, dims);

	// Next free index of the cellPairOffsets array
	int index = 0;

	// ----------------------------------------------------
	// Opposite cells with midpoint in origin
	// Opposite cells are set up to have index ((i+(n/2)) % n)
	// ----------------------------------------------------

	// process only half of the centers to get no pair twice
	for(int i=0; i<3; ++i){ // centers
		int j = (i+3)%6;
		pairCells(_centers[i], _centers[j], index, dims);
	}
	// process only half of the edges to get no pair twice
	for(int i=0; i<6; ++i){ // edges
		int j = (i+6)%12;
		pairCells(_edges[i], _edges[j], index, dims);
	}
	// process only half of the corners to get no pair twice
	for(int i=0; i<4; ++i){ // corners
		int j = (i+4)%8;
		pairCells(_corners[i], _corners[j], index, dims);
	}

	// ----------------------------------------------------
	// Forward neighbors of origin (similar to half shell)
	// ----------------------------------------------------

	pairOriginWithForewardNeighbors(index, dims);


	// ----------------------------------------------------
	// 'Cone' for all centers
	// ----------------------------------------------------
	for(int i=0; i<6; ++i){
		auto oc = _centers[(i+3)%6]; // opposite center
		auto cc = _centers[i]; // current center

		// Create pairs for each pair cc <--> oc + offset where offset is not in the direction of [CC Origin OC]
		// Therefore with every cell (except OC) in the plane with normal vector [CC Origin OC] that contains OC.
		pairCellsWithPlane(cc, oc, index, dims);
	}

	// ----------------------------------------------------
	// 'Cone' for all edges
	// ----------------------------------------------------
	for(int i=0; i<12; ++i){
		auto oe = _edges[(i+6)%12]; // opposite edge
		auto ce = _edges[i]; // current edge

		// Create pairs for each pair ce <--> (corners adjacent to oe)
		pairCellsWithAdjacentCorners(ce, oe, index, dims);

	}

}


template<class CellTemplate>
void MidpointTraversal<CellTemplate>::pairCells(std::tuple<long, long, long>& a,
		std::tuple<long, long, long>& b, int& index, std::array<long, 3>& dims){

	using threeDimensionalMapping::threeToOneD;

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

	// store offset pair
	_cellPairOffsets[index] = std::make_pair(aIndex, bIndex);
	// increment index for the next pair
	++index;
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::pairOriginWithForewardNeighbors(int& index, std::array<long, 3>& dims){
	//TODO ____ Implement
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::pairCellsWithPlane(std::tuple<long, long, long>& cc,
		std::tuple<long, long, long>& oc, int& index, std::array<long, 3>& dims){
	//TODO ____ Implement
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::pairCellsWithAdjacentCorners(std::tuple<long, long, long>& ce,
		std::tuple<long, long, long>& oe, int& index, std::array<long, 3>& dims){
	//TODO ____ Implement
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::rebuild(std::vector<CellTemplate> &cells,
		const std::array<unsigned long, 3> &dims, CellPairTraversalData *data){
	//TODO ____ Implement
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::traverseCellPairs(CellProcessor& cellProcessor){
	//TODO ____ Implement
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::traverseCellPairsOuter(CellProcessor& cellProcessor){
	//TODO ____ Implement
}

template<class CellTemplate>
void MidpointTraversal<CellTemplate>::traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount){
	//TODO ____ Implement
}


template<class CellTemplate>
void MidpointTraversal<CellTemplate>::processBaseCell(CellProcessor& cellProcessor, unsigned long cellIndex) const{
	//TODO ____ Implement
}

#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_MIDPOINTTRAVERSAL_H_ */
