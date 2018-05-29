#ifndef CELL_H_
#define CELL_H_

class Cell {
public:
	Cell() : haloCell(false), boundaryCell(false), innerCell(false), innerMostCell(false), _cellIndex(0) {
		for (int d = 0; d < 3; ++d) {
			_boxMin[d] = 0.0;
			_boxMax[d] = 0.0;
		}
	}
	
	virtual ~Cell() {}

	void assignCellToHaloRegion() { haloCell = true; }
	void assignCellToBoundaryRegion() { boundaryCell = true; }
	void assignCellToInnerRegion() { innerCell = true; }
	void assignCellToInnerMostAndInnerRegion() { innerCell = true; innerMostCell = true; }
	
	void skipCellFromHaloRegion() { haloCell = false; }
	void skipCellFromBoundaryRegion() { boundaryCell = false; }
	void skipCellFromInnerRegion() { innerCell = false; }
	void skipCellFromInnerMostRegion() { innerMostCell = false; }
	
	bool isHaloCell() const { return haloCell; }
	bool isBoundaryCell() const { return boundaryCell; }
	bool isInnerCell() const { return innerCell; }
	bool isInnerMostCell() const { return innerMostCell; }

	unsigned long getCellIndex() const { return _cellIndex; }
	void setCellIndex(unsigned long cellIndex) { _cellIndex = cellIndex; }

	double getBoxMin(int d) const { return _boxMin[d]; }
	double getBoxMax(int d) const { return _boxMax[d]; }
	void setBoxMin(const double b[3]) { for (int d = 0; d < 3; ++d) { _boxMin[d] = b[d]; } }
	void setBoxMax(const double b[3]) { for (int d = 0; d < 3; ++d) { _boxMax[d] = b[d]; } }

protected:
	//! true when the cell is in the halo region 
	bool haloCell;
	//! true when the cell is in the boundary region
	bool boundaryCell;
	//! true when the cell is in the inner region. Innermost cells are always also innerCells.
	bool innerCell;
	//! true when the cell is in the innermost region (does not have neighbors, that are boundary cells)
	bool innerMostCell;
	//! the index of a cell. On one process every index must be unique.
	unsigned long _cellIndex;
	//! lower left front corner
	double _boxMin[3];
	//! upper right back corner
	double _boxMax[3];
};

#endif /* CELL_H_ */
