#ifndef CELL_H_
#define CELL_H_

class Cell {
public:
	Cell() : haloCell(false), boundaryCell(false), innerCell(false), innerMostCell(false) {}
	
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

protected:
	//! true when the cell is in the halo region 
	bool haloCell;
	//! true when the cell is in the boundary region
	bool boundaryCell;
	//! true when the cell is in the inner region. Innermost cells are always also innerCells.
	bool innerCell;
	//! true when the cell is in the innermost region (does not have neighbors, that are boundary cells)
	bool innerMostCell;
};

#endif /* CELL_H_ */
