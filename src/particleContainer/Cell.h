#ifndef CELL_H_
#define CELL_H_

class Cell {
public:
	Cell() : haloCell(false), boundaryCell(false), innerCell(false) {}
	
	virtual ~Cell() {}

	void assignCellToHaloRegion() { haloCell = true; }
	void assignCellToBoundaryRegion() { boundaryCell = true; }
	void assignCellToInnerRegion() { innerCell = true; }
	
	void skipCellFromHaloRegion() { haloCell = false; }
	void skipCellFromBoundaryRegion() { boundaryCell = false; }
	void skipCellFromInnerRegion() { innerCell = false; }
	
	bool isHaloCell() const { return haloCell; }
	bool isBoundaryCell() const { return boundaryCell; }
	bool isInnerCell() const { return innerCell; }

protected:
	//! true when the cell is in the halo region 
	bool haloCell;
	//! true when the cell is in the boundary region
	bool boundaryCell;
	//! true when the cell is in the inner region
	bool innerCell;
};

#endif /* CELL_H_ */
