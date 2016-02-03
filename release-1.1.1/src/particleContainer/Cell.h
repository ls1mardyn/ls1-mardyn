#ifndef CELL_H_
#define CELL_H_

class Cell {
public:
	Cell() : haloCell(false), boundaryCell(false), innerCell(false), inActiveWindow(false) {}
	
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

	void setInActiveWindow() { inActiveWindow = true;}
	void clearInActiveWindow() { inActiveWindow = false;}
	bool isInActiveWindow() { return inActiveWindow; }

protected:
	//! true when the cell is in the halo region 
	bool haloCell;
	//! true when the cell is in the boundary region
	bool boundaryCell;
	//! true when the cell is in the inner region
	bool innerCell;

	//! for debugging to check that all cells are processed by the cell handler
	bool inActiveWindow;
};

#endif /* CELL_H_ */
