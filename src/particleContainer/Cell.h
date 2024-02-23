#ifndef CELL_H_
#define CELL_H_

class Cell {
public:
	Cell() : _cellIndex(0) {}

	virtual ~Cell() {}

	virtual void assignCellToHaloRegion() = 0;
	virtual void assignCellToBoundaryRegion() = 0;
	virtual void assignCellToInnerRegion() = 0;
	virtual void assignCellToInnerMostAndInnerRegion() = 0;

	virtual void skipCellFromHaloRegion() = 0;
	virtual void skipCellFromBoundaryRegion() = 0;
	virtual void skipCellFromInnerRegion() = 0;
	virtual void skipCellFromInnerMostRegion() = 0;

	virtual bool isHaloCell() const = 0;
	virtual bool isBoundaryCell() const = 0;
	virtual bool isInnerCell() const = 0;
	virtual bool isInnerMostCell() const = 0;

	virtual double getBoxMin(int d) const = 0;
	virtual double getBoxMax(int d) const = 0;
	virtual void setBoxMin(const double b[3]) = 0;
	virtual void setBoxMax(const double b[3]) = 0;

	unsigned long getCellIndex() const { return _cellIndex; }
	void setCellIndex(unsigned long cellIndex) { _cellIndex = cellIndex; }

protected:
	//! the index of a cell. On one process every index must be unique.
	unsigned long _cellIndex;
};

#endif /* CELL_H_ */
