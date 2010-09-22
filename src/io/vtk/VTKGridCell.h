/*
 * VTKGridCell.h
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#ifndef VTKGRIDCELL_H_
#define VTKGRIDCELL_H_

class VTKGridVertex;

/**
 * Represents a grid cell in a vtk unstructured grid.
 */
class VTKGridCell {

private:

	/**
	 * The 8 vertices by which a cell is defined
	 */
	VTKGridVertex* _vertices[8];

	//! the index of the cell (should correspond to the index in the linked-cell data structure).
	unsigned int _index;

	int _numberOfMolecules;

	/*! indicates if this cell contains valid data, which should be plotted.
	 *  This is probably not true for halo cells. Nevertheless it might be convenient
	 *  to use also VTKGridCells for each halo cell to get rid of nasty index calculations.
	 */
//	bool _isValidDataCell;

public:

	VTKGridCell();

	virtual ~VTKGridCell();

	void setVertex(int vertexIndex, VTKGridVertex* vertex);

	VTKGridVertex* const * getVertices() const;

	void setIndex(int index);

	unsigned int getIndex() const;

//	bool isValidDataCell() const;

	/**
	 * set all the data fields.
	 */
	void setCellData(int numberOfMolecules);

	int getNumberOfMolecules();
};

#endif /* VTKGRIDCELL_H_ */
