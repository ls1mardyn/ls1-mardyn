/*
 * VTKGridVertex.h
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#ifndef VTKGRIDVERTEX_H_
#define VTKGRIDVERTEX_H_

/**
 * Represents a vertex of an unstructured vtk grid.
 */
class VTKGridVertex {

private:
	//! the coordinates of this vertex
	double _coordinates[3];

	/*!
	 * @short index of this vertex
	 * When the vertex first is plotted, it is assigned an index. Later on, when
	 * adjacent cells are plotted, either the vertex is plotted again, too, (thus
	 * we would end up with 8 vertices with exactly the same coordinates plotted) or
	 * it is referenced via the index.
	 */
	int _index;

public:

	VTKGridVertex();

	virtual ~VTKGridVertex();

	const double* const getCoordinates() const;

	void setIndex(int index);

	int getIndex() const;

	void setCoordinates(double x, double y, double z);
};

#endif /* VTKGRIDVERTEX_H_ */
