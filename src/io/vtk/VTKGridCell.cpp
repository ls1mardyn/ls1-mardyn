/*
 * VTKGridCell.cpp
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#include "VTKGridCell.h"

VTKGridCell::VTKGridCell() { }

VTKGridCell::~VTKGridCell() { }

void VTKGridCell::setVertex(int index, VTKGridVertex* vertex) {
	_vertices[index] = vertex;
}


VTKGridVertex* const * VTKGridCell::getVertices() const {
	return _vertices;
}

void VTKGridCell::setIndex(int index) {
	_index = index;
}

unsigned int VTKGridCell::getIndex() const {
	return _index;
}

//bool VTKGridCell::isValidDataCell() const {
//	return _isValidDataCell;
//}


void VTKGridCell::setCellData(int numberOfMolecules) {
	_numberOfMolecules = numberOfMolecules;
}


int VTKGridCell::getNumberOfMolecules() {
	return _numberOfMolecules;
}
