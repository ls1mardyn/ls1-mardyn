/*
 * VTKGridCell.cpp
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#include "VTKGridCell.h"

VTKGridCell::VTKGridCell() :
		_index(0), _numberOfMolecules(0), _load(0.), _level(0) {
}

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


void VTKGridCell::setCellData(int numberOfMolecules, double load, int level) {
	_numberOfMolecules = numberOfMolecules;
	_load = load;
	_level = level;
}


int VTKGridCell::getNumberOfMolecules() const {
	return _numberOfMolecules;
}


double VTKGridCell::getLoad() const {
	return _load;
}


int VTKGridCell::getLevel() const {
	return _level;
}
