/*
 * VTKGridVertex.cpp
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#include "VTKGridVertex.h"

VTKGridVertex::VTKGridVertex()
: _index(-1) { }


VTKGridVertex::~VTKGridVertex() {
}


const double* const VTKGridVertex::getCoordinates() const {
	return _coordinates;
}


void VTKGridVertex::setIndex(int index) {
	_index = index;
}


int VTKGridVertex::getIndex() const {
	return _index;
}


void VTKGridVertex::setCoordinates(double x, double y, double z) {
	_coordinates[0] = x;
	_coordinates[1] = y;
	_coordinates[2] = z;
}
