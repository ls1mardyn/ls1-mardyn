/*
 * DrawableMolecule.cpp
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt
 */

#include "DrawableMolecule.h"

DrawableMolecule::DrawableMolecule(const Molecule& molecule, int numComponents)
: _x(molecule.r(0),molecule.r(1), molecule.r(2)),
  _v(molecule.v(0), molecule.v(1), molecule.v(2)), _id(molecule.getID()), _numComponents(numComponents), _cid(molecule.componentid()) {
}

DrawableMolecule::DrawableMolecule()
: _x(0, 0, 0),
  _v(0, 0, 0), _id(0), _numComponents(1), _cid(0) {
}

DrawableMolecule::~DrawableMolecule() {
}


std::vector<std::string> DrawableMolecule::getDrawableValues() const {
	std::vector<std::string> v;
	v.push_back("Molecule");
	v.push_back("Molecule ID");
	v.push_back("Component ID");
	v.push_back("Velocity");
	return v;
}

vtkSmartPointer<vtkActor> DrawableMolecule::draw(std::string valueName){
	if (valueName == "Molecule") {
		return drawValue(_x,_id, 0, _numObjects, false);
	}
	else if (valueName == "Molecule ID"){
		return drawValue(_x, _id, 0, _numObjects, true);
	}
	else if (valueName == "Component ID"){
		return drawValue(_x, _cid, 0, _numComponents, true);
	}
	else if (valueName == "Velocity" ){
		return drawVector(_x,_v);
	}
	else {
		std::cout<<"Invalid value name"<<std::endl;
		return 0;
	}
}

