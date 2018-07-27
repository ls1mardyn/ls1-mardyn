#include "Coordinate3D.h"

Coordinate3D::Coordinate3D(XMLfileUnits& xmlconfig, std::string nodename) : _vec{0.0, 0.0, 0.0} {
	std::string oldpath = xmlconfig.getcurrentnodepath();
	xmlconfig.changecurrentnode(nodename);
	readXML(xmlconfig);
	xmlconfig.changecurrentnode(oldpath);
}


void Coordinate3D::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValueReduced("x", _vec[0]);
	xmlconfig.getNodeValueReduced("y", _vec[1]);
	xmlconfig.getNodeValueReduced("z", _vec[2]);
}

void Coordinate3D::get(double vec[3]) {
	for(int d = 0; d < 3; ++d) {
		vec[d] = _vec[d];
	}
}

