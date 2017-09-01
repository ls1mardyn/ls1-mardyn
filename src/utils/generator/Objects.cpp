#include "utils/generator/Objects.h"
#include "utils/generator/ObjectFactory.h"

#include "utils/Logger.h"

using namespace std;
using Log::global_log;
/** Class implemenationso ofcuboid */

Cuboid::Cuboid::Cuboid() : _lowerCorner{0, 0, 0}, _upperCorner{0, 0, 0} {}
Cuboid::Cuboid(double lower[3], double upper[3]) {
	for(int d = 0; d < 3; d++) {
		_lowerCorner[d] = lower[d];
		_upperCorner[d] = upper[d];
	}
}

void Cuboid::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValueReduced("lower/x", _lowerCorner[0]);
	xmlconfig.getNodeValueReduced("lower/y", _lowerCorner[1]);
	xmlconfig.getNodeValueReduced("lower/z", _lowerCorner[2]);
	global_log->info() << "lower corner: " << _lowerCorner[0] << ", " << _lowerCorner[1] << ", " << _lowerCorner[2] << endl;
	xmlconfig.getNodeValueReduced("upper]/x", _upperCorner[0]);
	xmlconfig.getNodeValueReduced("upper]/y", _upperCorner[1]);
	xmlconfig.getNodeValueReduced("upper]/z", _upperCorner[2]);
	global_log->info() << "upper corner: " << _upperCorner[0] << ", " << _upperCorner[1] << ", " << _upperCorner[2] << endl;
}

bool Cuboid::isInside(double r[3]) {
	return (_lowerCorner[0] <= r[0] && r[0] <= _upperCorner[0])
		&& (_lowerCorner[1] <= r[1] && r[1] <= _upperCorner[1])
		&& (_lowerCorner[2] <= r[2] && r[2] <= _upperCorner[2]);
}

bool Cuboid::isInsideNoBorder(double r[3]) {
	return (_lowerCorner[0] < r[0] && r[0] < _upperCorner[0])
		&& (_lowerCorner[1] < r[1] && r[1] < _upperCorner[1])
		&& (_lowerCorner[2] < r[2] && r[2] < _upperCorner[2]);
}

void Cuboid::getBboxMin(double rmin[3]) {
	for(int d = 0; d < 3; d++) {
		rmin[d] = _lowerCorner[d];
	}
}
void Cuboid::getBboxMax(double rmax[3]) {
	for(int d = 0; d < 3; d++) {
		rmax[d] = _upperCorner[d];
	}
}

/** Class implementing a sphere */

Sphere::Sphere() : _center{0, 0, 0}, _radius(0), _radiusSquare(_radius*_radius) {}
Sphere::Sphere(double center[3], double r) : _radius(r), _radiusSquare(r*r) {
	for(int d = 0; d < 3; d++) {
		_center[d] = center[d];
	}
}

void Sphere::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValueReduced("center/x", _center[0]);
	xmlconfig.getNodeValueReduced("center/y", _center[1]);
	xmlconfig.getNodeValueReduced("center/z", _center[2]);
	global_log->info() << "center coordinate: "<< _center[0] << ", " << _center[1] << ", " << _center[2] << endl;
	xmlconfig.getNodeValueReduced("radius", _radius);
	_radiusSquare = _radius*_radius;
	global_log->info() << "Radius: " << _radius << endl;
}


bool Sphere::isInside(double r[3]) {
	double dr2 = 0.0;
	double dr[3];
	for(int d = 0; d < 3; d++) {
		dr[d] = r[d] - _center[d];
		dr2 += dr[d] * dr[d];
	}
	return (dr2 <= _radiusSquare);
}

bool Sphere::isInsideNoBorder(double r[3]) {
	double dr2 = 0.0;
	double dr[3];
	for(int d = 0; d < 3; d++) {
		dr[d] = r[d] - _center[d];
		dr2 += dr[d] * dr[d];
	}
	return (dr2 < _radiusSquare);
}

void Sphere::getBboxMin(double rmin[3]) {
	for(int d = 0; d < 3; d++) {
		rmin[d] = _center[d] - _radius;
	}
}
void Sphere::getBboxMax(double rmax[3]) {
	for(int d = 0; d < 3; d++) {
		rmax[d] = _center[d] + _radius;
	}
}

Cylinder::Cylinder::Cylinder() : _radius(0), _height(0), _centerBase{0, 0, 0}, _radiusSquare(_radius*_radius)  {}
Cylinder::Cylinder(double centerBase[3], double radius, double height) : _radius(radius), _height(height), _radiusSquare(radius*radius) {
	for(int d = 0; d < 3; d++) {
		_centerBase[d] = centerBase[d];
	}
}

void Cylinder::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValueReduced("centerBase/x", _centerBase[0]);
	xmlconfig.getNodeValueReduced("centerBase/y", _centerBase[1]);
	xmlconfig.getNodeValueReduced("centerBase/z", _centerBase[2]);
	global_log->info() << "center base coordinate: "<< _centerBase[0] << ", " << _centerBase[1] << ", " << _centerBase[2] << endl;
	xmlconfig.getNodeValueReduced("radius", _radius);
	_radiusSquare = _radius*_radius;
	global_log->info() << "Radius: " << _radius << endl;
	xmlconfig.getNodeValueReduced("height", _height);
	global_log->info() << "Height: " << _height << endl;
}

bool Cylinder::isInside(double r[3]) {
	double dr[2];
	dr[0] = r[0] - _centerBase[0];
	dr[1] = r[1] - _centerBase[1];
	return (dr[0]*dr[0] + dr[1]*dr[1] <= _radiusSquare) && (r[2] >= _centerBase[2]) && (r[2] <= _centerBase[2] + _height);
}

bool Cylinder::isInsideNoBorder(double r[3]) {
	double dr[2];
	dr[0] = r[0] - _centerBase[0];
	dr[1] = r[1] - _centerBase[1];
	return (dr[0]*dr[0] + dr[1]*dr[1] < _radiusSquare) && (r[2] > _centerBase[2]) && (r[2] < _centerBase[2] + _height);
}

void Cylinder::getBboxMin(double rmin[3]) {
	rmin[0] = _centerBase[0] - _radius;
	rmin[1] = _centerBase[1] - _radius;
	rmin[2] = _centerBase[2];
}
void Cylinder::getBboxMax(double rmax[3]) {
	rmax[0] = _centerBase[0] + _radius;
	rmax[1] = _centerBase[1] + _radius;
	rmax[2] = _centerBase[2] + _height;
}


ObjectUnification::ObjectUnification() : _ob1(nullptr), _ob2(nullptr) {}

void ObjectUnification::readXML(XMLfileUnits& xmlconfig) {
	ObjectFactory object_factory;
	if(xmlconfig.changecurrentnode("object1")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob1 = object_factory.create(object_type);
		_ob1->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("object2")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob2 = object_factory.create(object_type);
		_ob2->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
}


ObjectSubtractor::ObjectSubtractor() : _ob1(nullptr), _ob2(nullptr) {}

void ObjectSubtractor::readXML(XMLfileUnits& xmlconfig) {
	ObjectFactory object_factory;
	if(xmlconfig.changecurrentnode("object1")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob1 = object_factory.create(object_type);
		_ob1->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("object2")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob2 = object_factory.create(object_type);
		_ob2->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
}


ObjectIntersection::ObjectIntersection() : _ob1(nullptr), _ob2(nullptr) {}

void ObjectIntersection::readXML(XMLfileUnits& xmlconfig) {
	ObjectFactory object_factory;
	if(xmlconfig.changecurrentnode("object1")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob1 = object_factory.create(object_type);
		_ob1->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("object2")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob2 = object_factory.create(object_type);
		_ob2->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
}
