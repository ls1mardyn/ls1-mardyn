#include "utils/generator/Objects.h"
#include "utils/generator/ObjectFactory.h"

#include "utils/Coordinate3D.h"
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
	Coordinate3D lowerCorner(xmlconfig, "lower");
	Coordinate3D upperCorner(xmlconfig, "upper");
	lowerCorner.get(_lowerCorner);
	upperCorner.get(_upperCorner);
	global_log->info() << "lower corner: " << _lowerCorner[0] << ", " << _lowerCorner[1] << ", " << _lowerCorner[2] << endl;
	global_log->info() << "upper corner: " << _upperCorner[0] << ", " << _upperCorner[1] << ", " << _upperCorner[2] << endl;
}

bool Cuboid::isInside(double r[3]) {
	return (lowerCorner(0) <= r[0] && r[0] <= upperCorner(0))
		&& (lowerCorner(1) <= r[1] && r[1] <= upperCorner(1))
		&& (lowerCorner(2) <= r[2] && r[2] <= upperCorner(2));
}

bool Cuboid::isInsideNoBorder(double r[3]) {
	return (lowerCorner(0) < r[0] && r[0] < upperCorner(0))
		&& (lowerCorner(1) < r[1] && r[1] < upperCorner(1))
		&& (lowerCorner(2) < r[2] && r[2] < upperCorner(2));
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
	Coordinate3D center(xmlconfig, "center");
	center.get(_center);
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
	Coordinate3D centerBase(xmlconfig, "centerBase");
	centerBase.get(_centerBase);
	global_log->info() << "center base coordinate: "<< _centerBase[0] << ", " << _centerBase[1] << ", " << _centerBase[2] << endl;
	xmlconfig.getNodeValueReduced("radius", _radius);
	_radiusSquare = _radius*_radius;
	global_log->info() << "Radius: " << _radius << endl;
	xmlconfig.getNodeValueReduced("height", _height);
	global_log->info() << "Height: " << _height << endl;
	int axis = 2;
	xmlconfig.getNodeValue("direction", axis);
	_direction.axis = static_cast<ObjectAxis>(axis);
	global_log->info() << "Direction: " << _direction.axis << endl;
	switch(_direction.axis) {
		case OBJ_AXIS_X:
			_direction.base1 = 1; _direction.base2 = 2; _direction.height = 0; break;
		case OBJ_AXIS_Y:
			_direction.base1 = 0; _direction.base2 = 2; _direction.height = 1; break;
		case OBJ_AXIS_Z:
		default:
			_direction.base1 = 0; _direction.base2 = 1; _direction.height = 2;
	}
}

bool Cylinder::isInside(double r[3]) {
	double dr[2];
	uint8_t b1, b2, h;
	b1 = _direction.base1; b2 = _direction.base2; h = _direction.height;
	dr[0] = r[b1] - _centerBase[b1];
	dr[1] = r[b2] - _centerBase[b2];
	return (dr[0]*dr[0] + dr[1]*dr[1] <= _radiusSquare) && (r[h] >= _centerBase[h]) && (r[h] <= _centerBase[h] + _height);
}

bool Cylinder::isInsideNoBorder(double r[3]) {
	double dr[2];
	uint8_t b1, b2, h;
	b1 = _direction.base1; b2 = _direction.base2; h = _direction.height;
	dr[0] = r[b1] - _centerBase[b1];
	dr[1] = r[b2] - _centerBase[b2];
	return (dr[0]*dr[0] + dr[1]*dr[1] < _radiusSquare) && (r[h] > _centerBase[h]) && (r[h] < _centerBase[h] + _height);
}

void Cylinder::getBboxMin(double rmin[3]) {
	uint8_t b1, b2, h;
	b1 = _direction.base1; b2 = _direction.base2; h = _direction.height;
	rmin[b1] = _centerBase[b1] - _radius;
	rmin[b2] = _centerBase[b2] - _radius;
	rmin[h] = _centerBase[h];
}
void Cylinder::getBboxMax(double rmax[3]) {
	uint8_t b1, b2, h;
	b1 = _direction.base1; b2 = _direction.base2; h = _direction.height;
	rmax[b1] = _centerBase[b1] + _radius;
	rmax[b2] = _centerBase[b2] + _radius;
	rmax[h] = _centerBase[h] + _height;
}


ObjectUnification::ObjectUnification() : _ob1(nullptr), _ob2(nullptr) {}

void ObjectUnification::readXML(XMLfileUnits& xmlconfig) {
	ObjectFactory object_factory;
	if(xmlconfig.changecurrentnode("object1")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob1 = std::shared_ptr<Object>(object_factory.create(object_type));
		_ob1->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("object2")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob2 = std::shared_ptr<Object>(object_factory.create(object_type));
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
		_ob1 = std::shared_ptr<Object>(object_factory.create(object_type));
		_ob1->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("object2")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob2 = std::shared_ptr<Object>(object_factory.create(object_type));
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
		_ob1 = std::shared_ptr<Object>(object_factory.create(object_type));
		_ob1->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("object2")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_ob2 = std::shared_ptr<Object>(object_factory.create(object_type));
		_ob2->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
}


ObjectShifter::ObjectShifter() : _obj(nullptr) {
	_shift[0] = 0.0;
	_shift[1] = 0.0;
	_shift[2] = 0.0;
}

void ObjectShifter::readXML(XMLfileUnits& xmlconfig) {
	ObjectFactory object_factory;
	if(xmlconfig.changecurrentnode("object")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		_obj = std::shared_ptr<Object>(object_factory.create(object_type));
		_obj->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	Coordinate3D shift(xmlconfig, "shift");
	shift.get(_shift);
	global_log->info() << "shift vector: "<< _shift[0] << ", " << _shift[1] << ", " << _shift[2] << endl;
}
