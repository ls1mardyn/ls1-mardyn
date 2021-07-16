/*
 * Copyright (c) 2014-2017 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef OBJECTS_H
#define OBJECTS_H

#include <memory>
#include <cstdint>

#include "utils/xmlfileUnits.h"

enum ObjectAxis : uint8_t {
	OBJ_AXIS_X = 0,
	OBJ_AXIS_Y = 1,
	OBJ_AXIS_Z = 2,
};

class Object {
public:
	Object() {}
	virtual void readXML(XMLfileUnits& xmlconfig) {}
	virtual ~Object() {}
	/** Determines if the given point is inside the object */
	virtual bool isInside(double r[3]) = 0;
	/** Determines if the given point is inside the object excluding it's border  */
	virtual bool isInsideNoBorder(double r[3]) = 0;
	/** Get lower corner of a bounding box around the object */
	virtual void getBboxMin(double rmin[3]) = 0;
	/** Get upper corner of a bounding box around the object */
	virtual void getBboxMax(double rmax[3]) = 0;
	/** Get name of object */
	virtual std::string getName() = 0;

	/* Forward getName to getPluginName required by pluginFactory template */
	virtual std::string getPluginName() final { return getName(); };
};

/** Class implementing a cuboid */
class Cuboid : public Object {
public:
	Cuboid();
	Cuboid(double lower[3], double upper[3]);

	/** @brief Read in XML configuration for Cuboid and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <object>
	     <lower> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </lower>
	     <upper> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </upper>
	   </object>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	std::string getName() { return std::string("Cuboid"); }
	static Object* createInstance() { return new Cuboid(); }

	bool isInside(double r[3]);
	bool isInsideNoBorder(double r[3]);
	void getBboxMin(double rmin[3]);
	void getBboxMax(double rmax[3]);

	double upperCorner(int d) { return _upperCorner[d]; }
	double lowerCorner(int d) { return _lowerCorner[d]; }

private:
	double _lowerCorner[3];
	double _upperCorner[3];
};

/** Class implementing a sphere */
class Sphere : public Object {
public:
	Sphere();
	Sphere(double center[3], double r);

	/** @brief Read in XML configuration for Sphere and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <object>
	     <center> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </center>
	     <radius>DOUBLE</radius>
	   </object>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	std::string getName() { return std::string("Sphere"); }
	static Object* createInstance() { return new Sphere(); }

	bool isInside(double r[3]);
	bool isInsideNoBorder(double r[3]);
	void getBboxMin(double rmin[3]);
	void getBboxMax(double rmax[3]);

private:
	double _center[3];
	double _radius;
	double _radiusSquare;
};

/** Class implementing a cyliner */
class Cylinder : public Object {
public:
	Cylinder();
	/** Constructor
	 * @param[in]  centerBase   Center of the circle of the lower base of the cylinder.
	 * @param[in]  radius       Raius of the cylinder (x-y-axis)
	 * @param[in]  height       Height of the cylinder (z-axis)
	 */
	Cylinder(double centerBase[3], double radius, double height);

	/** @brief Read in XML configuration for Cylinder and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <object>
	     <centerBase> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </centerBase>
	     <radius>DOUBLE</radius> <!-- x-y plane -->
	     <height>DOUBLE</height> <!-- z-direction -->
	   </object>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	std::string getName() { return std::string("Cylinder"); }
	static Object* createInstance() { return new Cylinder(); }

	bool isInside(double r[3]);
	bool isInsideNoBorder(double r[3]);
	void getBboxMin(double rmin[3]);
	void getBboxMax(double rmax[3]);

private:
	double _radius;
	double _height;
	double _centerBase[3];
	double _radiusSquare;
	struct CylinderDirection {
		ObjectAxis axis;
		uint8_t base1, base2, height;
	} _direction;
};

/** Abstract class to unify two objects */
class ObjectUnification : public Object {
public:
	ObjectUnification();
	/** Constructor
	 * @param[in]  obj1  First object.
	 * @param[in]  obj2  Second object.
	 */
	ObjectUnification(std::shared_ptr<Object> obj1, std::shared_ptr<Object> obj2) : _ob1(obj1), _ob2(obj2) {}

	/** @brief Read in XML configuration for ObjectUnification and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <object>
	     <object1 type="..."> <!-- ... --> </object1>
	     <object2 type="..."> <!-- ... --> </object2>
	   </object>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	std::string getName() { return std::string("ObjectUnification"); }
	static Object* createInstance() { return new ObjectUnification(); }

	bool isInside(double r[3]) {
		return _ob1->isInside(r) || _ob2->isInside(r);
	}

	bool isInsideNoBorder(double r[3]) {
		/* either inside one of the objects or common border (intersection) */
		return _ob1->isInsideNoBorder(r) || _ob2->isInsideNoBorder(r) || (_ob1->isInside(r) && _ob2->isInside(r));
	}

	void getBboxMin(double rmin[3]) {
		double rmin1[3], rmin2[3];
		_ob1->getBboxMin(rmin1);
		_ob2->getBboxMin(rmin2);
		for(int d = 0; d < 3; d++) {
			rmin[d] = (rmin1[d] < rmin2[d]) ? rmin1[d] : rmin2[d] ;
		}
	}

	void getBboxMax(double rmax[3]) {
		double rmax1[3], rmax2[3];
		_ob1->getBboxMax(rmax1);
		_ob2->getBboxMax(rmax2);
		for(int d = 0; d < 3; d++) {
			rmax[d] = (rmax1[d] > rmax2[d]) ? rmax1[d] : rmax2[d] ;
		}
	}

private:
	std::shared_ptr<Object> _ob1;
	std::shared_ptr<Object> _ob2;
};

/** Abstract class to subtract one object from another
 *
 * @note Boundaries of the subtrahend object lying within the minuend object are included in the resulting object.
 */
class ObjectSubtractor : public Object {
public:
	ObjectSubtractor();
	/** Constructor
	 * @param[in]  original_ob  The original object.
	 * @param[in]  subtract_ob  The object which shall be subtract from the original object.
	 */
	ObjectSubtractor(std::shared_ptr<Object> original_ob, std::shared_ptr<Object> subtract_ob) : _ob1(original_ob), _ob2(subtract_ob) {}

	/** @brief Read in XML configuration for ObjectSubtractor and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <object>
	     <object1 type="..."> <!-- (minuend) ... --> </object1>
	     <object2 type="..."> <!-- (subtrahend) ... --> </object2>
	   </object>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	std::string getName() { return std::string("ObjectSubtractor"); }
	static Object* createInstance() { return new ObjectSubtractor(); }

	bool isInside(double r[3]) {
		return _ob1->isInside(r) && (!_ob2->isInsideNoBorder(r));
	}

	bool isInsideNoBorder(double r[3]) {
		return _ob1->isInsideNoBorder(r) && (!_ob2->isInside(r));
	}

	void getBboxMin(double rmin[3]) {
		_ob1->getBboxMin(rmin);
	}

	void getBboxMax(double rmax[3]) {
		_ob1->getBboxMax(rmax);
	}

private:
	std::shared_ptr<Object> _ob1;
	std::shared_ptr<Object> _ob2;
};

/** Abstract class for the intersection of two objects */
class ObjectIntersection : public Object {
public:
	ObjectIntersection();
	/** Constructor
	 * @param[in]  obj1  First object.
	 * @param[in]  obj2  Second object.
	 */
	ObjectIntersection(std::shared_ptr<Object> obj1, std::shared_ptr<Object> obj2) : _ob1(obj1), _ob2(obj2) {}

	/** @brief Read in XML configuration for ObjectIntersection and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <object>
	     <object1 type="..."> <!-- (summand) ... --> </object1>
	     <object2 type="..."> <!-- (summand) ... --> </object2>
	   </object>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	std::string getName() { return std::string("ObjectIntersection"); }
	static Object* createInstance() { return new ObjectIntersection(); }

	bool isInside(double r[3]) {
		return _ob1->isInside(r) && _ob2->isInside(r);
	}

	bool isInsideNoBorder(double r[3]) {
		return _ob1->isInsideNoBorder(r) && _ob2->isInsideNoBorder(r);
	}

	void getBboxMin(double rmin[3]) {
		double rmin1[3], rmin2[3];
		_ob1->getBboxMin(rmin1);
		_ob2->getBboxMin(rmin2);
		for(int d = 0; d < 3; d++) {
			rmin[d] = (rmin1[d] < rmin2[d]) ? rmin2[d] : rmin1[d] ;
		}
	}

	void getBboxMax(double rmax[3]) {
		double rmax1[3], rmax2[3];
		_ob1->getBboxMax(rmax1);
		_ob2->getBboxMax(rmax2);
		for(int d = 0; d < 3; d++) {
			rmax[d] = (rmax1[d] > rmax2[d]) ? rmax2[d] : rmax1[d] ;
		}
	}

private:
	std::shared_ptr<Object> _ob1;
	std::shared_ptr<Object> _ob2;
};

/** Abstract class to shift an object */
class ObjectShifter : public Object {
public:
	ObjectShifter();
	ObjectShifter(std::shared_ptr<Object> obj, double shift[3]) : _obj(obj) {
		_shift[0] = shift[0];
		_shift[1] = shift[1];
		_shift[2] = shift[2];
	}

	/** @brief Read in XML configuration for ObjectIntersection and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <object>
	     <object type="..."> <!-- ... --> </object>
	     <shift> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </shift>
	   </object>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	std::string getName() { return std::string("ObjectShifter"); }
	static Object* createInstance() { return new ObjectShifter(); }


	bool isInside(double r[3]) {
		double r_transformed[3] = {
			r[0] - _shift[0],
			r[1] - _shift[1],
			r[2] - _shift[2]
		};
		return _obj->isInside(r_transformed);
	}

	bool isInsideNoBorder(double r[3]) {
		double r_transformed[3] = {
			r[0] - _shift[0],
			r[1] - _shift[1],
			r[2] - _shift[2]
		};
		return _obj->isInsideNoBorder(r_transformed);
	}

	void getBboxMin(double rmin[3]) {
		_obj->getBboxMin(rmin);
		for(int d = 0; d < 3; d++) {
			rmin[d] += _shift[d];
		}
	}

	void getBboxMax(double rmax[3]) {
		_obj->getBboxMax(rmax);
		for(int d = 0; d < 3; d++) {
			rmax[d] += _shift[d];
		}
	}

private:
	std::shared_ptr<Object> _obj;
	double _shift[3];
};

#endif /* OBJECTS_H */
