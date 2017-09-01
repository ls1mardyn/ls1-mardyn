/*
 * Copyright (c) 2013-2014 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "Generator.h"

#include "LA.h"
#include "LesSolver.h"
#include "utils/generator/Objects.h"
#include "utils/Coordinate3D.h"
#include "molecules/Molecule.h"

#include <iostream>
#include <cmath>


void Generator::init(Lattice& lattice, Basis& basis, double origin[3], Object *object) {
    _lattice = lattice;
    _basis = basis;
    for(int d = 0; d < 3; d++) {
        _origin[d] = origin[d];
    }
    _object = object;
    init();
}

void Generator::setBoudingBox(double bBoxMin[3], double bBoxMax[3]) {
	Object *bBox = new Cuboid(bBoxMin, bBoxMax);
	Object *boundedObject = new ObjectIntersection(bBox, _object);
	_object =  boundedObject;
}

void Generator::init() {
    _object->getBboxMin(_bBoxMin);
    _object->getBboxMax(_bBoxMax);
    _baseCount = 0;

    /* transpose (a,b,c) */
    double *A[3];
    for(int i = 0; i < 3; i++) {
        A[i] = new double[3];
        A[i][0] = _lattice.a()[i];
        A[i][1] = _lattice.b()[i];
        A[i][2] = _lattice.c()[i];
    }

    /* As we use grid coordinates the origin has to be subtracted.
     * This displacement will be added later on again in getMolecule. */
    double boxMin[3], boxMax[3];
    for(int d = 0; d < 3; d++) {
        boxMin[d] = _bBoxMin[d] - _origin[d];
        boxMax[d] = _bBoxMax[d] - _origin[d];
    }
    double x0[3];
    double x1[3];
    LesSolve(A, (double*) boxMin, 3, x0);
    LesSolve(A, (double*) boxMax, 3, x1);

    long startdims[3];
    long enddims[3];
    for(int d = 0; d < 3; d++) {
        startdims[d] = floor(x0[d]);
        enddims[d] = ceil(x1[d]);
    }

    double rtmp[3];
    double x[3];
    for(int j = 0; j < 3; j++) {
        for(int d = 0; d < 3; d++) {
            rtmp[d] = boxMin[d];
        }
        rtmp[j] = boxMax[j];
        LesSolve(A, rtmp, 3, x);
        for(int d = 0; d < 3; d++) {
            long coord = floor(x[d]);
            if(coord < startdims[d]) { startdims[d] = coord; }
            coord = ceil(x[d]);
            if(coord > enddims[d]) { enddims[d] = coord; }
        }
    }

    _lattice.setDimsMin(startdims);
    _lattice.setDimsMax(enddims);
}

void Generator::readXML(XMLfileUnits& xmlconfig) {
	using std::endl;
	if(xmlconfig.changecurrentnode("lattice")) {
		_lattice.readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("basis")) {
		_basis.readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	Coordinate3D origin(xmlconfig, "latticeOrigin");
	origin.get(_origin);
	global_log->info() << "Origin: " << _origin[0] << ", " << _origin[1] << ", " << _origin[2] << endl;

	if(xmlconfig.changecurrentnode("object")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		ObjectFactory object_factory;
		global_log->debug() << "Obj name: " << object_type << endl;
		_object = object_factory.create(object_type);
		if(_object == nullptr) {
			global_log->debug() << "Unknown object type: " << object_type << endl;
		}
		global_log->error() << "Created object of type: " << _object->getName() << endl;
		_object->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
}


int Generator::getMolecule(Molecule* molecule) {
	for(;;) {
		if(_baseCount == 0) {
			if(_lattice.getPoint(_lattice_point) == 0) {
				return 0;
			}
		}

		Molecule molecule_base;
		molecule_base = _basis.getMolecule(_baseCount);
		double r[3];
		for(int d = 0; d < 3; d++) {
			r[d] = _origin[d] + _lattice_point[d] + molecule_base.r(d);
		}

		_baseCount = (_baseCount + 1) % _basis.numMolecules();

		/* _bBoxMin[0] != _bBoxMax[0] means "no bounding box mode */
		if(_object->isInside(r)) {
			for(int d = 0; d < 3; d++) {
				molecule->setr(d, r[d]);
				molecule->setComponent(molecule_base.component());
			}
			break;
		}
	}
	return 1;
}
