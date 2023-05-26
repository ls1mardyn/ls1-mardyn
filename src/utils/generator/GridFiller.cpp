/*
 * Copyright (c) 2013-2017 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "GridFiller.h"

#include "LA.h"
#include "LesSolver.h"
#include "utils/generator/Objects.h"
#include "utils/Coordinate3D.h"
#include "molecules/Molecule.h"

#include <iostream>
#include <cmath>


void GridFiller::init(Lattice& lattice, Basis& basis, double origin[3]) {
    _lattice = lattice;
    _basis = basis;
    for(int d = 0; d < 3; d++) {
        _origin[d] = origin[d];
    }
    init();
}

void GridFiller::init() {
	double bBoxMin[3];
	double bBoxMax[3];
    _object->getBboxMin(bBoxMin);
    _object->getBboxMax(bBoxMax);
    _baseCount = 0;

    /* transpose (a,b,c) */
    double *A[3];
    for(int i = 0; i < 3; i++) {
        A[i] = new double[3];
        A[i][0] = _lattice.a()[i];
        A[i][1] = _lattice.b()[i];
        A[i][2] = _lattice.c()[i];
    }    /* As we use grid coordinates the origin has to be subtracted.
     * This displacement will be added later on again in getMolecule. */
    double boxMin[3], boxMax[3];
    for(int d = 0; d < 3; d++) {
        boxMin[d] = bBoxMin[d] - _origin[d];
        boxMax[d] = bBoxMax[d] - _origin[d];
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

    for(int i = 0; i < 3; i++) {
    	delete[] A[i];
    }

    _lattice.setDimsMin(startdims);
    _lattice.setDimsMax(enddims);
}

void GridFiller::readXML(XMLfileUnits& xmlconfig) {
	using std::endl;
	if(xmlconfig.changecurrentnode("lattice")) {
		_lattice.readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("basis")) {
		_basis.readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	xmlconfig.getNodeValue("latticeOccupancy", _latticeOccupancy);
	Log::global_log->info() << "Lattice occupancy: " << _latticeOccupancy << std::endl;
	double rho = 0.0;
	if(xmlconfig.getNodeValueReduced("density", rho)) {
		Log::global_log->info() << "Density: " << rho << std::endl;
		rho /= _latticeOccupancy;
		Log::global_log->info() << "Initializing cubic lattice with density (fully occupied): " << rho << std::endl;
		Log::global_log->warning() << "Initializing cubic lattice with density overwrites previously set lattice system and vectors." << std::endl;
		long num_molecules_per_crystal_cell = _basis.numMolecules() * _lattice.numCenters();
		Log::global_log->debug() << "Number of molecules per crystall cell: " << num_molecules_per_crystal_cell << std::endl;
		double crystal_cell_volume = num_molecules_per_crystal_cell / rho;
		double l = pow(crystal_cell_volume, 1./3.);
		double a[3], b[3], c[3];
		a[0] = l; a[1] = 0; a[2] = 0;
		b[0] = 0; b[1] = l; b[2] = 0;
		c[0] = 0; c[1] = 0; c[2] = l;
		_lattice.init(cubic, _lattice.centering(), a, b, c);
	}

	Coordinate3D origin(xmlconfig, "latticeOrigin");
	origin.get(_origin.data());
	Log::global_log->info() << "Origin: " << _origin[0] << ", " << _origin[1] << ", " << _origin[2] << std::endl;
}


int GridFiller::getMolecule(Molecule* molecule) {
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

		if(_dis(_gen) > _latticeOccupancy) {
			continue;
		}

		if(_object->isInside(r)) {
			for(int d = 0; d < 3; d++) {
				molecule->setr(d, r[d]);
				molecule->setv(d, molecule_base.v(d));
				molecule->setD(d, molecule_base.D(d));
				molecule->setComponent(molecule_base.component());
			}
			molecule->setq(molecule_base.q());
			break;
		}
	}
	return 1;
}
