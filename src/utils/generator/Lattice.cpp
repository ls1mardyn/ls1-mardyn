/*
 * Copyright (c) 2012-2014 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "Lattice.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <string>

#include "utils/Logger.h"



/** List of the names of the 7 Bravais lattices */
static const char* LatticeSystemNames[] = {
	"triclinic",
	"monoclinic",
	"orthorombic",
	"tetragonal",
	"rhomboedral",
	"hexagonal",
	"cubic"
};

/** List with the number of cetneres in the different lattice centering types */
static int LatticeCenteringNums[6] = { 1, 2, 4, 2, 2, 2 };

/** List of the names of the centerings */
static const char* LatticeCenteringNames[] = {
	"primitive",
	"body",
	"face",
	"base A",
	"base B",
	"base C"
};

/** Array holding the relative coordinates of the lattice centers in the a,b,c system */
static const double LatticeCenteringCoords[6][4][3] = {
	{ /* primitive */
		{0.0, 0.0, 0.0}
	},
	{ /* body */
		{0.0, 0.0, 0.0},
		{0.5, 0.5, 0.5}
	},
	{ /* face */
		{0.0, 0.0, 0.0},
		{0.5, 0.5, 0.0},
		{0.5, 0.0, 0.5},
		{0.0, 0.5, 0.5}
	},
	{ /* base A */
		{0.0, 0.0, 0.0},
		{0.0, 0.5, 0.5}
	},
	{ /* base B */
		{0.0, 0.0, 0.0},
		{0.5, 0.0, 0.5}
	},
	{ /* base C */
		{0.0, 0.0, 0.0},
		{0.5, 0.5, 0.0}
	}
};

void Lattice::readXML(XMLfileUnits& xmlconfig) {
	std::string latticeSystem;
	xmlconfig.getNodeValue("@system", latticeSystem);
	Log::global_log->info() << "Lattice system type: " << latticeSystem << std::endl;

	std::string latticeCentering;
	xmlconfig.getNodeValue("@centering", latticeCentering);
	Log::global_log->info() << "Lattice centering: " << latticeCentering << std::endl;

	double a[3];
	double b[3];
	double c[3];
	xmlconfig.getNodeValueReduced("vec[@id='a']/x", a[0]);
	xmlconfig.getNodeValueReduced("vec[@id='a']/y", a[1]);
	xmlconfig.getNodeValueReduced("vec[@id='a']/z", a[2]);
	Log::global_log->info() << "Vec a: " << a[0] << ", " << a[1] << ", " << a[2] << std::endl;
	xmlconfig.getNodeValueReduced("vec[@id='b']/x", b[0]);
	xmlconfig.getNodeValueReduced("vec[@id='b']/y", b[1]);
	xmlconfig.getNodeValueReduced("vec[@id='b']/z", b[2]);
	Log::global_log->info() << "Vec b: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
	xmlconfig.getNodeValueReduced("vec[@id='c']/x", c[0]);
	xmlconfig.getNodeValueReduced("vec[@id='c']/y", c[1]);
	xmlconfig.getNodeValueReduced("vec[@id='c']/z", c[2]);
	Log::global_log->info() << "Vec c: " << c[0] << ", " << c[1] << ", " << c[2] << std::endl;

	init( Lattice::system(latticeSystem), Lattice::centering(latticeCentering), a, b, c);
}

void Lattice::init(LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3]) {
	_system    = system;
	_centering = centering;
	for(int d = 0; d < 3; d++) {
		_a[d] = a[d];
		_b[d] = b[d];
		_c[d] = c[d];
	}
	_centeringCounter = 0;
}

int Lattice::getPoint(double* r) {
    /* after iterating over all centering positions for current cell, move to next cell */
	if(_centeringCounter >= LatticeCenteringNums[_centering]) {
		_centeringCounter = 0;
		_pos[0]++;
		/* for hexagonal lattic we have to skip mid points of hexagons */
		if( (_system == hexagonal) && ( (_pos[0] - _pos[1] + 2) % 3 == 0) ) {
			_pos[0]++;
		}
		if(_pos[0] >= _dimsMax[0]) {
			_pos[0] = _dimsMin[0];
			_pos[1]++;
			/* for hexagonal lattic we have to skip mid points of hexagons */
			if( (_system == hexagonal) && ( (_pos[0] - _pos[1] + 2) % 3 == 0) ) {
				_pos[0]++;
			}
			if(_pos[1] >= _dimsMax[1]) {
				_pos[1] = _dimsMin[1];
				_pos[2]++;
				if(_pos[2] >= _dimsMax[2]) {
					return 0;
				}
			}
		}
	}

	double ia = _pos[0] + LatticeCenteringCoords[_centering][_centeringCounter][0];
	double ib = _pos[1] + LatticeCenteringCoords[_centering][_centeringCounter][1];
	double ic = _pos[2] + LatticeCenteringCoords[_centering][_centeringCounter][2];
	for(int d = 0; d < 3; d++) {
		r[d] = ia * _a[d] + ib * _b[d] + ic * _c[d];
	}
	_centeringCounter++;
	return 1;
}

bool Lattice::checkValidity(){
	/* Checks for validity of system centering combination */
	switch(_system) {
		case triclinic:
			switch(_centering) {
				case primitive:
					break;
				default:
					return false;
			}
			break;
		case monoclinic:
			switch(_centering) {
				case primitive:	case base_A: case base_B: case base_C:
					break;
				default:
					return false;
			}
			break;
		case orthorombic:
			switch(_centering) {
				case primitive: case base_A: case base_B: case base_C: case body: case face:
					break;
				default:
					return false;
			}
			break;
		case tetragonal:
			switch(_centering) {
				case primitive: case body:
					break;
				default:
					return false;
			}
			break;
		case rhomboedral:
			switch(_centering) {
				case primitive:
					break;
				default:
					return false;
			}
			break;
		case hexagonal:
			switch(_centering) {
				case primitive:
					break;
				default:
					return false;
			}
			break;
		case cubic:
			switch(_centering) {
				case primitive: case body: case face:
					break;
				default:
					return false;
			}
			break;
		default:
			return 0;
	}
	/** @todo have to implement checking of lattice vectors, as well. */
	return true;
}

void Lattice::setDimsMin(long dimsMin[3]) {
	for(int d = 0; d < 3; d++) {
		_dimsMin[d] = dimsMin[d];
		_pos[d] = dimsMin[d];
	}
}

void Lattice::setDimsMax(long dimsMax[3]) {
	for(int d = 0; d < 3; d++) {
		_dimsMax[d] = dimsMax[d];
	}
}

const char* Lattice::systemName() {
	return LatticeSystemNames[_system];
}
const char* Lattice::centeringName() {
	return LatticeCenteringNames[_centering];
}

int Lattice::numCenters(LatticeCentering centering) {
	return LatticeCenteringNums[centering];
}

LatticeSystem Lattice::system(std::string name) {
	const char* str = name.c_str();
	if( strcmp( LatticeSystemNames[triclinic], str ) == 0 ) { return triclinic; }
	else if( strcmp( LatticeSystemNames[monoclinic], str ) == 0 ) { return monoclinic; }
	else if( strcmp( LatticeSystemNames[orthorombic], str ) == 0 ) { return orthorombic; }
	else if( strcmp( LatticeSystemNames[tetragonal], str ) == 0 ) { return tetragonal; }
	else if( strcmp( LatticeSystemNames[rhomboedral], str ) == 0 ) { return rhomboedral; }
	else if( strcmp( LatticeSystemNames[hexagonal], str ) == 0 ) { return hexagonal; }
	else if( strcmp( LatticeSystemNames[cubic], str ) == 0 ) { return cubic; }
	else { return unknownSystem; }
}

LatticeCentering Lattice::centering(std::string name) {
	const char* str = name.c_str();
	if( strcmp( LatticeCenteringNames[primitive], str ) == 0 ) { return primitive; }
	else if( strcmp( LatticeCenteringNames[body], str ) == 0 ) { return body; }
	else if( strcmp( LatticeCenteringNames[face], str ) == 0 ) { return face; }
	else if( strcmp( LatticeCenteringNames[base_A], str ) == 0 ) { return base_A; }
	else if( strcmp( LatticeCenteringNames[base_B], str ) == 0 ) { return base_B; }
	else if( strcmp( LatticeCenteringNames[base_C], str ) == 0 ) { return base_C; }
	else { return unknownCentering; }
}
