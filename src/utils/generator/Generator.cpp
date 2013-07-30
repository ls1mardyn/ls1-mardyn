/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "Generator.h"

#include <cmath>


void Generator::init(Lattice& lattice, Basis& basis, double origin[3]) {
	_lattice = lattice;
	_basis = basis;
	for(int d = 0; d < 3; d++) {
		_origin[d] = origin[d];
	}
	_baseCount = 0;
}

int Generator::getMolecule(molecule_t* molecule) {
	if(_baseCount == 0) {
		if(_lattice.getPoint(_lattice_point) == 0) {
			return 0;
		}
	}

	molecule_t molecule_base;
	molecule_base = _basis.getMolecule(_baseCount);

	for(int d = 0; d < 3; d++) {
		molecule->r[d] = _origin[d] + _lattice_point[d] + molecule_base.r[d];
	}
	molecule->cid = molecule_base.cid;

	_baseCount = (_baseCount + 1) % _basis.numMolecules();

	return 1;
}

generator_t* generator_create() {
    return new generator_t();
}

void generator_destroy(generator_t* generator) {
    delete generator;
}

void generator_init(generator_t* generator, lattice_t* lattice, basis_t* basis, double origin[3]) {
    generator->init(*lattice, *basis, origin);
}

int generator_getMolecule(generator_t* generator, molecule_t* molecule) {
    return generator->getMolecule(molecule);
}
