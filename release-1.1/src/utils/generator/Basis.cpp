/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "Basis.h"

void Basis::addMolecule(molecule_t molecule) {
	_molecules.push_back(molecule);
}

molecule_t Basis::getMolecule(int i) {
	return _molecules[i];
}


int Basis::numMolecules(){
	return _molecules.size();
}


basis_t* basis_create() {
    return new basis_t();
}

void basis_addMolecule(basis_t *basis, molecule_t molecule) {
    basis->addMolecule(molecule);
}

int basis_numMolecules(basis_t* basis) {
    return basis->numMolecules();
}

void basis_destroy(basis_t *basis) {
    delete basis;
}
