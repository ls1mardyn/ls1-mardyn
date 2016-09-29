/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef BASIS_H
#define BASIS_H

#include "molecule.h"


#ifdef __cplusplus
extern "C" {
#endif

struct basis_t;
typedef struct basis_t basis_t;


/* C interface */

basis_t* basis_create();
void basis_destroy(basis_t* basis);
void basis_addMolecule(basis_t* basis, molecule_t molecule);
int basis_numMolecules(basis_t* basis);
molecule_t basis_getMolecule(basis_t* basis, int i);

#ifdef __cplusplus
}
#endif


/* C++ interface */
#ifdef __cplusplus

#include <vector>

/** Structure holding the basis used within a unit cell */
class Basis {
public:
	Basis(){}
	~Basis(){}

	/** Add molecule to basis
	 * @param[in]  molecule  Molecule to be added to the basis
	 */
	void addMolecule(molecule_t molecule);

	/** Number of molecules of the basis
	 * @return  number of molecules in the basis
	 */
	int numMolecules();

	/** Obtain molecule from basis
	 * @param[in]  i  Position of molecule to be returned
	 * @return  Molecule at position i
	 */
	molecule_t getMolecule(int i);

private:
    std::vector<molecule_t> _molecules;
};

struct basis_t : Basis {};

#endif /* C++ interface */

#endif /* BASIS_H */
