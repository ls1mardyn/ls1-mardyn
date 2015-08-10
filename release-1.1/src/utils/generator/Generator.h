/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include "Basis.h"
#include "Lattice.h"

#ifdef __cplusplus
extern "C" {
#endif

/* C interface */

struct generator_t;
typedef struct generator_t generator_t;

generator_t* generator_create();
void generator_destroy(generator_t* generator);
void generator_init(generator_t* generator, lattice_t* lattice, basis_t* basis, double origin[3]);
int generator_getMolecule(generator_t* generator, molecule_t *molecule);

#ifdef __cplusplus
}
#endif

/* C++ interface */
#ifdef __cplusplus

/** Lattice generator */
class Generator {
public:
    Generator(){}
    ~Generator(){}

	/** Initialize the generator
	 * @param[in]  lattice  The underlying point lattice to be used
	 * @param[in]  basis    The molecular basis to be put on each lattice point
	 * @param[in]  origin   The origin for the lattice
	 */
    void init(Lattice& lattice, Basis& basis, double origin[3]);

	/** Get a single molecule
	 * By subsequent calls all molecules will be returned, one by one.
	 * @param[out] molecule  Pointer to molecule data structure where to store the molecule data (coordinate and component id)
	 * @return     0 if no more molecules can be returned
	 */
    int getMolecule(molecule_t *molecule);

private:
    Lattice _lattice;
    Basis _basis;
    double _origin[3];

	/* Internal values/counters used during the creation by getMolecule */
	long _baseCount;
	double _lattice_point[3];
};

struct generator_t : Generator {};

#endif /* C++ interface */

#endif /* GENERATOR_H */
