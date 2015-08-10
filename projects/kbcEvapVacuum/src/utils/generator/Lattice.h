/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef LATTICE_H
#define LATTICE_H

#ifdef __cplusplus
extern "C" {
#endif

/** Enum with the 7 Bravais lattice types. */
enum LatticeSystem {
	triclinic = 0,
	monoclinic,
	orthorombic,
	tetragonal,
	rhomboedral,
	hexagonal,
	cubic
};
typedef enum LatticeSystem LatticeSystem;

/** Enum with the four lattice centering types. */
enum LatticeCentering {
	primitive = 0,
	body,
	face,
	base_A,
	base_B,
	base_C
};
typedef enum LatticeCentering LatticeCentering;


/* C interface */

struct lattice_t;
typedef struct lattice_t lattice_t;

lattice_t* lattice_create();
void lattice_destroy(lattice_t* lattice);
void lattice_init(lattice_t* lattice, LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3], long dims[3]);
int lattice_getPoint(lattice_t* lattice, double* r);
int lattice_checkValidity(lattice_t* lattice);
const char* lattice_systemName(lattice_t* lattice);
const char* lattice_centeringName(lattice_t* lattice);


#ifdef __cplusplus
}
#endif


/* C++ interface */
#ifdef __cplusplus

class Lattice {
public:
	/** Lattice constructor */
	Lattice(){}
	/** Lattice destructor */
	~Lattice(){}

    /** Initialize lattice
	 * @param[in]  system     Lattice system
	 * @param[in]  centering  Lattice centering
	 * @param[in]  a          1st lattice vector
	 * @param[in]  b          2nd lattice vector
	 * @param[in]  c          3rd lattice vector
	 * @param[in]  dims       lattice dimension given in multiples of lattice vectors a, b and c
	 */

	void init(LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3], long dims[3]);
	/** Get a point in the grid.
	 * Returns successively all points in the specified lattice
	 * @param[out]  r  pointer where to store the 3 coordinates of the point's coordinates
	 * @return      0 if no point returned, 1 otherwise
	 */
	int getPoint(double* r);

	/** Check if lattice specifications represent a valid Bravais lattice
	 * @return 1 if valid Bravais lattice, 0 otherwise
	 */
	int checkValidity();

    /** Get the name of the lattice system. */
    const char* systemName();
    /** Get the name of the lattice centering. */
    const char* centeringName();


private:
	LatticeSystem _system;
	LatticeCentering _centering;
    /* Lattice vectors */
    double _a[3];
    double _b[3];
    double _c[3];
	/* Lattice dimensions in multiples of a,b,c */
	long _dims[3];

	/* internal counter for current point output */
	long _pos[3]; /* Current lattice cell */
	int _centeringCounter; /* Next centering in lattice cell */
};


struct lattice_t : Lattice {};

#endif /* C++ interface */

#endif /* LATTICE_H */
