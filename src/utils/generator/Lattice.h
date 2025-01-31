/*
 * Copyright (c) 2012-2017 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef LATTICE_H
#define LATTICE_H

#include <string>

#include "utils/xmlfileUnits.h"

/** Enum with the 7 Bravais lattice types. */
enum LatticeSystem {
    unknownSystem = -1,
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
    unknownCentering = -1,
    primitive = 0,
    body,
    face,
    base_A,
    base_B,
    base_C
};
typedef enum LatticeCentering LatticeCentering;


class Lattice {
public:
	/** Lattice constructor */
	Lattice(){}
	/** Lattice destructor */
	~Lattice(){}

	/** @brief Read in XML configuration for Lattice and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <lattice system="..." centering="...">
	     <vec id='a'> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </vec>
	     <vec id='b'> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </vec>
	     <vec id='c'> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </vec>
	   </lattice>
	   \endcode
	 * where system can be one of the values "triclinic", "monoclinic", "orthorombic", "tetragonal",
	 * "rhomboedral", and "hexagonal", "cubic" (see Bravais lattices) and centering can be one of the values "primitive",
	 * "body", "face", "base A", "base B", and "base C".
	 */
	void readXML(XMLfileUnits& xmlconfig);

	/** Initialize lattice
		* @param[in]  system     Lattice system
		* @param[in]  centering  Lattice centering
		* @param[in]  a          1st lattice vector
		* @param[in]  b          2nd lattice vector
		* @param[in]  c          3rd lattice vector
		*/
	void init(LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3]);

	/** Set lower corner of lattice given in multiples of the lattice vectors and reset current position marker.
		* @param[in]  dimsMin  coordinates of lower corner
		*/
	void setDimsMin(long dimsMin[3]);

	/** Set upper corner of lattice given in multiples of the lattice vectors
		* @param[in]  dimsMax  coordinate of upper corner
		*/
	void setDimsMax(long dimsMax[3]);

	/** Get a point in the grid.
		* Returns successively all points in the specified lattice
		* @param[out]  r  pointer where to store the 3 coordinates of the point's coordinates
		* @return      0 if no point returned, 1 otherwise
		*/
	int getPoint(double* r);

	/** Check if lattice specifications represent a valid Bravais lattice
		* @return true if valid Bravais lattice, false otherwise
		*/
	bool checkValidity();

	void setSystem(LatticeSystem system) { _system = system; }
	LatticeSystem system() { return _system; }

	void setCentering(LatticeCentering centering) { _centering = centering; }
	LatticeCentering centering() { return _centering; }
	int numCenters() { return numCenters(centering()); }

	/** Get the name of the lattice system. */
	const char* systemName();
	/** Get the name of the lattice centering. */
	const char* centeringName();
	/** Get the number of centers for given centering type */
	static int numCenters(LatticeCentering centering);
	/** Get the lattice centering to a given name */
	static LatticeCentering centering(std::string name);
	/** Get the lattice system to a given name */
	static LatticeSystem system(std::string name);
	/** Get pointer to lattice vector a */
	inline const double* a() { return _a; }
	/** Get pointer to lattice vector b */
	inline const double* b() { return _b; }
	/** Get pointer to lattice vector c */
	inline const double* c() { return _c; }
	/** Set i-th element of lattice vector a */
	void seta(short i, double a_i) { _a[i] = a_i; }
	/** Set i-th element of lattice vector b */
	void setb(short i, double b_i) { _b[i] = b_i; }
	/** Set i-th element of lattice vector c */
	void setc(short i, double c_i) { _c[i] = c_i; }


private:
	LatticeSystem _system;
	LatticeCentering _centering;
	/* Lattice vectors */
	double _a[3];
	double _b[3];
	double _c[3];
	/* Lattice dimensions in multiples of a,b,c */
	long _dimsMin[3];
	long _dimsMax[3];

	/* internal counter for current point output */
	long _pos[3]; /* Current lattice cell */
	int _centeringCounter; /* Next centering in lattice cell */
};

#endif /* LATTICE_H */
