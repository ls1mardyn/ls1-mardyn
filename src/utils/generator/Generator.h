/*
 * Copyright (c) 2013-2014 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include <array>
#include <random>

#include "Basis.h"
#include "Lattice.h"
#include "Objects.h"
#include "ObjectFactory.h"
#include "molecules/Molecule.h"

class Object;

/** Lattice generator */
class Generator {
public:
	Generator() : _lattice(), _basis(), _origin{{0.0, 0.0, 0.0}}, _object(nullptr), _latticeOccupancy(1.0), _dis(0.0, 1.0), _gen(0) {}
	~Generator(){}

	/** @brief Read in XML configuration for Generator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * @note This structure is not fixed yet and may see changes
	 * \code{.xml}
	   <generator>
	     <lattice><!-- ... --></lattice>
	     <basis><!-- ... --></basis>
	     <latticeOrigin>
	         <x>DOUBLE</x>
	         <y>DOUBLE</y>
	         <z>DOUBLE</z>
	     </latticeOrigin>
	     <densit>DOUBLE</density>
	     <latticeOccupancy>DOUBLE</latticeOccupancy>
	     <object type=""><!-- ... --></object>
	   </generator>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	/** Initialize the generator
		* @param[in]  lattice  The underlying point lattice to be used
		* @param[in]  basis    The molecular basis to be put on each lattice point
		* @param[in]  origin   The origin for the lattice
		* @param[in]  object   Volume object to be filled
		*/
	void init(Lattice& lattice, Basis& basis, double origin[3], Object *object);

	/** Initialize the generator with current internal state */
	void init();

	/* Set outer bounding box for the generator */
	void setBoudingBox(double bBoxMin[3], double bBoxMax[3]);

	/** Get a single molecule
	 * By subsequent calls all molecules will be returned, one by one.
	 * @param[out] molecule  Pointer to molecule data structure where to store the molecule data (coordinate and component id)
	 * @return     0 if no more molecules can be returned
	 */
	int getMolecule(Molecule *molecule);

private:
	bool isInsideBox(double r[3]);

	Lattice _lattice;
	Basis _basis;
	std::array<double, 3> _origin;
	Object *_object;
	double _latticeOccupancy;

	std::uniform_real_distribution<> _dis;
	std::mt19937 _gen;

	/* Internal values/counters used during the creation by getMolecule */
	long _baseCount;
	double _lattice_point[3];
};

#endif /* GENERATOR_H */
