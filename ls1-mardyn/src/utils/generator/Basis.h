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

#include <vector>

#include "molecules/Molecule.h"
#include "utils/xmlfileUnits.h"

/** Structure holding the basis used within a unit cell */
class Basis {
public:
	Basis(){}
	~Basis(){}

	/** @brief Read in XML configuration for Basis and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <basis>
	     <site>
	       <componentid>INT</componentid>
	       <coordinate> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </coordinate>
	     </site>
	     ...
	   </basis>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	/** Add molecule to basis
	 * @param[in]  molecule  Molecule to be added to the basis
	 */
	void addMolecule(Molecule molecule);

	/** Number of molecules of the basis
	 * @return  number of molecules in the basis
	 */
	int numMolecules();

	/** Obtain molecule from basis
	 * @param[in]  i  Position of molecule to be returned
	 * @return  Molecule at position i
	 */
	Molecule getMolecule(int i);

private:
    std::vector<Molecule> _molecules;
};

#endif /* BASIS_H */
