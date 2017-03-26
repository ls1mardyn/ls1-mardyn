/*
 * Molecule.h
 *
 *  Created on: Jan 28, 2017
 *      Author: tchipevn
 */

#ifndef SRC_MOLECULES_MOLECULE_H_
#define SRC_MOLECULES_MOLECULE_H_

#ifndef MARDYN_WR
	#include "FullMolecule.h"
	typedef FullMolecule Molecule;
#else
	#include "Molecule_WR.h"
	typedef Molecule_WR Molecule;
#endif

#endif /* SRC_MOLECULES_MOLECULE_H_ */
