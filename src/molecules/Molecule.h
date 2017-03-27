/*
 * Molecule.h
 *
 *  Created on: Jan 28, 2017
 *      Author: tchipevn
 */

#ifndef SRC_MOLECULES_MOLECULE_H_
#define SRC_MOLECULES_MOLECULE_H_

/**
 * the old class Molecule is now called FullMolecule and it implements MoleculeInterface.
 * Please bear with us and introduce the necessary changes in MoleculeInterface
 * and provide a stub at least for compiling Molecule_WR
 */

#ifndef MARDYN_WR
	#include "FullMolecule.h"
	typedef FullMolecule Molecule;
#else
	#include "Molecule_WR.h"
	typedef Molecule_WR Molecule;
#endif

#endif /* SRC_MOLECULES_MOLECULE_H_ */
