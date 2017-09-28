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

#include "FullMolecule.h"
#include "Molecule_WR.h"

#ifndef ENABLE_REDUCED_MEMORY_MODE
	typedef FullMolecule Molecule;
#else
	typedef Molecule_WR Molecule;
#endif

inline FullMolecule* downcastMoleculePointerFull(MoleculeInterface* c) {
	mardyn_assert(static_cast<FullMolecule*>(c) == dynamic_cast<FullMolecule*>(c));
	return static_cast<FullMolecule*>(c);
}

inline FullMolecule& downcastMoleculeReferenceFull(MoleculeInterface& c) {
	mardyn_assert(&static_cast<FullMolecule&>(c) == &dynamic_cast<FullMolecule&>(c));
	return static_cast<FullMolecule&>(c);
}

inline Molecule_WR* downcastMoleculePointerWR(MoleculeInterface* c) {
	mardyn_assert(static_cast<Molecule_WR*>(c) == dynamic_cast<Molecule_WR*>(c));
	return static_cast<Molecule_WR*>(c);
}

inline Molecule_WR& downcastMoleculeReferenceWR(MoleculeInterface& c) {
	mardyn_assert(&static_cast<Molecule_WR&>(c) == &dynamic_cast<Molecule_WR&>(c));
	return static_cast<Molecule_WR&>(c);
}

#endif /* SRC_MOLECULES_MOLECULE_H_ */
