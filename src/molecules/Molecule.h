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

#ifndef MARDYN_WR
	typedef FullMolecule Molecule;
#else
	typedef Molecule_WR Molecule;
#endif

inline FullMolecule* downcastPointerFull(MoleculeInterface* c) {
	assert(static_cast<FullMolecule*>(c) == dynamic_cast<FullMolecule*>(c));
	return static_cast<FullMolecule*>(c);
}

inline FullMolecule& downcastReferenceFull(MoleculeInterface& c) {
	assert(&static_cast<FullMolecule&>(c) == &dynamic_cast<FullMolecule&>(c));
	return static_cast<FullMolecule&>(c);
}

inline Molecule_WR* downcastPointerWR(MoleculeInterface* c) {
	assert(static_cast<Molecule_WR*>(c) == dynamic_cast<Molecule_WR*>(c));
	return static_cast<Molecule_WR*>(c);
}

inline Molecule_WR& downcastReferenceWR(MoleculeInterface& c) {
	assert(&static_cast<Molecule_WR&>(c) == &dynamic_cast<Molecule_WR&>(c));
	return static_cast<Molecule_WR&>(c);
}

#endif /* SRC_MOLECULES_MOLECULE_H_ */
