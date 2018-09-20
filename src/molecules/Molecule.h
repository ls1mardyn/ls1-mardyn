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
 * and provide a stub at least for compiling MoleculeRMM
 */

#include "FullMolecule.h"
#include "MoleculeRMM.h"

#ifndef ENABLE_REDUCED_MEMORY_MODE
#ifdef MARDYN_AUTOPAS
	#include "AutoPasFullMolecule.h"
	typedef AutoPasFullMolecule Molecule;
#else
	typedef FullMolecule Molecule;
#endif
#else
	typedef MoleculeRMM Molecule;
#endif

inline FullMolecule* downcastMoleculePointerFull(MoleculeInterface* c) {
	mardyn_assert(static_cast<FullMolecule*>(c) == dynamic_cast<FullMolecule*>(c));
	return static_cast<FullMolecule*>(c);
}

inline FullMolecule& downcastMoleculeReferenceFull(MoleculeInterface& c) {
	mardyn_assert(&static_cast<FullMolecule&>(c) == &dynamic_cast<FullMolecule&>(c));
	return static_cast<FullMolecule&>(c);
}

inline MoleculeRMM* downcastMoleculePointerRMM(MoleculeInterface* c) {
	mardyn_assert(static_cast<MoleculeRMM*>(c) == dynamic_cast<MoleculeRMM*>(c));
	return static_cast<MoleculeRMM*>(c);
}

inline MoleculeRMM& downcastMoleculeReferenceRMM(MoleculeInterface& c) {
	mardyn_assert(&static_cast<MoleculeRMM&>(c) == &dynamic_cast<MoleculeRMM&>(c));
	return static_cast<MoleculeRMM&>(c);
}

#endif /* SRC_MOLECULES_MOLECULE_H_ */
