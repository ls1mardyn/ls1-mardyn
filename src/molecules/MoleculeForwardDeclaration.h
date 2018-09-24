/*
 * MoleculeForwardDeclaration.h
 *
 *  Created on: Jan 28, 2017
 *      Author: tchipevn
 */

#ifndef SRC_MOLECULES_MOLECULEFORWARDDECLARATION_H_
#define SRC_MOLECULES_MOLECULEFORWARDDECLARATION_H_

#ifndef ENABLE_REDUCED_MEMORY_MODE
#ifdef MARDYN_AUTOPAS
	class AutoPasSimpleMolecule;
	typedef AutoPasSimpleMolecule Molecule;
#else
	class FullMolecule;
	typedef FullMolecule Molecule;
#endif

#else
class MoleculeRMM;
typedef MoleculeRMM Molecule;
#endif

#endif /* SRC_MOLECULES_MOLECULEFORWARDDECLARATION_H_ */
