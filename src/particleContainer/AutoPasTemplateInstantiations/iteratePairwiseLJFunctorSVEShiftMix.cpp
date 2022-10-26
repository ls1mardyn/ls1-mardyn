/**
 * @file iteratePairwiseAVXLJFunctorShiftMix.cpp
 * @author F. Gratl
 * @date 25.06.22
 */

#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "autopas/molecularDynamics/LJFunctorSVE.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::iteratePairwise(
		autopas::LJFunctorSVE<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond