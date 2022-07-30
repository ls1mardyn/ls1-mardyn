/**
 * @file iteratePairwiseLJFunctorShiftNoMix.cpp
 * @author F. Gratl
 * @date 25.06.22
 */

#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "autopas/molecularDynamics/LJFunctor.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::iteratePairwise(
		autopas::LJFunctor<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
