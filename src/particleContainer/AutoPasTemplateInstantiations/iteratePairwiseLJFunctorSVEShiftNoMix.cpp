/**
 * @file iteratePairwiseAVXLJFunctorShiftNoMix.cpp
 * @author F. Gratl
 * @date 25.06.22
 */

#ifdef __ARM_FEATURE_SVE
#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "autopas/molecularDynamics/LJFunctorSVE.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::iteratePairwise(
		autopas::LJFunctorSVE<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
#endif