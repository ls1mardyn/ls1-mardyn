/**
 * @file iteratePairwiseAVXLJFunctorNoShiftMix.cpp
 * @author F. Gratl
 * @date 25.06.22
 */

#ifdef __ARM_FEATURE_SVE
#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "molecularDynamicsLibrary/LJFunctorSVE.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::iteratePairwise(
		mdLib::LJFunctorSVE<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
#endif
