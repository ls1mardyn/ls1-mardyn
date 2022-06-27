/**
 * @file iteratePairwiseAVXLJFunctorNoShiftNoMix.cpp
 * @author F. Gratl
 * @date 25.06.22
 */

#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::iteratePairwise(
		autopas::LJFunctorAVX<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
