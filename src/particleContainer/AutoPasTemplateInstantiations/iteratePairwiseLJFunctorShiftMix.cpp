/**
 * @file iteratePairwiseLJFunctorShiftMix.cpp
 * @author F. Gratl
 * @date 25.06.22
 */

#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "molecularDynamicsLibrary/LJFunctor.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::iteratePairwise(
		mdLib::LJFunctor<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
