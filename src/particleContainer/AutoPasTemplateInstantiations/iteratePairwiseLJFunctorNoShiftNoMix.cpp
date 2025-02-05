/**
 * @file iteratePairwiseLJFunctorNoShiftNoMix.cpp
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
				/*applyShift*/ false,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
