/**
 * @file iteratePairwiseAVXLJFunctorShiftMix.cpp
 * @author F. Gratl
 * @date 25.06.22
 */

#ifdef __AVX__
#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "molecularDynamicsLibrary/LJFunctorAVX.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::computeInteractions(
		mdLib::LJFunctorAVX<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
#endif
