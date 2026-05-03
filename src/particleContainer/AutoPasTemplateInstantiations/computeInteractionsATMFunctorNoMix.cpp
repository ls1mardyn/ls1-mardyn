/**
 * @file computeInteractionsATMFunctorNoMix.cpp
 * @author muehlhaeusser
 * @date 20.03.25
 */

#ifdef ENABLE_THREE_BODY
#include "autopas/AutoPasImpl.h"
#include "molecules/Molecule.h"
#include "molecularDynamicsLibrary/AxilrodTellerFunctor.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<AutoPasSimpleMolecule>::computeInteractions(
		mdLib::AxilrodTellerFunctor<
				Molecule,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true>
		*);
//! @endcond
#endif
