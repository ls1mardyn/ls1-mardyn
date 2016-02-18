/*
 * VectorizedChargeP2PCellProcessor.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: tchipev
 */

#include "bhfmm/cellProcessors/VectorizedChargeP2PCellProcessor.h"
//#include "bhfmm/cellProcessors/CellDataSoA.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"
#include "Domain.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "parallel/DomainDecompBase.h"
#include <iostream>

namespace bhfmm {

using namespace Log;

VectorizedChargeP2PCellProcessor::VectorizedChargeP2PCellProcessor(Domain & domain) :
		CellProcessor(0.0, 0.0),
		_domain(domain), _upotXpoles(0.0), _virial(0.0), _P_xx(0.0), _P_yy(0.0), _P_zz(0.0), _centers_dist_lookup(128) {
#ifdef ENABLE_MPI
	_timer.set_sync(false);
#endif
}

VectorizedChargeP2PCellProcessor :: ~VectorizedChargeP2PCellProcessor () {
}


void VectorizedChargeP2PCellProcessor::initTraversal(const size_t numCells) {
//	_upotXpoles = 0.0;
//	_virial = 0.0;
//	_P_xx = 0.0;
//	_P_yy = 0.0;
//	_P_zz = 0.0;

	global_log->debug() << "VectorizedLJCellProcessor::initTraversal() to " << numCells << " cells." << std::endl;

}


void VectorizedChargeP2PCellProcessor::endTraversal() {
	// divide virial and P_** values by half, as of Rajat's correction in Revision 2674 of potforce.h
//	_virial *= 0.5;
//	_P_xx *= 0.5;
//	_P_yy *= 0.5;
//	_P_zz *= 0.5;
//
//	_domain.addLocalVirial(_virial);
//	_domain.addLocalUpot(_upotXpoles);
//	_domain.addLocalP_xx(_P_xx);
//	_domain.addLocalP_yy(_P_yy);
//	_domain.addLocalP_zz(_P_zz);
}


void VectorizedChargeP2PCellProcessor::preprocessCell(ParticleCell & c) {
}


void VectorizedChargeP2PCellProcessor::postprocessCell(ParticleCell & c) {
}

template<class MacroPolicy>
inline
void VectorizedChargeP2PCellProcessor :: _loopBodyNovecCharges (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
{
}


template<class ForcePolicy>
	unsigned long
inline VectorizedChargeP2PCellProcessor::calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers,
		double* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z
#if VCCP_VEC_TYPE==VCCP_VEC_SSE3
		, size_t end_j, const __m128d m1_r_x, const __m128d m1_r_y, const __m128d m1_r_z
#elif VCCP_VEC_TYPE==VCCP_VEC_AVX
		, size_t end_j, const __m256d m1_r_x, const __m256d m1_r_y, const __m256d m1_r_z
#endif
		) {
	return 0;
}

template<class ForcePolicy, class MacroPolicy>
void VectorizedChargeP2PCellProcessor :: _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
} // void VectorizedChargeP2PCellProcessor::_calculatePairs

void VectorizedChargeP2PCellProcessor::processCell(ParticleCell & c) {
//	_timer.start();
//	_calculatePairs<SingleCellPolicy_, AllMacroPolicy_>(*(c.getCellDataSoACharge()), *(c.getCellDataSoACharge()));
//	_timer.stop();
}

void VectorizedChargeP2PCellProcessor::processCellPair(ParticleCell & c1,
		ParticleCell & c2) {
}

void VectorizedChargeP2PCellProcessor::printTimers() {
//	for (int i = 0; i < numprocs; i++) {
//		if (i == myrank) {
//			std::cout << "rank: " << myrank << std::endl;
//			std::cout << "\t\t" << _timer.get_etime() << "\t\t" << "s in Electro P2P (VectorizedCellProcessor)" << std::endl;
//		}
//		domainDecomp.barrier();
//	}
}

} /* namespace bhfmm */
