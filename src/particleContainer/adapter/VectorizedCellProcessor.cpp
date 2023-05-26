/**
 * \file
 * \brief VectorizedCellProcessor.cpp
 */

#include "VectorizedCellProcessor.h"
#include "CellDataSoA.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"
#include "particleContainer/FullParticleCell.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include <algorithm>
#include "vectorization/MaskGatherChooser.h"

using namespace Log;
using namespace std;

VectorizedCellProcessor::VectorizedCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		CellProcessor(cutoffRadius, LJcutoffRadius), _domain(domain),
		// maybe move the following to somewhere else:
		_epsRFInvrc3(2. * (domain.getepsilonRF() - 1.) / ((cutoffRadius * cutoffRadius * cutoffRadius) * (2. * domain.getepsilonRF() + 1.))), 
		_eps_sig(), _shift6(), _upot6lj(0.0), _upotXpoles(0.0), _virial(0.0), _myRF(0.0){

#if VCP_VEC_TYPE==VCP_NOVEC
	global_log->info() << "VectorizedCellProcessor: using no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	global_log->info() << "VectorizedCellProcessor: using SSE3 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	global_log->info() << "VectorizedCellProcessor: using AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	global_log->info() << "VectorizedCellProcessor: using AVX2 intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_KNL) || (VCP_VEC_TYPE==VCP_VEC_KNL_GATHER)
	global_log->info() << "VectorizedCellProcessor: using KNL intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_AVX512F) || (VCP_VEC_TYPE==VCP_VEC_AVX512F_GATHER)
	global_log->info() << "VectorizedCellProcessor: using SKX intrinsics." << std::endl;
#endif

	ComponentList components = *(_simulation.getEnsemble()->getComponents());
	// Get the maximum Component ID.
	size_t maxID = 0;
	const ComponentList::const_iterator end = components.end();
	for (ComponentList::const_iterator c = components.begin(); c != end; ++c)
		maxID = std::max(maxID, static_cast<size_t>(c->ID()));

	// Assign a center list start index for each component.
	std::vector<size_t> compIDs;
	compIDs.resize(maxID + 1, 0);
	size_t centers = 0;
	for (ComponentList::const_iterator c = components.begin(); c != end; ++c) {
		compIDs[c->ID()] = centers;
		centers += c->numLJcenters();
	}

	// One row for each LJ Center, one pair (epsilon*24, sigma^2) for each LJ Center in each row.
	_eps_sig.resize(centers, AlignedArray<vcp_real_calc>(centers * 2));
	_shift6.resize(centers, AlignedArray<vcp_real_calc>(centers));

	// Construct the parameter tables.
	for (size_t comp_i = 0; comp_i < components.size(); ++comp_i) {
		for (size_t comp_j = 0; comp_j < components.size(); ++comp_j) {
			ParaStrm & p = _domain.getComp2Params()(components[comp_i].ID(),
					components[comp_j].ID());
			p.reset_read();
			for (size_t center_i = 0;
					center_i < components[comp_i].numLJcenters(); ++center_i) {
				for (size_t center_j = 0;
						center_j < components[comp_j].numLJcenters();
						++center_j) {
					// Extract epsilon*24.0, sigma^2 and shift*6.0 from paramStreams.
					double eps, sig, shift;
					p >> eps;
					p >> sig;
					p >> shift;
					_eps_sig[compIDs[comp_i] + center_i][2 * (compIDs[comp_j] + center_j)] = static_cast <vcp_real_calc>(eps);
					_eps_sig[compIDs[comp_i] + center_i][2 * (compIDs[comp_j] + center_j) + 1] = static_cast<vcp_real_calc>(sig);
					_shift6[compIDs[comp_i] + center_i][compIDs[comp_j] + center_j] = static_cast<vcp_real_calc>(shift);
				}
			}
		}
	}

	// initialize thread data
	_numThreads = mardyn_get_max_threads();
	global_log->info() << "VectorizedCellProcessor: allocate data for " << _numThreads << " threads." << std::endl;
	_threadData.resize(_numThreads);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		VLJCPThreadData * myown = new VLJCPThreadData();
		const int myid = mardyn_get_thread_num();
		_threadData[myid] = myown;
	} // end pragma omp parallel
}

VectorizedCellProcessor :: ~VectorizedCellProcessor () {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int myid = mardyn_get_thread_num();
		delete _threadData[myid];
	}
}


void VectorizedCellProcessor::initTraversal() {
	#if defined(_OPENMP)
	#pragma omp master
	#endif
	{
		_virial = 0.0;
		_upot6lj = 0.0;
		_upotXpoles = 0.0;
		_myRF = 0.0;
	} // end pragma omp master
}


void VectorizedCellProcessor::endTraversal() {
	vcp_real_accum glob_upot6lj = 0.0;
	vcp_real_accum glob_upotXpoles = 0.0;
	vcp_real_accum glob_virial = 0.0;
	vcp_real_accum glob_myRF = 0.0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:glob_upot6lj, glob_upotXpoles, glob_virial, glob_myRF)
	#endif
	{
		const int tid = mardyn_get_thread_num();

		// reduce vectors and clear local variable
		vcp_real_accum thread_upot = 0.0, thread_upotXpoles = 0.0, thread_virial = 0.0, thread_myRF = 0.0;

		load_hSum_Store_Clear(&thread_upot, _threadData[tid]->_upot6ljV);
		load_hSum_Store_Clear(&thread_upotXpoles, _threadData[tid]->_upotXpolesV);
		load_hSum_Store_Clear(&thread_virial, _threadData[tid]->_virialV);
		load_hSum_Store_Clear(&thread_myRF, _threadData[tid]->_myRFV);

		// add to global sum
		glob_upot6lj += thread_upot;
		glob_upotXpoles += thread_upotXpoles;
		glob_virial += thread_virial;
		glob_myRF += thread_myRF;
	} // end pragma omp parallel reduction

	_upot6lj = glob_upot6lj;
	_upotXpoles = glob_upotXpoles;
	_virial = glob_virial;
	_myRF = glob_myRF;
	_domain.setLocalVirial(_virial + 3.0 * _myRF + _domain.getLocalVirial());
	_domain.setLocalUpot(_upot6lj / 6.0 + _upotXpoles + _myRF + _domain.getLocalUpot());
}

	//const DoubleVec minus_one = DoubleVec::set1(-1.0); //currently not used, would produce warning
	const RealCalcVec zero = RealCalcVec::zero();
	const RealCalcVec one = RealCalcVec::set1(1.0);
	const RealCalcVec two = RealCalcVec::set1(2.0);
	const RealCalcVec three = RealCalcVec::set1(3.0);
	const RealCalcVec four = RealCalcVec::set1(4.0);
	const RealCalcVec five = RealCalcVec::set1(5.0);
	const RealCalcVec six = RealCalcVec::set1(6.0);
	const RealCalcVec ten = RealCalcVec::set1(10.0);
	const RealCalcVec _05 = RealCalcVec::set1(0.5);
	const RealCalcVec _075 = RealCalcVec::set1(0.75);
	const RealCalcVec _1pt5 = RealCalcVec::set1(1.5);
	const RealCalcVec _15 = RealCalcVec::set1(15.0);

	template<bool calculateMacroscopic>
	vcp_inline void VectorizedCellProcessor :: _loopBodyLJ(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& sum_upot6lj, RealAccumVec& sum_virial,
			const MaskCalcVec& forceMask,
			const RealCalcVec& eps_24, const RealCalcVec& sig2,
			const RealCalcVec& shift6)
	{
		const RealCalcVec c_dx = r1_x - r2_x;
		const RealCalcVec c_dy = r1_y - r2_y;
		const RealCalcVec c_dz = r1_z - r2_z;

		const RealCalcVec c_r2 = RealCalcVec::scal_prod(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const RealCalcVec r2_inv = RealCalcVec::fastReciprocal_mask(c_r2, forceMask);

		const RealCalcVec lj2 = sig2 * r2_inv;//1FP (scale)
		const RealCalcVec lj4 = lj2 * lj2;//1FP (scale)
		const RealCalcVec lj6 = lj4 * lj2;//1FP (scale)
		const RealCalcVec lj12 = lj6 * lj6;//1FP (scale)
		const RealCalcVec lj12m6 = lj12 - lj6;//1FP (scale)

		const RealCalcVec eps24r2inv = eps_24 * r2_inv;//1FP (scale)
		const RealCalcVec lj12lj12m6 = lj12 + lj12m6;//1FP (scale)
		const RealCalcVec scale = eps24r2inv * lj12lj12m6;//1FP (scale)

		f_x = c_dx * scale;//1FP (apply scale)
		f_y = c_dy * scale;//1FP (apply scale)
		f_z = c_dz * scale;//1FP (apply scale)
		const RealCalcVec m_dx = m1_r_x - m2_r_x;//1FP (virial) (does not count)
		const RealCalcVec m_dy = m1_r_y - m2_r_y;//1FP (virial) (does not count)
		const RealCalcVec m_dz = m1_r_z - m2_r_z;//1FP (virial) (does not count)

		V_x = RealAccumVec::convertCalcToAccum(m_dx * f_x);//1FP (virial)
		V_y = RealAccumVec::convertCalcToAccum(m_dy * f_y);//1FP (virial)
		V_z = RealAccumVec::convertCalcToAccum(m_dz * f_z);//1FP (virial)

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {

			const RealCalcVec upot_sh = RealCalcVec::fmadd(eps_24, lj12m6, shift6); //2 FP upot				//shift6 is not masked -> we have to mask upot_shifted
			const RealCalcVec upot_masked = RealCalcVec::apply_mask(upot_sh, forceMask); //mask it
			const RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(upot_masked);

			sum_upot6lj = sum_upot6lj + upot_accum;//1FP (sum macro)

			sum_virial = sum_virial + V_x + V_y + V_z;//1 FP (sum macro) + 2 FP (virial)
		}
	}


	template<bool calculateMacroscopic>
	vcp_inline void VectorizedCellProcessor :: _loopBodyCharge(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& qii,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			const RealCalcVec& qjj,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
			const MaskCalcVec& forceMask)
	{
		const RealCalcVec c_dx = r1_x - r2_x;
		const RealCalcVec c_dy = r1_y - r2_y;
		const RealCalcVec c_dz = r1_z - r2_z;//fma not possible since they will be reused...

		const RealCalcVec c_dr2 = RealCalcVec::scal_prod(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const RealCalcVec c_dr2_inv = RealCalcVec::fastReciprocal_mask(c_dr2, forceMask);//masked
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
	    const RealCalcVec c_dr_inv = RealCalcVec::fastReciprocSqrt_mask(c_dr2, forceMask);//masked
#else
	    const RealCalcVec c_dr_inv = RealCalcVec::sqrt(c_dr2_inv);//masked
#endif

		const RealCalcVec q1q2per4pie0 = qii * qjj;
		const RealCalcVec upot = q1q2per4pie0 * c_dr_inv;//masked
		const RealCalcVec fac = upot * c_dr2_inv;//masked

		f_x = c_dx * fac;
		f_y = c_dy * fac;
		f_z = c_dz * fac;
		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		V_x = RealAccumVec::convertCalcToAccum(m_dx * f_x);
		V_y = RealAccumVec::convertCalcToAccum(m_dy * f_y);
		V_z = RealAccumVec::convertCalcToAccum(m_dz * f_z);
		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {
			RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(upot);
			sum_upotXpoles = sum_upotXpoles + upot_accum;
			sum_virial = sum_virial + V_x + V_y + V_z;//DoubleVec::scal_prod(m_dx, m_dy, m_dz, f_x, f_y, f_z);
		}
	}

	template<bool calculateMacroscopic>
	vcp_inline void VectorizedCellProcessor :: _loopBodyChargeDipole(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& q,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			const RealCalcVec& e_x, const RealCalcVec& e_y, const RealCalcVec& e_z,
			const RealCalcVec& p,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& M_x, RealAccumVec& M_y, RealAccumVec& M_z,
			RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
			const MaskCalcVec& forceMask)
	{
		const RealCalcVec dx = r1_x - r2_x;
		const RealCalcVec dy = r1_y - r2_y;
		const RealCalcVec dz = r1_z - r2_z;

		const RealCalcVec dr2 = RealCalcVec::scal_prod(dx, dy, dz, dx, dy, dz);

		const RealCalcVec dr2_inv = RealCalcVec::fastReciprocal_mask(dr2, forceMask);
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
		const RealCalcVec dr_inv = RealCalcVec::fastReciprocSqrt_mask(dr2, forceMask);
#else
	    const RealCalcVec dr_inv = RealCalcVec::sqrt(dr2_inv);
#endif
		const RealCalcVec dr3_inv = dr2_inv * dr_inv;

		const RealCalcVec re = RealCalcVec::scal_prod(dx, dy, dz, e_x, e_y, e_z);

		const RealCalcVec qpper4pie0 = q * p;
		const RealCalcVec qpper4pie0dr3 = qpper4pie0 * dr3_inv;

		const RealCalcVec fac = dr2_inv * three * re;

		f_x = qpper4pie0dr3 * RealCalcVec::fnmadd(dx, fac, e_x);
		f_y = qpper4pie0dr3 * RealCalcVec::fnmadd(dy, fac, e_y);
		f_z = qpper4pie0dr3 * RealCalcVec::fnmadd(dz, fac, e_z);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		V_x = RealAccumVec::convertCalcToAccum(m_dx * f_x);
		V_y = RealAccumVec::convertCalcToAccum(m_dy * f_y);
		V_z = RealAccumVec::convertCalcToAccum(m_dz * f_z);

		// Check if we have to add the macroscopic values up.
		if (calculateMacroscopic)
		{
			const RealCalcVec minusUpot =  qpper4pie0dr3 * re;//already masked
			RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(minusUpot);
			sum_upotXpoles = sum_upotXpoles - upot_accum;
			sum_virial = sum_virial + V_x + V_y + V_z; //DoubleVec::scal_prod(m_dx, m_dy, m_dz, f_x, f_y, f_z);//already masked
		}

		const RealCalcVec e_x_dy_minus_e_y_dx = RealCalcVec::fmsub(e_x, dy, e_y * dx);
		const RealCalcVec e_y_dz_minus_e_z_dy = RealCalcVec::fmsub(e_y, dz, e_z * dy);
		const RealCalcVec e_z_dx_minus_e_x_dz = RealCalcVec::fmsub(e_z, dx, e_x * dz);

		M_x = RealAccumVec::convertCalcToAccum(qpper4pie0dr3 * e_y_dz_minus_e_z_dy);
		M_y = RealAccumVec::convertCalcToAccum(qpper4pie0dr3 * e_z_dx_minus_e_x_dz);
		M_z = RealAccumVec::convertCalcToAccum(qpper4pie0dr3 * e_x_dy_minus_e_y_dx);
	}

	template<bool calculateMacroscopic>
	vcp_inline void VectorizedCellProcessor :: _loopBodyDipole(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& eii_x, const RealCalcVec& eii_y, const RealCalcVec& eii_z,
			const RealCalcVec& pii,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
			const RealCalcVec& pjj,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& M1_x, RealAccumVec& M1_y, RealAccumVec& M1_z,
			RealAccumVec& M2_x, RealAccumVec& M2_y, RealAccumVec& M2_z,
			RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial, RealAccumVec& sum_myRF,
			const MaskCalcVec& forceMask,
			const RealCalcVec& epsRFInvrc3)
	{
		const RealCalcVec dx = r1_x - r2_x;
		const RealCalcVec dy = r1_y - r2_y;
		const RealCalcVec dz = r1_z - r2_z;

		const RealCalcVec dr2 = RealCalcVec::scal_prod(dx, dy, dz, dx, dy, dz);

		const RealCalcVec dr2_inv = RealCalcVec::fastReciprocal_mask(dr2, forceMask);
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
	    const RealCalcVec dr_inv = RealCalcVec::fastReciprocSqrt_mask(dr2, forceMask);
#else
	    const RealCalcVec dr_inv = RealCalcVec::sqrt(dr2_inv);
#endif
		const RealCalcVec dr2three_inv = three * dr2_inv;

		const RealCalcVec p1p2 = RealCalcVec::apply_mask(pii * pjj, forceMask);
		const RealCalcVec p1p2per4pie0 = p1p2;
		const RealCalcVec rffac = p1p2 * epsRFInvrc3;

		const RealCalcVec p1p2per4pie0r3 = p1p2per4pie0 * dr_inv * dr2_inv;
		const RealCalcVec p1p2threeper4pie0r5 = p1p2per4pie0r3 * dr2three_inv;

		const RealCalcVec e1e2 = RealCalcVec::scal_prod(eii_x, eii_y, eii_z, ejj_x, ejj_y, ejj_z);
		const RealCalcVec re1 = RealCalcVec::scal_prod(dx, dy, dz, eii_x, eii_y, eii_z);
		const RealCalcVec re2 = RealCalcVec::scal_prod(dx, dy, dz, ejj_x, ejj_y, ejj_z);

		const RealCalcVec re1threeperr2 = re1 * dr2three_inv;
		const RealCalcVec re2threeperr2 = re2 * dr2three_inv;
		const RealCalcVec re1re2perr2 = dr2_inv * re1 * re2;

		const RealCalcVec e1e2minus5re1re2perr2 = RealCalcVec::fnmadd(five, re1re2perr2, e1e2);//-five*re1+e1e2


		f_x = p1p2threeper4pie0r5 * RealCalcVec::scal_prod(dx, eii_x, ejj_x, e1e2minus5re1re2perr2, re2, re1);
		f_y = p1p2threeper4pie0r5 * RealCalcVec::scal_prod(dy, eii_y, ejj_y, e1e2minus5re1re2perr2, re2, re1);
		f_z = p1p2threeper4pie0r5 * RealCalcVec::scal_prod(dz, eii_z, ejj_z, e1e2minus5re1re2perr2, re2, re1);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		V_x = RealAccumVec::convertCalcToAccum(m_dx * f_x);
		V_y = RealAccumVec::convertCalcToAccum(m_dy * f_y);
		V_z = RealAccumVec::convertCalcToAccum(m_dz * f_z);

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {

			// can we precompute some of this?
			const RealCalcVec upot = p1p2per4pie0r3 * RealCalcVec::fnmadd(three, re1re2perr2, e1e2);//already masked
			const RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(upot);

			sum_upotXpoles = sum_upotXpoles + upot_accum;

			//const DoubleVec virial = DoubleVec::scal_prod(m_dx, m_dy, m_dz, f_x, f_y, f_z);//already masked
			sum_virial = sum_virial + V_x + V_y + V_z;

			RealAccumVec rf_accum = RealAccumVec::convertCalcToAccum(rffac * e1e2);
			sum_myRF = sum_myRF + rf_accum;
		}

		const RealCalcVec e1_x_e2_y_minus_e1_y_e2_x = RealCalcVec::fmsub(eii_x, ejj_y, eii_y * ejj_x);
		const RealCalcVec e1_y_e2_z_minus_e1_z_e2_y = RealCalcVec::fmsub(eii_y, ejj_z, eii_z * ejj_y);
		const RealCalcVec e1_z_e2_x_minus_e1_x_e2_z = RealCalcVec::fmsub(eii_z, ejj_x, eii_x * ejj_z);

		M1_x = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(p1p2per4pie0r3, RealCalcVec::fmsub(re2threeperr2, RealCalcVec::fmsub(eii_y, dz, eii_z * dy), e1_y_e2_z_minus_e1_z_e2_y), rffac * e1_y_e2_z_minus_e1_z_e2_y));
		M1_y = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(p1p2per4pie0r3, RealCalcVec::fmsub(re2threeperr2, RealCalcVec::fmsub(eii_z, dx, eii_x * dz), e1_z_e2_x_minus_e1_x_e2_z), rffac * e1_z_e2_x_minus_e1_x_e2_z));
		M1_z = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(p1p2per4pie0r3, RealCalcVec::fmsub(re2threeperr2, RealCalcVec::fmsub(eii_x, dy, eii_y * dx), e1_x_e2_y_minus_e1_y_e2_x), rffac * e1_x_e2_y_minus_e1_y_e2_x));

		M2_x = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(p1p2per4pie0r3, RealCalcVec::fmadd(re1threeperr2, RealCalcVec::fmsub(ejj_y, dz, ejj_z * dy), e1_y_e2_z_minus_e1_z_e2_y), rffac * e1_y_e2_z_minus_e1_z_e2_y));
		M2_y = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(p1p2per4pie0r3, RealCalcVec::fmadd(re1threeperr2, RealCalcVec::fmsub(ejj_z, dx, ejj_x * dz), e1_z_e2_x_minus_e1_x_e2_z), rffac * e1_z_e2_x_minus_e1_x_e2_z));
		M2_z = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(p1p2per4pie0r3, RealCalcVec::fmadd(re1threeperr2, RealCalcVec::fmsub(ejj_x, dy, ejj_y * dx), e1_x_e2_y_minus_e1_y_e2_x), rffac * e1_x_e2_y_minus_e1_y_e2_x));
	}

	template<bool calculateMacroscopic>
	vcp_inline void VectorizedCellProcessor :: _loopBodyChargeQuadrupole(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& q,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
			const RealCalcVec& m,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& M_x, RealAccumVec& M_y, RealAccumVec& M_z,
			RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
			const MaskCalcVec& forceMask) {

		const RealCalcVec c_dx = r1_x - r2_x;
		const RealCalcVec c_dy = r1_y - r2_y;
		const RealCalcVec c_dz = r1_z - r2_z;//fma not possible since they will be reused...

		const RealCalcVec c_dr2 = RealCalcVec::scal_prod(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const RealCalcVec invdr2 = RealCalcVec::fastReciprocal_mask(c_dr2, forceMask);
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
	    const RealCalcVec invdr = RealCalcVec::fastReciprocSqrt_mask(c_dr2, forceMask);
#else
	    const RealCalcVec invdr = RealCalcVec::sqrt(invdr2);
#endif

		const RealCalcVec qQ05per4pie0 = _05 * q * m;

		RealCalcVec costj = RealCalcVec::scal_prod(ejj_x, ejj_y, ejj_z, c_dx, c_dy, c_dz);
		costj = costj * invdr;

		const RealCalcVec qQinv4dr3 = qQ05per4pie0 * invdr * invdr2;
		const RealCalcVec part1 = three * costj * costj;
		const RealCalcVec upot = qQinv4dr3 * (part1 - one);

		/**********
		 * Force
		 **********/
		const RealCalcVec minus_partialRijInvdr = three * upot * invdr2;
		const RealCalcVec partialTjInvdr = six * costj * qQinv4dr3 * invdr;

		const RealCalcVec fac = RealCalcVec::fmadd(costj * partialTjInvdr, invdr, minus_partialRijInvdr);

		f_x = RealCalcVec::fmsub(fac, c_dx, partialTjInvdr * ejj_x);
		f_y = RealCalcVec::fmsub(fac, c_dy, partialTjInvdr * ejj_y);
		f_z = RealCalcVec::fmsub(fac, c_dz, partialTjInvdr * ejj_z);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		V_x = RealAccumVec::convertCalcToAccum(m_dx * f_x);
		V_y = RealAccumVec::convertCalcToAccum(m_dy * f_y);
		V_z = RealAccumVec::convertCalcToAccum(m_dz * f_z);

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {
			const RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(upot);

			sum_upotXpoles = sum_upotXpoles + upot_accum;

			sum_virial = sum_virial + V_x + V_y + V_z;
		}

		/**********
		 * Torque
		 **********/
		const RealCalcVec minuseXrij_x = RealCalcVec::fmsub(ejj_z, c_dy, ejj_y * c_dz);
		const RealCalcVec minuseXrij_y = RealCalcVec::fmsub(ejj_x, c_dz, ejj_z * c_dx);
		const RealCalcVec minuseXrij_z = RealCalcVec::fmsub(ejj_y, c_dx, ejj_x * c_dy);

		M_x = RealAccumVec::convertCalcToAccum(partialTjInvdr * minuseXrij_x);
		M_y = RealAccumVec::convertCalcToAccum(partialTjInvdr * minuseXrij_y);
		M_z = RealAccumVec::convertCalcToAccum(partialTjInvdr * minuseXrij_z);
	}

	template<bool calculateMacroscopic>
	vcp_inline void VectorizedCellProcessor :: _loopBodyDipoleQuadrupole(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& eii_x, const RealCalcVec& eii_y, const RealCalcVec& eii_z,
			const RealCalcVec& p,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
			const RealCalcVec& m,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& M1_x, RealAccumVec& M1_y, RealAccumVec& M1_z,
			RealAccumVec& M2_x, RealAccumVec& M2_y, RealAccumVec& M2_z,
			RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
			const MaskCalcVec& forceMask) {


		const RealCalcVec c_dx = r1_x - r2_x;
		const RealCalcVec c_dy = r1_y - r2_y;
		const RealCalcVec c_dz = r1_z - r2_z;

		const RealCalcVec c_dr2 = RealCalcVec::scal_prod(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const RealCalcVec invdr2 = RealCalcVec::fastReciprocal_mask(c_dr2, forceMask);
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
	    const RealCalcVec invdr = RealCalcVec::fastReciprocSqrt_mask(c_dr2, forceMask);
#else
	    const RealCalcVec invdr = RealCalcVec::sqrt(invdr2);
#endif

		const RealCalcVec myqfac = _1pt5 * p * m * invdr2 * invdr2;

		RealCalcVec costi = RealCalcVec::scal_prod(eii_x, eii_y, eii_z, c_dx, c_dy, c_dz);
		costi = costi * invdr;

		RealCalcVec costj = RealCalcVec::scal_prod(ejj_x, ejj_y, ejj_z, c_dx, c_dy, c_dz);
		costj = costj * invdr;

		const RealCalcVec cos2tj = costj * costj;

		const RealCalcVec cosgij = RealCalcVec::scal_prod(eii_x, eii_y, eii_z, ejj_x, ejj_y, ejj_z);

		/************
		 * Potential
		 ************/
		// TODO: Check if upot has to be multiplied by -1 according to DISS_STOLL S.178.
		// This affects also the implementation in potforce.h
		const RealCalcVec _5cos2tjminus1 = RealCalcVec::fmsub(five, cos2tj, one);
		const RealCalcVec _2costj = two * costj;

		RealCalcVec part1 = costi * _5cos2tjminus1;
		RealCalcVec part2 = _2costj * cosgij;

		RealCalcVec const upot = myqfac * (part2 - part1);

		const RealCalcVec myqfacXinvdr = myqfac * invdr;
		const RealCalcVec minus_partialRijInvdr = four * upot * invdr2;
		const RealCalcVec minus_partialTiInvdr = myqfacXinvdr * _5cos2tjminus1;

		part1 = RealCalcVec::fmsub(five, costi * costj, cosgij); // *-1!

		const RealCalcVec minus_partialTjInvdr = myqfacXinvdr * two * part1;
		const RealCalcVec partialGij = myqfac * _2costj;

		RealCalcVec part3 = RealCalcVec::fmadd(costi, minus_partialTiInvdr, costj * minus_partialTjInvdr);
		const RealCalcVec fac = RealCalcVec::fnmadd(part3, invdr, minus_partialRijInvdr);//minus_pRI - part3*infdr

		// Force components
		f_x = RealCalcVec::scal_prod(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dx, eii_x, ejj_x);

		f_y = RealCalcVec::scal_prod(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dy, eii_y, ejj_y);

		f_z = RealCalcVec::scal_prod(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dz, eii_z, ejj_z);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		V_x = RealAccumVec::convertCalcToAccum(m_dx * f_x);
		V_y = RealAccumVec::convertCalcToAccum(m_dy * f_y);
		V_z = RealAccumVec::convertCalcToAccum(m_dz * f_z);

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {
			const RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(upot);

			sum_upotXpoles = sum_upotXpoles + upot_accum;

			sum_virial = sum_virial + V_x + V_y + V_z;
		}

		/**********
		 * Torque
		 **********/
		const RealCalcVec eii_x_ejj_y_minus_eii_y_ejj_x = RealCalcVec::fmsub(eii_x, ejj_y, eii_y * ejj_x);
		const RealCalcVec eii_y_ejj_z_minus_eii_z_ejj_y = RealCalcVec::fmsub(eii_y, ejj_z, eii_z * ejj_y);
		const RealCalcVec eii_z_ejj_x_minus_eii_x_ejj_z = RealCalcVec::fmsub(eii_z, ejj_x, eii_x * ejj_z);

		const RealCalcVec partialGij_eiXej_x = partialGij * eii_y_ejj_z_minus_eii_z_ejj_y;
		const RealCalcVec partialGij_eiXej_y = partialGij * eii_z_ejj_x_minus_eii_x_ejj_z;
		const RealCalcVec partialGij_eiXej_z = partialGij * eii_x_ejj_y_minus_eii_y_ejj_x;

		RealCalcVec eXrij_x = RealCalcVec::fmsub(eii_y, c_dz, eii_z * c_dy);
		RealCalcVec eXrij_y = RealCalcVec::fmsub(eii_z, c_dx, eii_x * c_dz);
		RealCalcVec eXrij_z = RealCalcVec::fmsub(eii_x, c_dy, eii_y * c_dx);

		M1_x = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(minus_partialTiInvdr, eXrij_x, partialGij_eiXej_x));
		M1_y = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(minus_partialTiInvdr, eXrij_y, partialGij_eiXej_y));
		M1_z = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(minus_partialTiInvdr, eXrij_z, partialGij_eiXej_z));

		eXrij_x = RealCalcVec::fmsub(ejj_y, c_dz, ejj_z * c_dy);
		eXrij_y = RealCalcVec::fmsub(ejj_z, c_dx, ejj_x * c_dz);
		eXrij_z = RealCalcVec::fmsub(ejj_x, c_dy, ejj_y * c_dx);

		M2_x = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(minus_partialTjInvdr, eXrij_x, partialGij_eiXej_x));
		M2_y = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(minus_partialTjInvdr, eXrij_y, partialGij_eiXej_y));
		M2_z = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(minus_partialTjInvdr, eXrij_z, partialGij_eiXej_z));
	}

	template<bool calculateMacroscopic>
	vcp_inline void VectorizedCellProcessor :: _loopBodyQuadrupole(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& eii_x, const RealCalcVec& eii_y, const RealCalcVec& eii_z,
			const RealCalcVec& mii,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
			const RealCalcVec& mjj,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& Mii_x, RealAccumVec& Mii_y, RealAccumVec& Mii_z,
			RealAccumVec& Mjj_x, RealAccumVec& Mjj_y, RealAccumVec& Mjj_z,
			RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
			const MaskCalcVec& forceMask)
	{
		const RealCalcVec c_dx = r1_x - r2_x;
		const RealCalcVec c_dy = r1_y - r2_y;
		const RealCalcVec c_dz = r1_z - r2_z;

		const RealCalcVec c_dr2 = RealCalcVec::scal_prod(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const RealCalcVec invdr2 = RealCalcVec::fastReciprocal_mask(c_dr2, forceMask);
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
	    const RealCalcVec invdr = RealCalcVec::fastReciprocSqrt_mask(c_dr2, forceMask);
#else
	    const RealCalcVec invdr = RealCalcVec::sqrt(invdr2);
#endif

		RealCalcVec qfac = _075 * invdr;
		qfac = qfac * (mii * mjj);
		qfac = qfac * (invdr2 * invdr2);

		RealCalcVec costi = RealCalcVec::scal_prod(eii_x, eii_y, eii_z, c_dx, c_dy, c_dz);
		costi = costi * invdr;

		RealCalcVec costj = RealCalcVec::scal_prod(ejj_x, ejj_y, ejj_z, c_dx, c_dy, c_dz);
		costj = costj * invdr;

		const RealCalcVec cos2ti = costi * costi;
		const RealCalcVec cos2tj = costj * costj;

		const RealCalcVec cosgij = RealCalcVec::scal_prod(eii_x, eii_y, eii_z, ejj_x, ejj_y, ejj_z);

		RealCalcVec term = five * (costi * costj);
		term = cosgij - term;

		/************
		 * Potential
		 ************/
		RealCalcVec part2 = _15 * cos2ti * cos2tj;
		RealCalcVec part3 = two * term * term;
		RealCalcVec upot = RealCalcVec::fmadd(five, (cos2ti + cos2tj), part2);
		upot = (one + part3) - upot;
		upot = qfac * upot;

		/**********
		 * Force
		 **********/
		const RealCalcVec minus_partialRijInvdr = five * upot * invdr2;

		// partialTiInvdr & partialTjInvdr
		RealCalcVec part1 = qfac * ten * invdr;
		part2 = two * term;

		// partialTiInvdr only
		part3 = three * costi * cos2tj;
		RealCalcVec part4 = costi + RealCalcVec::fmadd(part2, costj, part3);
		const RealCalcVec minus_partialTiInvdr = part1 * part4;

		// partialTjInvdr only
		part3 = three * costj * cos2ti;
		part4 = costj + RealCalcVec::fmadd(part2, costi, part3);
		const RealCalcVec minus_partialTjInvdr = part1 * part4;

		const RealCalcVec partialGij = qfac * four * term;

		// fac
		part1 = minus_partialTiInvdr * costi;
		part2 = minus_partialTjInvdr * costj;
		//part3 = (part1 + part2) * invdr;
		const RealCalcVec fac = RealCalcVec::fnmadd((part1 + part2), invdr, minus_partialRijInvdr); //min - part3 = -part3 + min

		// Force components


		f_x = RealCalcVec::scal_prod(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dx, eii_x, ejj_x);

		f_y = RealCalcVec::scal_prod(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dy, eii_y, ejj_y);

		f_z = RealCalcVec::scal_prod(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dz, eii_z, ejj_z);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		V_x = RealAccumVec::convertCalcToAccum(m_dx * f_x);
		V_y = RealAccumVec::convertCalcToAccum(m_dy * f_y);
		V_z = RealAccumVec::convertCalcToAccum(m_dz * f_z);

		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (calculateMacroscopic) {
			const RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(upot);

			sum_upotXpoles = sum_upotXpoles + upot_accum;


			sum_virial = sum_virial + V_x + V_y + V_z;
		}

		/**********
		 * Torque
		 **********/

		const RealCalcVec eii_x_ejj_y_minus_eii_y_ejj_x = RealCalcVec::fmsub(eii_x, ejj_y, eii_y * ejj_x);
		const RealCalcVec eii_y_ejj_z_minus_eii_z_ejj_y = RealCalcVec::fmsub(eii_y, ejj_z, eii_z * ejj_y);
		const RealCalcVec eii_z_ejj_x_minus_eii_x_ejj_z = RealCalcVec::fmsub(eii_z, ejj_x, eii_x * ejj_z);

		const RealCalcVec partialGij_eiXej_x = partialGij * eii_y_ejj_z_minus_eii_z_ejj_y;
		const RealCalcVec partialGij_eiXej_y = partialGij * eii_z_ejj_x_minus_eii_x_ejj_z;
		const RealCalcVec partialGij_eiXej_z = partialGij * eii_x_ejj_y_minus_eii_y_ejj_x;

		RealCalcVec eXrij_x = RealCalcVec::fmsub(eii_y, c_dz, eii_z * c_dy);
		RealCalcVec eXrij_y = RealCalcVec::fmsub(eii_z, c_dx, eii_x * c_dz);
		RealCalcVec eXrij_z = RealCalcVec::fmsub(eii_x, c_dy, eii_y * c_dx);

		Mii_x = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(minus_partialTiInvdr, eXrij_x, partialGij_eiXej_x));
		Mii_y = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(minus_partialTiInvdr, eXrij_y, partialGij_eiXej_y));
		Mii_z = RealAccumVec::convertCalcToAccum(RealCalcVec::fmsub(minus_partialTiInvdr, eXrij_z, partialGij_eiXej_z));

		eXrij_x = RealCalcVec::fmsub(ejj_y, c_dz, ejj_z * c_dy);
		eXrij_y = RealCalcVec::fmsub(ejj_z, c_dx, ejj_x * c_dz);
		eXrij_z = RealCalcVec::fmsub(ejj_x, c_dy, ejj_y * c_dx);

		Mjj_x = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(minus_partialTjInvdr, eXrij_x, partialGij_eiXej_x));
		Mjj_y = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(minus_partialTjInvdr, eXrij_y, partialGij_eiXej_y));
		Mjj_z = RealAccumVec::convertCalcToAccum(RealCalcVec::fmadd(minus_partialTjInvdr, eXrij_z, partialGij_eiXej_z));
	}

template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
void VectorizedCellProcessor::_calculatePairs(CellDataSoA & soa1, CellDataSoA & soa2) {
	const int tid = mardyn_get_thread_num();
	VLJCPThreadData &my_threadData = *_threadData[tid];

	// initialize dist lookups
	soa2.initDistLookupPointers(my_threadData._centers_dist_lookup,
			my_threadData._ljc_dist_lookup, my_threadData._charges_dist_lookup,
			my_threadData._dipoles_dist_lookup,
			my_threadData._quadrupoles_dist_lookup);

	// Pointer for molecules
	const vcp_real_calc * const soa1_mol_pos_x = soa1._mol_pos.xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = soa1._mol_pos.yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = soa1._mol_pos.zBegin();

	//for better readability:
	typedef ConcSites::SiteType SiteType;
	typedef ConcSites::CoordinateType Coordinate;
	typedef CellDataSoA::QuantityType QuantityType;

	// Pointer for LJ centers
	const vcp_real_calc * const soa1_ljc_r_x = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::LJC, Coordinate::X);
	const vcp_real_calc * const soa1_ljc_r_y = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::LJC, Coordinate::Y);
	const vcp_real_calc * const soa1_ljc_r_z = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::LJC, Coordinate::Z);
		 vcp_real_accum * const soa1_ljc_f_x = soa1.getBeginAccum(QuantityType::FORCE, SiteType::LJC, Coordinate::X);
		 vcp_real_accum * const soa1_ljc_f_y = soa1.getBeginAccum(QuantityType::FORCE, SiteType::LJC, Coordinate::Y);
		 vcp_real_accum * const soa1_ljc_f_z = soa1.getBeginAccum(QuantityType::FORCE, SiteType::LJC, Coordinate::Z);
		 vcp_real_accum * const soa1_ljc_V_x = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::LJC, Coordinate::X);
		 vcp_real_accum * const soa1_ljc_V_y = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::LJC, Coordinate::Y);
		 vcp_real_accum * const soa1_ljc_V_z = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::LJC, Coordinate::Z);
	const int * const soa1_mol_ljc_num = soa1._mol_ljc_num;
	const vcp_ljc_id_t * const soa1_ljc_id = soa1._ljc_id;


	const vcp_real_calc * const soa2_ljc_m_r_x = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::LJC, Coordinate::X);
	const vcp_real_calc * const soa2_ljc_m_r_y = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::LJC, Coordinate::Y);
	const vcp_real_calc * const soa2_ljc_m_r_z = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::LJC, Coordinate::Z);
	const vcp_real_calc * const soa2_ljc_r_x = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::LJC, Coordinate::X);
	const vcp_real_calc * const soa2_ljc_r_y = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::LJC, Coordinate::Y);
	const vcp_real_calc * const soa2_ljc_r_z = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::LJC, Coordinate::Z);
		 vcp_real_accum * const soa2_ljc_f_x = soa2.getBeginAccum(QuantityType::FORCE, SiteType::LJC, Coordinate::X);
		 vcp_real_accum * const soa2_ljc_f_y = soa2.getBeginAccum(QuantityType::FORCE, SiteType::LJC, Coordinate::Y);
		 vcp_real_accum * const soa2_ljc_f_z = soa2.getBeginAccum(QuantityType::FORCE, SiteType::LJC, Coordinate::Z);
		 vcp_real_accum * const soa2_ljc_V_x = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::LJC, Coordinate::X);
		 vcp_real_accum * const soa2_ljc_V_y = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::LJC, Coordinate::Y);
		 vcp_real_accum * const soa2_ljc_V_z = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::LJC, Coordinate::Z);
	const vcp_ljc_id_t * const soa2_ljc_id = soa2._ljc_id;

	vcp_lookupOrMask_single* const soa2_ljc_dist_lookup = my_threadData._ljc_dist_lookup;

	// Pointer for charges
	const vcp_real_calc * const soa1_charges_r_x = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::CHARGE, Coordinate::X);
	const vcp_real_calc * const soa1_charges_r_y = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::CHARGE, Coordinate::Y);
	const vcp_real_calc * const soa1_charges_r_z = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa1_charges_f_x = soa1.getBeginAccum(QuantityType::FORCE, SiteType::CHARGE, Coordinate::X);
		 vcp_real_accum * const soa1_charges_f_y = soa1.getBeginAccum(QuantityType::FORCE, SiteType::CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa1_charges_f_z = soa1.getBeginAccum(QuantityType::FORCE, SiteType::CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa1_charges_V_x = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::CHARGE, Coordinate::X);
		 vcp_real_accum * const soa1_charges_V_y = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa1_charges_V_z = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::CHARGE, Coordinate::Z);
	const vcp_real_calc * const soa1_charges_q = soa1._charges_q;
	const int * const soa1_mol_charges_num = soa1._mol_charges_num;

	const vcp_real_calc * const soa2_charges_m_r_x = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::CHARGE, Coordinate::X);
	const vcp_real_calc * const soa2_charges_m_r_y = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::CHARGE, Coordinate::Y);
	const vcp_real_calc * const soa2_charges_m_r_z = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::CHARGE, Coordinate::Z);
	const vcp_real_calc * const soa2_charges_r_x = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::CHARGE, Coordinate::X);
	const vcp_real_calc * const soa2_charges_r_y = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::CHARGE, Coordinate::Y);
	const vcp_real_calc * const soa2_charges_r_z = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa2_charges_f_x = soa2.getBeginAccum(QuantityType::FORCE, SiteType::CHARGE, Coordinate::X);
		 vcp_real_accum * const soa2_charges_f_y = soa2.getBeginAccum(QuantityType::FORCE, SiteType::CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa2_charges_f_z = soa2.getBeginAccum(QuantityType::FORCE, SiteType::CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa2_charges_V_x = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::CHARGE, Coordinate::X);
		 vcp_real_accum * const soa2_charges_V_y = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa2_charges_V_z = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::CHARGE, Coordinate::Z);
	const vcp_real_calc * const soa2_charges_q = soa2._charges_q;

	vcp_lookupOrMask_single* const soa2_charges_dist_lookup = my_threadData._charges_dist_lookup;

	// Pointer for dipoles
	const vcp_real_calc * const soa1_dipoles_r_x = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::DIPOLE, Coordinate::X);
	const vcp_real_calc * const soa1_dipoles_r_y = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::DIPOLE, Coordinate::Y);
	const vcp_real_calc * const soa1_dipoles_r_z = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::DIPOLE, Coordinate::Z);
		 vcp_real_accum * const soa1_dipoles_f_x = soa1.getBeginAccum(QuantityType::FORCE, SiteType::DIPOLE, Coordinate::X);
		 vcp_real_accum * const soa1_dipoles_f_y = soa1.getBeginAccum(QuantityType::FORCE, SiteType::DIPOLE, Coordinate::Y);
		 vcp_real_accum * const soa1_dipoles_f_z = soa1.getBeginAccum(QuantityType::FORCE, SiteType::DIPOLE, Coordinate::Z);
		 vcp_real_accum * const soa1_dipoles_V_x = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::DIPOLE, Coordinate::X);
		 vcp_real_accum * const soa1_dipoles_V_y = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::DIPOLE, Coordinate::Y);
		 vcp_real_accum * const soa1_dipoles_V_z = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::DIPOLE, Coordinate::Z);
	const vcp_real_calc * const soa1_dipoles_p = soa1._dipoles_p;
	const vcp_real_calc * const soa1_dipoles_e_x = soa1._dipoles_e.xBegin();
	const vcp_real_calc * const soa1_dipoles_e_y = soa1._dipoles_e.yBegin();
	const vcp_real_calc * const soa1_dipoles_e_z = soa1._dipoles_e.zBegin();
	vcp_real_accum * const soa1_dipoles_M_x = soa1._dipoles_M.xBegin();
	vcp_real_accum * const soa1_dipoles_M_y = soa1._dipoles_M.yBegin();
	vcp_real_accum * const soa1_dipoles_M_z = soa1._dipoles_M.zBegin();
	const int * const soa1_mol_dipoles_num = soa1._mol_dipoles_num;

	const vcp_real_calc * const soa2_dipoles_m_r_x = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::DIPOLE, Coordinate::X);
	const vcp_real_calc * const soa2_dipoles_m_r_y = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::DIPOLE, Coordinate::Y);
	const vcp_real_calc * const soa2_dipoles_m_r_z = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::DIPOLE, Coordinate::Z);
	const vcp_real_calc * const soa2_dipoles_r_x = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::DIPOLE, Coordinate::X);
	const vcp_real_calc * const soa2_dipoles_r_y = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::DIPOLE, Coordinate::Y);
	const vcp_real_calc * const soa2_dipoles_r_z = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::DIPOLE, Coordinate::Z);
		 vcp_real_accum * const soa2_dipoles_f_x = soa2.getBeginAccum(QuantityType::FORCE, SiteType::DIPOLE, Coordinate::X);
		 vcp_real_accum * const soa2_dipoles_f_y = soa2.getBeginAccum(QuantityType::FORCE, SiteType::DIPOLE, Coordinate::Y);
		 vcp_real_accum * const soa2_dipoles_f_z = soa2.getBeginAccum(QuantityType::FORCE, SiteType::DIPOLE, Coordinate::Z);
		 vcp_real_accum * const soa2_dipoles_V_x = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::DIPOLE, Coordinate::X);
		 vcp_real_accum * const soa2_dipoles_V_y = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::DIPOLE, Coordinate::Y);
		 vcp_real_accum * const soa2_dipoles_V_z = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::DIPOLE, Coordinate::Z);
	const vcp_real_calc * const soa2_dipoles_p = soa2._dipoles_p;
	const vcp_real_calc * const soa2_dipoles_e_x = soa2._dipoles_e.xBegin();
	const vcp_real_calc * const soa2_dipoles_e_y = soa2._dipoles_e.yBegin();
	const vcp_real_calc * const soa2_dipoles_e_z = soa2._dipoles_e.zBegin();
	vcp_real_accum * const soa2_dipoles_M_x = soa2._dipoles_M.xBegin();
	vcp_real_accum * const soa2_dipoles_M_y = soa2._dipoles_M.yBegin();
	vcp_real_accum * const soa2_dipoles_M_z = soa2._dipoles_M.zBegin();

	vcp_lookupOrMask_single* const soa2_dipoles_dist_lookup = my_threadData._dipoles_dist_lookup;

	// Pointer for quadrupoles
	const vcp_real_calc * const soa1_quadrupoles_r_x = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::QUADRUPOLE, Coordinate::X);
	const vcp_real_calc * const soa1_quadrupoles_r_y = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::QUADRUPOLE, Coordinate::Y);
	const vcp_real_calc * const soa1_quadrupoles_r_z = soa1.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::QUADRUPOLE, Coordinate::Z);
		 vcp_real_accum * const soa1_quadrupoles_f_x = soa1.getBeginAccum(QuantityType::FORCE, SiteType::QUADRUPOLE, Coordinate::X);
		 vcp_real_accum * const soa1_quadrupoles_f_y = soa1.getBeginAccum(QuantityType::FORCE, SiteType::QUADRUPOLE, Coordinate::Y);
		 vcp_real_accum * const soa1_quadrupoles_f_z = soa1.getBeginAccum(QuantityType::FORCE, SiteType::QUADRUPOLE, Coordinate::Z);
		 vcp_real_accum * const soa1_quadrupoles_V_x = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::QUADRUPOLE, Coordinate::X);
		 vcp_real_accum * const soa1_quadrupoles_V_y = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::QUADRUPOLE, Coordinate::Y);
		 vcp_real_accum * const soa1_quadrupoles_V_z = soa1.getBeginAccum(QuantityType::VIRIAL, SiteType::QUADRUPOLE, Coordinate::Z);
	const vcp_real_calc * const soa1_quadrupoles_m = soa1._quadrupoles_m;
	const vcp_real_calc * const soa1_quadrupoles_e_x = soa1._quadrupoles_e.xBegin();
	const vcp_real_calc * const soa1_quadrupoles_e_y = soa1._quadrupoles_e.yBegin();
	const vcp_real_calc * const soa1_quadrupoles_e_z = soa1._quadrupoles_e.zBegin();
	     vcp_real_accum * const soa1_quadrupoles_M_x = soa1._quadrupoles_M.xBegin();
	     vcp_real_accum * const soa1_quadrupoles_M_y = soa1._quadrupoles_M.yBegin();
	     vcp_real_accum * const soa1_quadrupoles_M_z = soa1._quadrupoles_M.zBegin();
	const int * const soa1_mol_quadrupoles_num = soa1._mol_quadrupoles_num;


	const vcp_real_calc * const soa2_quadrupoles_m_r_x = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::QUADRUPOLE, Coordinate::X);
	const vcp_real_calc * const soa2_quadrupoles_m_r_y = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::QUADRUPOLE, Coordinate::Y);
	const vcp_real_calc * const soa2_quadrupoles_m_r_z = soa2.getBeginCalc(QuantityType::MOL_POSITION, SiteType::QUADRUPOLE, Coordinate::Z);
	const vcp_real_calc * const soa2_quadrupoles_r_x = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::QUADRUPOLE, Coordinate::X);
	const vcp_real_calc * const soa2_quadrupoles_r_y = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::QUADRUPOLE, Coordinate::Y);
	const vcp_real_calc * const soa2_quadrupoles_r_z = soa2.getBeginCalc(QuantityType::CENTER_POSITION, SiteType::QUADRUPOLE, Coordinate::Z);
		 vcp_real_accum * const soa2_quadrupoles_f_x = soa2.getBeginAccum(QuantityType::FORCE, SiteType::QUADRUPOLE, Coordinate::X);
		 vcp_real_accum * const soa2_quadrupoles_f_y = soa2.getBeginAccum(QuantityType::FORCE, SiteType::QUADRUPOLE, Coordinate::Y);
		 vcp_real_accum * const soa2_quadrupoles_f_z = soa2.getBeginAccum(QuantityType::FORCE, SiteType::QUADRUPOLE, Coordinate::Z);
		 vcp_real_accum * const soa2_quadrupoles_V_x = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::QUADRUPOLE, Coordinate::X);
		 vcp_real_accum * const soa2_quadrupoles_V_y = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::QUADRUPOLE, Coordinate::Y);
		 vcp_real_accum * const soa2_quadrupoles_V_z = soa2.getBeginAccum(QuantityType::VIRIAL, SiteType::QUADRUPOLE, Coordinate::Z);
	const vcp_real_calc * const soa2_quadrupoles_m = soa2._quadrupoles_m;
	const vcp_real_calc * const soa2_quadrupoles_e_x = soa2._quadrupoles_e.xBegin();
	const vcp_real_calc * const soa2_quadrupoles_e_y = soa2._quadrupoles_e.yBegin();
	const vcp_real_calc * const soa2_quadrupoles_e_z = soa2._quadrupoles_e.zBegin();
	     vcp_real_accum * const soa2_quadrupoles_M_x = soa2._quadrupoles_M.xBegin();
	     vcp_real_accum * const soa2_quadrupoles_M_y = soa2._quadrupoles_M.yBegin();
	     vcp_real_accum * const soa2_quadrupoles_M_z = soa2._quadrupoles_M.zBegin();

	vcp_lookupOrMask_single* const soa2_quadrupoles_dist_lookup = my_threadData._quadrupoles_dist_lookup;




	RealAccumVec sum_upot6lj = RealAccumVec::zero();
	RealAccumVec sum_upotXpoles = RealAccumVec::zero();
	RealAccumVec sum_virial = RealAccumVec::zero();
	RealAccumVec sum_myRF = RealAccumVec::zero();

	const RealCalcVec ljrc2 = RealCalcVec::set1(static_cast<vcp_real_calc>(_LJCutoffRadiusSquare));
	const RealCalcVec cutoffRadiusSquare = RealCalcVec::set1(static_cast<vcp_real_calc>(_cutoffRadiusSquare));
	const RealCalcVec epsRFInvrc3 = RealCalcVec::set1(static_cast<vcp_real_calc>(_epsRFInvrc3));

	/*
	 *  Here different end values for the loops are defined. For loops, which do not vectorize over the last (possibly "uneven") amount of indices, the normal values are computed. These mark the end of the vectorized part.
	 *  The longloop values mark the end of the vectorized part, if vectorization is performed for all elements. For these, various conditions have to be fulfilled, to be sure, that no NaN values are stored and no segfaults arise:
	 *  * arrays have to be long enough
	 *  * arrays have to be filled with something
	 *  * _ljc_id has to be set to existing values for each index
	 *  All of these conditions are fulfilled by setting the non existing values within CellDataSoA.h and AlignedArray.h to zero.
	 */
	const size_t end_ljc_j = vcp_floor_to_vec_size(soa2._ljc_num);
	const size_t end_ljc_j_longloop = vcp_ceil_to_vec_size(soa2._ljc_num);//this is ceil _ljc_num, VCP_VEC_SIZE

	const size_t end_charges_j = vcp_floor_to_vec_size(soa2._charges_num);
	const size_t end_charges_j_longloop = vcp_ceil_to_vec_size(soa2._charges_num);//this is ceil _charges_num, VCP_VEC_SIZE

	const size_t end_dipoles_j = vcp_floor_to_vec_size(soa2._dipoles_num);
	const size_t end_dipoles_j_longloop = vcp_ceil_to_vec_size(soa2._dipoles_num);//this is ceil _dipoles_num, VCP_VEC_SIZE

	const size_t end_quadrupoles_j = vcp_floor_to_vec_size(soa2._quadrupoles_num);
	const size_t end_quadrupoles_j_longloop = vcp_ceil_to_vec_size(soa2._quadrupoles_num);//this is ceil _quadrupoles_num, VCP_VEC_SIZE

	size_t i_ljc_idx = 0;
	size_t i_charge_idx = 0;
	size_t i_charge_dipole_idx = 0;
	size_t i_charge_quadrupole_idx = 0;
	size_t i_dipole_charge_idx = 0;
	size_t i_dipole_idx = 0;
	size_t i_dipole_quadrupole_idx = 0;
	size_t i_quadrupole_charge_idx = 0;
	size_t i_quadrupole_dipole_idx = 0;
	size_t i_quadrupole_idx = 0;

	//if(soa1.getMolNum() < 8){
	//	printf("less than 8\n");
	//}

	// Iterate over each center in the first cell.
	const size_t soa1_mol_num = soa1.getMolNum();
	for (size_t i = 0; i < soa1_mol_num; ++i) {//over the molecules
		const RealCalcVec m1_r_x = RealCalcVec::broadcast(soa1_mol_pos_x + i);
		const RealCalcVec m1_r_y = RealCalcVec::broadcast(soa1_mol_pos_y + i);
		const RealCalcVec m1_r_z = RealCalcVec::broadcast(soa1_mol_pos_z + i);
		// Iterate over centers of second cell
		const countertype32 compute_molecule_ljc = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_ljc_idx, soa2._ljc_num,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				ljrc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);
		const countertype32 compute_molecule_charges = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_charge_idx, soa2._charges_num,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);
		const countertype32 compute_molecule_dipoles = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_dipole_idx, soa2._dipoles_num,
				soa2_dipoles_dist_lookup, soa2_dipoles_m_r_x, soa2_dipoles_m_r_y, soa2_dipoles_m_r_z,
				cutoffRadiusSquare,	end_dipoles_j, m1_r_x, m1_r_y, m1_r_z);
		const countertype32 compute_molecule_quadrupoles = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_quadrupole_idx, soa2._quadrupoles_num,
				soa2_quadrupoles_dist_lookup, soa2_quadrupoles_m_r_x, soa2_quadrupoles_m_r_y, soa2_quadrupoles_m_r_z,
				cutoffRadiusSquare, end_quadrupoles_j, m1_r_x, m1_r_y, m1_r_z);

		size_t end_ljc_loop = MaskGatherChooser::getEndloop(end_ljc_j_longloop, compute_molecule_ljc);
		size_t end_charges_loop = MaskGatherChooser::getEndloop(end_charges_j_longloop, compute_molecule_charges);
		size_t end_dipoles_loop = MaskGatherChooser::getEndloop(end_dipoles_j_longloop, compute_molecule_dipoles);
		size_t end_quadrupoles_loop = MaskGatherChooser::getEndloop(end_quadrupoles_j_longloop, compute_molecule_quadrupoles);

		if (compute_molecule_ljc==0) {
			i_ljc_idx += soa1_mol_ljc_num[i];
		}
		else {
			// LJ force computation
			for (int local_i = 0; local_i < soa1_mol_ljc_num[i]; local_i++) {//over the number of lj-centers in the molecule i
				RealAccumVec sum_fx1 = RealAccumVec::zero();
				RealAccumVec sum_fy1 = RealAccumVec::zero();
				RealAccumVec sum_fz1 = RealAccumVec::zero();

				RealAccumVec sum_Vx1 = RealAccumVec::zero();
				RealAccumVec sum_Vy1 = RealAccumVec::zero();
				RealAccumVec sum_Vz1 = RealAccumVec::zero();

				const RealCalcVec c_r_x1 = RealCalcVec::broadcast(soa1_ljc_r_x + i_ljc_idx);
				const RealCalcVec c_r_y1 = RealCalcVec::broadcast(soa1_ljc_r_y + i_ljc_idx);
				const RealCalcVec c_r_z1 = RealCalcVec::broadcast(soa1_ljc_r_z + i_ljc_idx);

				// Iterate over each pair of centers in the second cell.
				size_t j = ForcePolicy::InitJ2(i_ljc_idx);
				for (; j < end_ljc_loop; j += VCP_VEC_SIZE) {//over (all/some) lj-centers -- amount depends on ForcePolicy::InitJ (two cells: all, within one cell: only do i,j not j,i)
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_ljc_dist_lookup, j);
					// Only go on if at least 1 of the forces has to be calculated.
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const RealCalcVec c_r_x2 = MaskGatherChooser::load(soa2_ljc_r_x, j, lookupORforceMask);
						const RealCalcVec c_r_y2 = MaskGatherChooser::load(soa2_ljc_r_y, j, lookupORforceMask);
						const RealCalcVec c_r_z2 = MaskGatherChooser::load(soa2_ljc_r_z, j, lookupORforceMask);

						const RealCalcVec m_r_x2 = MaskGatherChooser::load(soa2_ljc_m_r_x, j, lookupORforceMask);
						const RealCalcVec m_r_y2 = MaskGatherChooser::load(soa2_ljc_m_r_y, j, lookupORforceMask);
						const RealCalcVec m_r_z2 = MaskGatherChooser::load(soa2_ljc_m_r_z, j, lookupORforceMask);

						const vcp_ljc_id_t id_i = soa1_ljc_id[i_ljc_idx];
						RealCalcVec fx, fy, fz;
						RealAccumVec Vx, Vy, Vz;

						RealCalcVec eps_24;
						RealCalcVec sig2;
						unpackEps24Sig2<MaskGatherChooser>(eps_24, sig2, _eps_sig[id_i], soa2_ljc_id, (vcp_ljc_id_t)j, lookupORforceMask);

						RealCalcVec shift6;
						unpackShift6<MaskGatherChooser>(shift6, _shift6[id_i], soa2_ljc_id, (vcp_ljc_id_t)j, lookupORforceMask);

						_loopBodyLJ<CalculateMacroscopic>(
							m1_r_x, m1_r_y, m1_r_z, c_r_x1, c_r_y1, c_r_z1,
							m_r_x2, m_r_y2, m_r_z2, c_r_x2, c_r_y2, c_r_z2,
							fx, fy, fz,
							Vx, Vy, Vz,
							sum_upot6lj, sum_virial,
							MaskGatherChooser::getForceMask(lookupORforceMask),
							eps_24, sig2,
							shift6);

						RealAccumVec a_fx = RealAccumVec::convertCalcToAccum(fx);
						RealAccumVec a_fy = RealAccumVec::convertCalcToAccum(fy);
						RealAccumVec a_fz = RealAccumVec::convertCalcToAccum(fz);

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_ljc_f_x, j, a_fx, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_ljc_f_y, j, a_fy, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_ljc_f_z, j, a_fz, lookupORforceMask);

						sum_fx1 = sum_fx1 + a_fx;
						sum_fy1 = sum_fy1 + a_fy;
						sum_fz1 = sum_fz1 + a_fz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_ljc_V_x, j, Vx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_ljc_V_y, j, Vy, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_ljc_V_z, j, Vz, lookupORforceMask);

						sum_Vx1 = sum_Vx1 + Vx;
						sum_Vy1 = sum_Vy1 + Vy;
						sum_Vz1 = sum_Vz1 + Vz;
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_ljc);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_ljc);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_ljc_dist_lookup, j, remainderM);

						const RealCalcVec c_r_x2 = MaskGatherChooser::load(soa2_ljc_r_x, j, lookupORforceMask);
						const RealCalcVec c_r_y2 = MaskGatherChooser::load(soa2_ljc_r_y, j, lookupORforceMask);
						const RealCalcVec c_r_z2 = MaskGatherChooser::load(soa2_ljc_r_z, j, lookupORforceMask);

						const RealCalcVec m_r_x2 = MaskGatherChooser::load(soa2_ljc_m_r_x, j, lookupORforceMask);
						const RealCalcVec m_r_y2 = MaskGatherChooser::load(soa2_ljc_m_r_y, j, lookupORforceMask);
						const RealCalcVec m_r_z2 = MaskGatherChooser::load(soa2_ljc_m_r_z, j, lookupORforceMask);

						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						RealCalcVec fx, fy, fz;
						RealAccumVec Vx, Vy, Vz;

						RealCalcVec eps_24;
						RealCalcVec sig2;
						unpackEps24Sig2<MaskGatherChooser>(eps_24, sig2, _eps_sig[id_i], soa2_ljc_id, (vcp_ljc_id_t)j, lookupORforceMask);

						RealCalcVec shift6;
						unpackShift6<MaskGatherChooser>(shift6, _shift6[id_i], soa2_ljc_id, (vcp_ljc_id_t)j, lookupORforceMask);

						_loopBodyLJ<CalculateMacroscopic>(
							m1_r_x, m1_r_y, m1_r_z, c_r_x1, c_r_y1, c_r_z1,
							m_r_x2, m_r_y2, m_r_z2, c_r_x2, c_r_y2, c_r_z2,
							fx, fy, fz,
							Vx, Vy, Vz,
							sum_upot6lj, sum_virial,
							remainderM,//use remainder mask as forcemask
							eps_24, sig2,
							shift6);

						RealAccumVec a_fx = RealAccumVec::convertCalcToAccum(fx);
						RealAccumVec a_fy = RealAccumVec::convertCalcToAccum(fy);
						RealAccumVec a_fz = RealAccumVec::convertCalcToAccum(fz);

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_ljc_f_x, j, a_fx, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_ljc_f_y, j, a_fy, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_ljc_f_z, j, a_fz, lookupORforceMask, remainderM);

						sum_fx1 = sum_fx1 + a_fx;
						sum_fy1 = sum_fy1 + a_fy;
						sum_fz1 = sum_fz1 + a_fz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_ljc_V_x, j, Vx, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_ljc_V_y, j, Vy, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_ljc_V_z, j, Vz, lookupORforceMask, remainderM);


						sum_Vx1 = sum_Vx1 + Vx;
						sum_Vy1 = sum_Vy1 + Vy;
						sum_Vz1 = sum_Vz1 + Vz;

					}
				}
#endif

				hSum_Add_Store(soa1_ljc_f_x + i_ljc_idx, sum_fx1);
				hSum_Add_Store(soa1_ljc_f_y + i_ljc_idx, sum_fy1);
				hSum_Add_Store(soa1_ljc_f_z + i_ljc_idx, sum_fz1);

				hSum_Add_Store(soa1_ljc_V_x + i_ljc_idx, sum_Vx1);
				hSum_Add_Store(soa1_ljc_V_y + i_ljc_idx, sum_Vy1);
				hSum_Add_Store(soa1_ljc_V_z + i_ljc_idx, sum_Vz1);

				i_ljc_idx++;
			}
		}
		// Computation of site interactions with charges

		if (compute_molecule_charges == 0) {
			i_charge_idx += soa1_mol_charges_num[i];
			i_dipole_charge_idx += soa1_mol_dipoles_num[i];
			i_quadrupole_charge_idx += soa1_mol_quadrupoles_num[i];
		}
		else {
			// Computation of charge-charge interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++) {

				const RealCalcVec q1 = RealCalcVec::broadcast(soa1_charges_q + i_charge_idx + local_i);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_charges_r_x + i_charge_idx + local_i);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_charges_r_y + i_charge_idx + local_i);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_charges_r_z + i_charge_idx + local_i);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_charge_idx + local_i);

				for (; j < end_charges_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_charges_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const RealCalcVec q2 = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_x, j, a_f_x, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_y, j, a_f_y, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_z, j, a_f_z, lookupORforceMask);

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask);
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_charges_dist_lookup, j, remainderM);

						const RealCalcVec q2 = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, a_f_x, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, a_f_y, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, a_f_z, lookupORforceMask, remainderM);

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask, remainderM);
					}
				}
#endif

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx + local_i, sum_f1_z);
				// Add old virial and summed calculated virials for center 1
				hSum_Add_Store(soa1_charges_V_x + i_charge_idx + local_i, sum_V1_x);
				hSum_Add_Store(soa1_charges_V_y + i_charge_idx + local_i, sum_V1_y);
				hSum_Add_Store(soa1_charges_V_z + i_charge_idx + local_i, sum_V1_z);

			}

			// Computation of dipole-charge interactions

			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++)
			{
				const RealCalcVec p = RealCalcVec::broadcast(soa1_dipoles_p + i_dipole_charge_idx);
				const RealCalcVec e_x = RealCalcVec::broadcast(soa1_dipoles_e_x + i_dipole_charge_idx);
				const RealCalcVec e_y = RealCalcVec::broadcast(soa1_dipoles_e_y + i_dipole_charge_idx);
				const RealCalcVec e_z = RealCalcVec::broadcast(soa1_dipoles_e_z + i_dipole_charge_idx);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_dipoles_r_x + i_dipole_charge_idx);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_dipoles_r_y + i_dipole_charge_idx);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_dipoles_r_z + i_dipole_charge_idx);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				RealAccumVec sum_M_x = RealAccumVec::zero();
				RealAccumVec sum_M_y = RealAccumVec::zero();
				RealAccumVec sum_M_z = RealAccumVec::zero();

				size_t j = ForcePolicy::InitJ2(i_charge_idx);
				for (; j < end_charges_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_charges_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const RealCalcVec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x - a_f_x;//negative, since dipole-charge, not charge-dipole -> direction inversed
						sum_f1_y = sum_f1_y - a_f_y;
						sum_f1_z = sum_f1_z - a_f_z;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_z, j, a_f_z, lookupORforceMask);//newton 3

						//store virials
						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask);//newton 3

						// Store torque

						sum_M_x = sum_M_x + M_x;
						sum_M_y = sum_M_y + M_y;
						sum_M_z = sum_M_z + M_z;
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_charges_dist_lookup, j, remainderM);
						const RealCalcVec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x - a_f_x;
						sum_f1_y = sum_f1_y - a_f_y;
						sum_f1_z = sum_f1_z - a_f_z;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3

						// Store torque

						sum_M_x = sum_M_x + M_x;
						sum_M_y = sum_M_y + M_y;
						sum_M_z = sum_M_z + M_z;
					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_charge_idx, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_charge_idx, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_charge_idx, sum_f1_z);

				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_dipoles_V_x + i_dipole_charge_idx, sum_V1_x);
				hSum_Add_Store(soa1_dipoles_V_y + i_dipole_charge_idx, sum_V1_y);
				hSum_Add_Store(soa1_dipoles_V_z + i_dipole_charge_idx, sum_V1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_charge_idx, sum_M_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_charge_idx, sum_M_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_charge_idx, sum_M_z);

				i_dipole_charge_idx++;
			}

			// Computation of quadrupole-charge interactions

			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++)
			{
				const RealCalcVec m = RealCalcVec::broadcast(soa1_quadrupoles_m + i_quadrupole_charge_idx);
				const RealCalcVec e_x = RealCalcVec::broadcast(soa1_quadrupoles_e_x + i_quadrupole_charge_idx);
				const RealCalcVec e_y = RealCalcVec::broadcast(soa1_quadrupoles_e_y + i_quadrupole_charge_idx);
				const RealCalcVec e_z = RealCalcVec::broadcast(soa1_quadrupoles_e_z + i_quadrupole_charge_idx);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_quadrupoles_r_x + i_quadrupole_charge_idx);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_quadrupoles_r_y + i_quadrupole_charge_idx);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_quadrupoles_r_z + i_quadrupole_charge_idx);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				RealAccumVec sum_M1_x = RealAccumVec::zero();
				RealAccumVec sum_M1_y = RealAccumVec::zero();
				RealAccumVec sum_M1_z = RealAccumVec::zero();

				size_t j = ForcePolicy::InitJ2(i_charge_idx);
				for (; j < end_charges_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_charges_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const RealCalcVec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x - a_f_x;
						sum_f1_y = sum_f1_y - a_f_y;
						sum_f1_z = sum_f1_z - a_f_z;


						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_z, j, a_f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask);//newton 3


						// Store torque

						sum_M1_x = sum_M1_x + M_x;
						sum_M1_y = sum_M1_y + M_y;
						sum_M1_z = sum_M1_z + M_z;
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_charges_dist_lookup, j, remainderM);
						const RealCalcVec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x - a_f_x;
						sum_f1_y = sum_f1_y - a_f_y;
						sum_f1_z = sum_f1_z - a_f_z;


						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;


						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = sum_M1_x + M_x;
						sum_M1_y = sum_M1_y + M_y;
						sum_M1_z = sum_M1_z + M_z;

					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_charge_idx, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_charge_idx, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_charge_idx, sum_f1_z);

				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_quadrupoles_V_x + i_quadrupole_charge_idx, sum_V1_x);
				hSum_Add_Store(soa1_quadrupoles_V_y + i_quadrupole_charge_idx, sum_V1_y);
				hSum_Add_Store(soa1_quadrupoles_V_z + i_quadrupole_charge_idx, sum_V1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_charge_idx, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_charge_idx, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_charge_idx, sum_M1_z);

				i_quadrupole_charge_idx++;
			}

			i_charge_idx += soa1_mol_charges_num[i];
		}

		// Computation of site interactions with dipoles

		// Continue with next molecule if no force has to be calculated
		if (compute_molecule_dipoles==0) {
			i_dipole_idx += soa1_mol_dipoles_num[i];
			i_charge_dipole_idx += soa1_mol_charges_num[i];
			i_quadrupole_dipole_idx += soa1_mol_quadrupoles_num[i];
		}
		else {
			// Computation of dipole-dipole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++) {

				const RealCalcVec p1 = RealCalcVec::broadcast(soa1_dipoles_p + i_dipole_idx + local_i);
				const RealCalcVec e1_x = RealCalcVec::broadcast(soa1_dipoles_e_x + i_dipole_idx + local_i);
				const RealCalcVec e1_y = RealCalcVec::broadcast(soa1_dipoles_e_y + i_dipole_idx + local_i);
				const RealCalcVec e1_z = RealCalcVec::broadcast(soa1_dipoles_e_z + i_dipole_idx + local_i);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_dipoles_r_x + i_dipole_idx + local_i);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_dipoles_r_y + i_dipole_idx + local_i);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_dipoles_r_z + i_dipole_idx + local_i);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				RealAccumVec sum_M1_x = RealAccumVec::zero();
				RealAccumVec sum_M1_y = RealAccumVec::zero();
				RealAccumVec sum_M1_z = RealAccumVec::zero();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_dipole_idx + local_i);
				for (; j < end_dipoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_dipoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const RealCalcVec p2 = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const RealCalcVec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const RealCalcVec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const RealCalcVec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyDipole<CalculateMacroscopic>(
							m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, p1,
							m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p2,
							f_x, f_y, f_z,
							Vx, Vy, Vz,
							M1_x, M1_y, M1_z,
							M2_x, M2_y, M2_z,
							sum_upotXpoles, sum_virial, sum_myRF,
							MaskGatherChooser::getForceMask(lookupORforceMask),
							epsRFInvrc3);

					RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
					RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
					RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_z, j, a_f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask);//newton 3


						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_x, j, M2_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_y, j, M2_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_z, j, M2_z, lookupORforceMask);//newton 3

					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_dipoles_dist_lookup, j, remainderM);
						const RealCalcVec p2 = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const RealCalcVec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const RealCalcVec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const RealCalcVec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyDipole<CalculateMacroscopic>(
							m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, p1,
							m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p2,
							f_x, f_y, f_z,
							Vx, Vy, Vz,
							M1_x, M1_y, M1_z,
							M2_x, M2_y, M2_z,
							sum_upotXpoles, sum_virial, sum_myRF,
							remainderM,
							epsRFInvrc3);

					RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
					RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
					RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_x, j, M2_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_y, j, M2_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_z, j, M2_z, lookupORforceMask, remainderM);//newton 3

					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_idx + local_i, sum_f1_z);

				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_dipoles_V_x + i_dipole_idx + local_i, sum_V1_x);
				hSum_Add_Store(soa1_dipoles_V_y + i_dipole_idx + local_i, sum_V1_y);
				hSum_Add_Store(soa1_dipoles_V_z + i_dipole_idx + local_i, sum_V1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_idx + local_i, sum_M1_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_idx + local_i, sum_M1_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_idx + local_i, sum_M1_z);

			}

			// Computation of charge-dipole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{

				const RealCalcVec q = RealCalcVec::broadcast(soa1_charges_q + i_charge_dipole_idx);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_charges_r_x + i_charge_dipole_idx);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_charges_r_y + i_charge_dipole_idx);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_charges_r_z + i_charge_dipole_idx);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				size_t j = ForcePolicy::InitJ2(i_dipole_idx);
				for (; j < end_dipoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_dipoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const RealCalcVec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);

						const RealCalcVec e_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const RealCalcVec e_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const RealCalcVec e_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_z, j, a_f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask);//newton 3

						// Store torque

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_x, j, M_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_y, j, M_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_z, j, M_z, lookupORforceMask);//newton 3

					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_dipoles_dist_lookup, j, remainderM);
						const RealCalcVec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);

						const RealCalcVec e_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const RealCalcVec e_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const RealCalcVec e_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3

						// Store torque

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_x, j, M_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_y, j, M_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_z, j, M_z, lookupORforceMask, remainderM);//newton 3
					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_dipole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_dipole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_dipole_idx, sum_f1_z);

				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_charges_V_x + i_charge_dipole_idx, sum_V1_x);
				hSum_Add_Store(soa1_charges_V_y + i_charge_dipole_idx, sum_V1_y);
				hSum_Add_Store(soa1_charges_V_z + i_charge_dipole_idx, sum_V1_z);

				i_charge_dipole_idx++;
			}

			// Computation of quadrupole-dipole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++) {

				const RealCalcVec m = RealCalcVec::broadcast(soa1_quadrupoles_m + i_quadrupole_dipole_idx);
				const RealCalcVec e1_x = RealCalcVec::broadcast(soa1_quadrupoles_e_x + i_quadrupole_dipole_idx);
				const RealCalcVec e1_y = RealCalcVec::broadcast(soa1_quadrupoles_e_y + i_quadrupole_dipole_idx);
				const RealCalcVec e1_z = RealCalcVec::broadcast(soa1_quadrupoles_e_z + i_quadrupole_dipole_idx);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_quadrupoles_r_x + i_quadrupole_dipole_idx);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_quadrupoles_r_y + i_quadrupole_dipole_idx);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_quadrupoles_r_z + i_quadrupole_dipole_idx);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				RealAccumVec sum_M1_x = RealAccumVec::zero();
				RealAccumVec sum_M1_y = RealAccumVec::zero();
				RealAccumVec sum_M1_z = RealAccumVec::zero();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_dipole_idx);
				for (; j < end_dipoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_dipoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const RealCalcVec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const RealCalcVec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const RealCalcVec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const RealCalcVec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M2_x, M2_y, M2_z, M1_x, M1_y, M1_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x - a_f_x;
						sum_f1_y = sum_f1_y - a_f_y;
						sum_f1_z = sum_f1_z - a_f_z;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_f_z, j, a_f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask);//newton 3

						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;


						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_x, j, M2_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_y, j, M2_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_z, j, M2_z, lookupORforceMask);//newton 3

					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_dipoles_dist_lookup, j, remainderM);
						const RealCalcVec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const RealCalcVec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const RealCalcVec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const RealCalcVec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M2_x, M2_y, M2_z, M1_x, M1_y, M1_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x - a_f_x;
						sum_f1_y = sum_f1_y - a_f_y;
						sum_f1_z = sum_f1_z - a_f_z;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3



						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;


						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_x, j, M2_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_y, j, M2_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_M_z, j, M2_z, lookupORforceMask, remainderM);//newton 3

					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_dipole_idx, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_dipole_idx, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_dipole_idx, sum_f1_z);

				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_quadrupoles_V_x + i_quadrupole_dipole_idx, sum_V1_x);
				hSum_Add_Store(soa1_quadrupoles_V_y + i_quadrupole_dipole_idx, sum_V1_y);
				hSum_Add_Store(soa1_quadrupoles_V_z + i_quadrupole_dipole_idx, sum_V1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_dipole_idx, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_dipole_idx, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_dipole_idx, sum_M1_z);

				i_quadrupole_dipole_idx++;

			}

			i_dipole_idx += soa1_mol_dipoles_num[i];
		}

		// Computation of site interactions with quadrupoles

		if (compute_molecule_quadrupoles==0) {
			i_quadrupole_idx += soa1_mol_quadrupoles_num[i];
			i_charge_quadrupole_idx += soa1_mol_charges_num[i];
			i_dipole_quadrupole_idx += soa1_mol_dipoles_num[i];
		}
		else {
			// Computation of quadrupole-quadrupole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++)
			{
				const RealCalcVec mii = RealCalcVec::broadcast(soa1_quadrupoles_m + i_quadrupole_idx + local_i);
				const RealCalcVec eii_x = RealCalcVec::broadcast(soa1_quadrupoles_e_x + i_quadrupole_idx + local_i);
				const RealCalcVec eii_y = RealCalcVec::broadcast(soa1_quadrupoles_e_y + i_quadrupole_idx + local_i);
				const RealCalcVec eii_z = RealCalcVec::broadcast(soa1_quadrupoles_e_z + i_quadrupole_idx + local_i);
				const RealCalcVec rii_x = RealCalcVec::broadcast(soa1_quadrupoles_r_x + i_quadrupole_idx + local_i);
				const RealCalcVec rii_y = RealCalcVec::broadcast(soa1_quadrupoles_r_y + i_quadrupole_idx + local_i);
				const RealCalcVec rii_z = RealCalcVec::broadcast(soa1_quadrupoles_r_z + i_quadrupole_idx + local_i);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				RealAccumVec sum_M1_x = RealAccumVec::zero();
				RealAccumVec sum_M1_y = RealAccumVec::zero();
				RealAccumVec sum_M1_z = RealAccumVec::zero();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_quadrupole_idx + local_i);
				for (; j < end_quadrupoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_quadrupoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const RealCalcVec mjj = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const RealCalcVec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const RealCalcVec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const RealCalcVec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const RealCalcVec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const RealCalcVec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const RealCalcVec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, mii,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, mjj,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_z, j, a_f_z, lookupORforceMask);//newton 3


						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask);//newton 3

						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M2_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M2_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M2_z, lookupORforceMask);//newton 3
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_quadrupoles_dist_lookup, j, remainderM);
						const RealCalcVec mjj = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const RealCalcVec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const RealCalcVec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const RealCalcVec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const RealCalcVec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const RealCalcVec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const RealCalcVec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, mii,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, mjj,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M2_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M2_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M2_z, lookupORforceMask, remainderM);//newton 3
					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_idx + local_i, sum_f1_z);

				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_quadrupoles_V_x + i_quadrupole_idx + local_i, sum_V1_x);
				hSum_Add_Store(soa1_quadrupoles_V_y + i_quadrupole_idx + local_i, sum_V1_y);
				hSum_Add_Store(soa1_quadrupoles_V_z + i_quadrupole_idx + local_i, sum_V1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_idx + local_i, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_idx + local_i, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_idx + local_i, sum_M1_z);

			}

			// Computation of charge-quadrupole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{
				const RealCalcVec q = RealCalcVec::broadcast(soa1_charges_q + i_charge_quadrupole_idx);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_charges_r_x + i_charge_quadrupole_idx);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_charges_r_y + i_charge_quadrupole_idx);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_charges_r_z + i_charge_quadrupole_idx);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				size_t j = ForcePolicy::InitJ2(i_quadrupole_idx);
				for (; j < end_quadrupoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_quadrupoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const RealCalcVec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const RealCalcVec e_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const RealCalcVec e_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const RealCalcVec e_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);


						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_z, j, a_f_z, lookupORforceMask);//newton 3

						// Store forces

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask);//newton 3


						// Store torque

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M_z, lookupORforceMask);//newton 3

					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_quadrupoles_dist_lookup, j, remainderM);
						const RealCalcVec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const RealCalcVec e_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const RealCalcVec e_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const RealCalcVec e_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);


						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M_x, M_y, M_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M_z, lookupORforceMask, remainderM);//newton 3

					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_quadrupole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_quadrupole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_quadrupole_idx, sum_f1_z);
				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_charges_V_x + i_charge_quadrupole_idx, sum_V1_x);
				hSum_Add_Store(soa1_charges_V_y + i_charge_quadrupole_idx, sum_V1_y);
				hSum_Add_Store(soa1_charges_V_z + i_charge_quadrupole_idx, sum_V1_z);

				i_charge_quadrupole_idx++;
			}

			// Computation of dipole-quadrupole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++)
			{
				const RealCalcVec p = RealCalcVec::broadcast(soa1_dipoles_p + i_dipole_quadrupole_idx);
				const RealCalcVec eii_x = RealCalcVec::broadcast(soa1_dipoles_e_x + i_dipole_quadrupole_idx);
				const RealCalcVec eii_y = RealCalcVec::broadcast(soa1_dipoles_e_y + i_dipole_quadrupole_idx);
				const RealCalcVec eii_z = RealCalcVec::broadcast(soa1_dipoles_e_z + i_dipole_quadrupole_idx);
				const RealCalcVec rii_x = RealCalcVec::broadcast(soa1_dipoles_r_x + i_dipole_quadrupole_idx);
				const RealCalcVec rii_y = RealCalcVec::broadcast(soa1_dipoles_r_y + i_dipole_quadrupole_idx);
				const RealCalcVec rii_z = RealCalcVec::broadcast(soa1_dipoles_r_z + i_dipole_quadrupole_idx);

				RealAccumVec sum_f1_x = RealAccumVec::zero();
				RealAccumVec sum_f1_y = RealAccumVec::zero();
				RealAccumVec sum_f1_z = RealAccumVec::zero();

				RealAccumVec sum_V1_x = RealAccumVec::zero();
				RealAccumVec sum_V1_y = RealAccumVec::zero();
				RealAccumVec sum_V1_z = RealAccumVec::zero();

				RealAccumVec sum_M1_x = RealAccumVec::zero();
				RealAccumVec sum_M1_y = RealAccumVec::zero();
				RealAccumVec sum_M1_z = RealAccumVec::zero();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_quadrupole_idx);
				for (; j < end_quadrupoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_quadrupoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const RealCalcVec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const RealCalcVec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const RealCalcVec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const RealCalcVec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const RealCalcVec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const RealCalcVec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const RealCalcVec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, p,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_x, j, a_f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_y, j, a_f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_z, j, a_f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask);//newton 3


						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M2_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M2_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M2_z, lookupORforceMask);//newton 3

					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
					const __mmask16 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x0000){
			#else /*VCP_DPDP */
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x00){
			#endif
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_quadrupoles_dist_lookup, j, remainderM);
						const RealCalcVec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const RealCalcVec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const RealCalcVec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const RealCalcVec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const RealCalcVec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const RealCalcVec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const RealCalcVec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealAccumVec M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						RealAccumVec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, p,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_f_x = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_f_y = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_f_z = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_f_x;
						sum_f1_y = sum_f1_y + a_f_y;
						sum_f1_z = sum_f1_z + a_f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_x, j, a_f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_y, j, a_f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_z, j, a_f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = sum_M1_x + M1_x;
						sum_M1_y = sum_M1_y + M1_y;
						sum_M1_z = sum_M1_z + M1_z;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M2_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M2_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M2_z, lookupORforceMask, remainderM);//newton 3

					}
				}
#endif
				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_quadrupole_idx, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_quadrupole_idx, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_quadrupole_idx, sum_f1_z);

				// Add old virials and summed calculated virials for center 1
				hSum_Add_Store(soa1_dipoles_V_x + i_dipole_quadrupole_idx, sum_V1_x);
				hSum_Add_Store(soa1_dipoles_V_y + i_dipole_quadrupole_idx, sum_V1_y);
				hSum_Add_Store(soa1_dipoles_V_z + i_dipole_quadrupole_idx, sum_V1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_quadrupole_idx, sum_M1_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_quadrupole_idx, sum_M1_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_quadrupole_idx, sum_M1_z);

				i_dipole_quadrupole_idx++;

			}

			i_quadrupole_idx += soa1_mol_quadrupoles_num[i];
		}
	}

	sum_upot6lj.aligned_load_add_store(&my_threadData._upot6ljV[0]);
	sum_upotXpoles.aligned_load_add_store(&my_threadData._upotXpolesV[0]);
	sum_virial.aligned_load_add_store(&my_threadData._virialV[0]);
	const RealAccumVec negative_sum_myRF = RealAccumVec::zero() - sum_myRF;
	negative_sum_myRF.aligned_load_add_store(&my_threadData._myRFV[0]);

} // void LennardJonesCellHandler::CalculatePairs_(LJSoA & soa1, LJSoA & soa2)

void VectorizedCellProcessor::processCell(ParticleCell & c) {
	FullParticleCell & full_c = downcastCellReferenceFull(c);

	CellDataSoA& soa = full_c.getCellDataSoA();
	if (c.isHaloCell() or soa.getMolNum() < 2) {
		return;
	}
	const bool CalculateMacroscopic = true;
	const bool ApplyCutoff = true;
	_calculatePairs<SingleCellPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa, soa);
}

void VectorizedCellProcessor::processCellPair(ParticleCell & c1, ParticleCell & c2, bool sumAll) {
	mardyn_assert(&c1 != &c2);
	FullParticleCell & full_c1 = downcastCellReferenceFull(c1);
	FullParticleCell & full_c2 = downcastCellReferenceFull(c2);

	CellDataSoA& soa1 = full_c1.getCellDataSoA();
	CellDataSoA& soa2 = full_c2.getCellDataSoA();

	// if one cell is empty, skip
	if (soa1.getMolNum() == 0 or soa2.getMolNum() == 0) {
		return;
	}

	const bool c1Halo = full_c1.isHaloCell();
	const bool c2Halo = full_c2.isHaloCell();

	// this variable determines whether
	// _calcPairs(soa1, soa2) or _calcPairs(soa2, soa1)
	// is more efficient
	const bool calc_soa1_soa2 = (soa1.getMolNum() <= soa2.getMolNum());

	
	if(sumAll) {

		// Macroscopic conditions: Compute always

		const bool ApplyCutoff = true;

		const bool CalculateMacroscopic = true;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}
	} else {
		// if one cell is empty, or both cells are Halo, skip
		if (c1Halo and c2Halo) {
			return;
		}

		// Macroscopic conditions:
		// if none of the cells is halo, then compute
		// if one of them is halo:
		// 		if full_c1-index < full_c2-index, then compute
		// 		else, then don't compute
		// This saves the Molecule::isLessThan checks
		// and works similar to the "Half-Shell" scheme

		const bool ApplyCutoff = true;

		if ((not c1Halo and not c2Halo) or						// no cell is halo or
				(full_c1.getCellIndex() < full_c2.getCellIndex())) 		// one of them is halo, but full_c1.index < full_c2.index
		{
			const bool CalculateMacroscopic = true;

			if (calc_soa1_soa2) {
				_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
			} else {
				_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
			}

		} else {
			mardyn_assert(c1Halo != c2Halo);							// one of them is halo and
			mardyn_assert(not (full_c1.getCellIndex() < full_c2.getCellIndex()));// full_c1.index not < full_c2.index

			const bool CalculateMacroscopic = false;

			if (calc_soa1_soa2) {
				_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
			} else {
				_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
			}
		}
	}
}

