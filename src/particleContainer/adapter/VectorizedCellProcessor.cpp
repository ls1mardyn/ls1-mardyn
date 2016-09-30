/**
 * \file
 * \brief VectorizedCellProcessor.cpp
 */

#include "VectorizedCellProcessor.h"
#include "CellDataSoA.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"
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
#elif VCP_VEC_TYPE==VCP_VEC_MIC
	global_log->info() << "VectorizedCellProcessor: using MIC intrinsics." << std::endl;
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
	_eps_sig.resize(centers, DoubleArray(centers * 2));
	_shift6.resize(centers, DoubleArray(centers));

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
					p >> _eps_sig[compIDs[comp_i] + center_i][2 * (compIDs[comp_j] + center_j)];
					p >> _eps_sig[compIDs[comp_i] + center_i][2 * (compIDs[comp_j] + center_j) + 1];
					p >> _shift6[compIDs[comp_i] + center_i][compIDs[comp_j] + center_j];
				}
			}
		}
	}

	// initialize thread data
	_numThreads = omp_get_max_threads();
	global_log->info() << "VectorizedCellProcessor: allocate data for " << _numThreads << " threads." << std::endl;
	_threadData.resize(_numThreads);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		VLJCPThreadData * myown = new VLJCPThreadData();
		const int myid = omp_get_thread_num();
		_threadData[myid] = myown;
	} // end pragma omp parallel
}

VectorizedCellProcessor :: ~VectorizedCellProcessor () {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int myid = omp_get_thread_num();
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
	double glob_upot6lj = 0.0;
	double glob_upotXpoles = 0.0;
	double glob_virial = 0.0;
	double glob_myRF = 0.0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:glob_upot6lj, glob_upotXpoles, glob_virial, glob_myRF)
	#endif
	{
		const int tid = omp_get_thread_num();

		// reduce vectors and clear local variable
		double thread_upot = 0.0, thread_upotXpoles = 0.0, thread_virial = 0.0, thread_myRF = 0.0;

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
	_domain.setLocalVirial(_virial + 3.0 * _myRF);
	_domain.setLocalUpot(_upot6lj / 6.0 + _upotXpoles + _myRF);
}

	//const vcp_double_vec minus_one = vcp_simd_set1(-1.0); //currently not used, would produce warning
	const vcp_double_vec zero = VCP_SIMD_ZEROV;
	const vcp_double_vec one = vcp_simd_set1(1.0);
	const vcp_double_vec two = vcp_simd_set1(2.0);
	const vcp_double_vec three = vcp_simd_set1(3.0);
	const vcp_double_vec four = vcp_simd_set1(4.0);
	const vcp_double_vec five = vcp_simd_set1(5.0);
	const vcp_double_vec six = vcp_simd_set1(6.0);
	const vcp_double_vec ten = vcp_simd_set1(10.0);
	const vcp_double_vec _05 = vcp_simd_set1(0.5);
	const vcp_double_vec _075 = vcp_simd_set1(0.75);
	const vcp_double_vec _1pt5 = vcp_simd_set1(1.5);
	const vcp_double_vec _15 = vcp_simd_set1(15.0);

	template<bool calculateMacroscopic>
	inline
	void VectorizedCellProcessor :: _loopBodyLJ(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& sum_upot6lj, vcp_double_vec& sum_virial,
			const vcp_mask_vec& forceMask,
			const vcp_double_vec& eps_24, const vcp_double_vec& sig2,
			const vcp_double_vec& shift6)
	{
		const vcp_double_vec c_dx = r1_x - r2_x;
		const vcp_double_vec c_dy = r1_y - r2_y;
		const vcp_double_vec c_dz = r1_z - r2_z;

		const vcp_double_vec c_r2 = vcp_simd_scalProd(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);
		const vcp_double_vec r2_inv_unmasked = one / c_r2;
		const vcp_double_vec r2_inv = vcp_simd_applymask(r2_inv_unmasked, forceMask);


		const vcp_double_vec lj2 = sig2 * r2_inv;//1FP (scale)
		const vcp_double_vec lj4 = lj2 * lj2;//1FP (scale)
		const vcp_double_vec lj6 = lj4 * lj2;//1FP (scale)
		const vcp_double_vec lj12 = lj6 * lj6;//1FP (scale)
		const vcp_double_vec lj12m6 = lj12 - lj6;//1FP (scale)

		const vcp_double_vec eps24r2inv = eps_24 * r2_inv;//1FP (scale)
		const vcp_double_vec lj12lj12m6 = lj12 + lj12m6;//1FP (scale)
		const vcp_double_vec scale = eps24r2inv * lj12lj12m6;//1FP (scale)

		f_x = c_dx * scale;//1FP (apply scale)
		f_y = c_dy * scale;//1FP (apply scale)
		f_z = c_dz * scale;//1FP (apply scale)
		const vcp_double_vec m_dx = m1_r_x - m2_r_x;//1FP (virial) (does not count)
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;//1FP (virial) (does not count)
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;//1FP (virial) (does not count)

		V_x = m_dx * f_x;//1FP (virial)
		V_y = m_dy * f_y;//1FP (virial)
		V_z = m_dz * f_z;//1FP (virial)

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {

			const vcp_double_vec upot_sh = vcp_simd_fma(eps_24, lj12m6, shift6); //2 FP upot				//shift6 is not masked -> we have to mask upot_shifted
			const vcp_double_vec upot_masked = vcp_simd_applymask(upot_sh, forceMask); //mask it

			sum_upot6lj = sum_upot6lj + upot_masked;//1FP (sum macro)

			sum_virial = sum_virial +  V_x + V_y + V_z;//1 FP (sum macro) + 2 FP (virial)
		}
	}


	template<bool calculateMacroscopic>
	inline void VectorizedCellProcessor :: _loopBodyCharge(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& qii,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			const vcp_double_vec& qjj,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
			const vcp_mask_vec& forceMask)
	{
		const vcp_double_vec c_dx = r1_x - r2_x;
		const vcp_double_vec c_dy = r1_y - r2_y;
		const vcp_double_vec c_dz = r1_z - r2_z;//fma not possible since they will be reused...

		const vcp_double_vec c_dr2 = vcp_simd_scalProd(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const vcp_double_vec c_dr2_inv_unmasked = one / c_dr2;
		const vcp_double_vec c_dr2_inv = vcp_simd_applymask(c_dr2_inv_unmasked, forceMask);//masked
	    const vcp_double_vec c_dr_inv = vcp_simd_sqrt(c_dr2_inv);//masked

		const vcp_double_vec q1q2per4pie0 = qii * qjj;
		const vcp_double_vec upot = q1q2per4pie0 * c_dr_inv;//masked
		const vcp_double_vec fac = upot * c_dr2_inv;//masked

		f_x = c_dx * fac;
		f_y = c_dy * fac;
		f_z = c_dz * fac;
		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		V_x = m_dx * f_x;
		V_y = m_dy * f_y;
		V_z = m_dz * f_z;
		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {
			sum_upotXpoles = sum_upotXpoles + upot;
			sum_virial = sum_virial + V_x + V_y + V_z;//vcp_simd_scalProd(m_dx, m_dy, m_dz, f_x, f_y, f_z);
		}
	}

	template<bool calculateMacroscopic>
	inline void VectorizedCellProcessor :: _loopBodyChargeDipole(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& q,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			const vcp_double_vec& e_x, const vcp_double_vec& e_y, const vcp_double_vec& e_z,
			const vcp_double_vec& p,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
			vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
			const vcp_mask_vec& forceMask)
	{
		const vcp_double_vec dx = r1_x - r2_x;
		const vcp_double_vec dy = r1_y - r2_y;
		const vcp_double_vec dz = r1_z - r2_z;

		const vcp_double_vec dr2 = vcp_simd_scalProd(dx, dy, dz, dx, dy, dz);

		const vcp_double_vec dr2_inv_unmasked = one / dr2;
		const vcp_double_vec dr2_inv = vcp_simd_applymask(dr2_inv_unmasked, forceMask);
		const vcp_double_vec dr_inv = vcp_simd_sqrt(dr2_inv);
		const vcp_double_vec dr3_inv = dr2_inv * dr_inv;

		const vcp_double_vec re = vcp_simd_scalProd(dx, dy, dz, e_x, e_y, e_z);

		const vcp_double_vec qpper4pie0 = q * p;
		const vcp_double_vec qpper4pie0dr3 = qpper4pie0 * dr3_inv;

		const vcp_double_vec fac = dr2_inv * three * re;

		f_x = qpper4pie0dr3 * vcp_simd_fnma(dx, fac, e_x);
		f_y = qpper4pie0dr3 * vcp_simd_fnma(dy, fac, e_y);
		f_z = qpper4pie0dr3 * vcp_simd_fnma(dz, fac, e_z);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		V_x = m_dx * f_x;
		V_y = m_dy * f_y;
		V_z = m_dz * f_z;

		// Check if we have to add the macroscopic values up.
		if (calculateMacroscopic)
		{
			const vcp_double_vec minusUpot =  qpper4pie0dr3 * re;//already masked
			sum_upotXpoles = sum_upotXpoles - minusUpot;
			sum_virial = sum_virial + V_x + V_y + V_z; //vcp_simd_scalProd(m_dx, m_dy, m_dz, f_x, f_y, f_z);//already masked
		}

		const vcp_double_vec e_x_dy_minus_e_y_dx = vcp_simd_fms(e_x, dy, e_y * dx);
		const vcp_double_vec e_y_dz_minus_e_z_dy = vcp_simd_fms(e_y, dz, e_z * dy);
		const vcp_double_vec e_z_dx_minus_e_x_dz = vcp_simd_fms(e_z, dx, e_x * dz);

		M_x = qpper4pie0dr3 * e_y_dz_minus_e_z_dy;
		M_y = qpper4pie0dr3 * e_z_dx_minus_e_x_dz;
		M_z = qpper4pie0dr3 * e_x_dy_minus_e_y_dx;
	}

	template<bool calculateMacroscopic>
	inline void VectorizedCellProcessor :: _loopBodyDipole(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
			const vcp_double_vec& pii,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
			const vcp_double_vec& pjj,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
			vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
			vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial, vcp_double_vec& sum_myRF,
			const vcp_mask_vec& forceMask,
			const vcp_double_vec& epsRFInvrc3)
	{
		const vcp_double_vec dx = r1_x - r2_x;
		const vcp_double_vec dy = r1_y - r2_y;
		const vcp_double_vec dz = r1_z - r2_z;

		const vcp_double_vec dr2 = vcp_simd_scalProd(dx, dy, dz, dx, dy, dz);

		const vcp_double_vec dr2_inv_unmasked = one / dr2;
		const vcp_double_vec dr2_inv = vcp_simd_applymask(dr2_inv_unmasked, forceMask);
		const vcp_double_vec dr_inv = vcp_simd_sqrt(dr2_inv);
		const vcp_double_vec dr2three_inv = three * dr2_inv;

		const vcp_double_vec p1p2 = vcp_simd_applymask(pii * pjj, forceMask);
		const vcp_double_vec p1p2per4pie0 = p1p2;
		const vcp_double_vec rffac = p1p2 * epsRFInvrc3;

		const vcp_double_vec p1p2per4pie0r3 = p1p2per4pie0 * dr_inv * dr2_inv;
		const vcp_double_vec p1p2threeper4pie0r5 = p1p2per4pie0r3 * dr2three_inv;

		const vcp_double_vec e1e2 = vcp_simd_scalProd(eii_x, eii_y, eii_z, ejj_x, ejj_y, ejj_z);
		const vcp_double_vec re1 = vcp_simd_scalProd(dx, dy, dz, eii_x, eii_y, eii_z);
		const vcp_double_vec re2 = vcp_simd_scalProd(dx, dy, dz, ejj_x, ejj_y, ejj_z);

		const vcp_double_vec re1threeperr2 = re1 * dr2three_inv;
		const vcp_double_vec re2threeperr2 = re2 * dr2three_inv;
		const vcp_double_vec re1re2perr2 = dr2_inv * re1 * re2;

		const vcp_double_vec e1e2minus5re1re2perr2 = vcp_simd_fnma(five, re1re2perr2, e1e2);//-five*re1+e1e2


		f_x = p1p2threeper4pie0r5 * vcp_simd_scalProd(dx, eii_x, ejj_x, e1e2minus5re1re2perr2, re2, re1);
		f_y = p1p2threeper4pie0r5 * vcp_simd_scalProd(dy, eii_y, ejj_y, e1e2minus5re1re2perr2, re2, re1);
		f_z = p1p2threeper4pie0r5 * vcp_simd_scalProd(dz, eii_z, ejj_z, e1e2minus5re1re2perr2, re2, re1);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		V_x = m_dx * f_x;
		V_y = m_dy * f_y;
		V_z = m_dz * f_z;

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {

			// can we precompute some of this?
			const vcp_double_vec upot = p1p2per4pie0r3 * vcp_simd_fnma(three, re1re2perr2, e1e2);//already masked
			sum_upotXpoles = sum_upotXpoles + upot;

			//const vcp_double_vec virial = vcp_simd_scalProd(m_dx, m_dy, m_dz, f_x, f_y, f_z);//already masked
			sum_virial = sum_virial + V_x + V_y + V_z;

			sum_myRF = vcp_simd_fma(rffac, e1e2, sum_myRF);
		}

		const vcp_double_vec e1_x_e2_y_minus_e1_y_e2_x = vcp_simd_fms(eii_x, ejj_y, eii_y * ejj_x);
		const vcp_double_vec e1_y_e2_z_minus_e1_z_e2_y = vcp_simd_fms(eii_y, ejj_z, eii_z * ejj_y);
		const vcp_double_vec e1_z_e2_x_minus_e1_x_e2_z = vcp_simd_fms(eii_z, ejj_x, eii_x * ejj_z);

		M1_x = vcp_simd_fma(p1p2per4pie0r3, vcp_simd_fms(re2threeperr2, vcp_simd_fms(eii_y, dz, eii_z * dy), e1_y_e2_z_minus_e1_z_e2_y), rffac * e1_y_e2_z_minus_e1_z_e2_y);
		M1_y = vcp_simd_fma(p1p2per4pie0r3, vcp_simd_fms(re2threeperr2, vcp_simd_fms(eii_z, dx, eii_x * dz), e1_z_e2_x_minus_e1_x_e2_z), rffac * e1_z_e2_x_minus_e1_x_e2_z);
		M1_z = vcp_simd_fma(p1p2per4pie0r3, vcp_simd_fms(re2threeperr2, vcp_simd_fms(eii_x, dy, eii_y * dx), e1_x_e2_y_minus_e1_y_e2_x), rffac * e1_x_e2_y_minus_e1_y_e2_x);

		M2_x = vcp_simd_fms(p1p2per4pie0r3, vcp_simd_fma(re1threeperr2, vcp_simd_fms(ejj_y, dz, ejj_z * dy), e1_y_e2_z_minus_e1_z_e2_y), rffac * e1_y_e2_z_minus_e1_z_e2_y);
		M2_y = vcp_simd_fms(p1p2per4pie0r3, vcp_simd_fma(re1threeperr2, vcp_simd_fms(ejj_z, dx, ejj_x * dz), e1_z_e2_x_minus_e1_x_e2_z), rffac * e1_z_e2_x_minus_e1_x_e2_z);
		M2_z = vcp_simd_fms(p1p2per4pie0r3, vcp_simd_fma(re1threeperr2, vcp_simd_fms(ejj_x, dy, ejj_y * dx), e1_x_e2_y_minus_e1_y_e2_x), rffac * e1_x_e2_y_minus_e1_y_e2_x);
	}

	template<bool calculateMacroscopic>
	inline void VectorizedCellProcessor :: _loopBodyChargeQuadrupole(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& q,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
			const vcp_double_vec& m,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
			vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
			const vcp_mask_vec& forceMask) {

		const vcp_double_vec c_dx = r1_x - r2_x;
		const vcp_double_vec c_dy = r1_y - r2_y;
		const vcp_double_vec c_dz = r1_z - r2_z;//fma not possible since they will be reused...

		const vcp_double_vec c_dr2 = vcp_simd_scalProd(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const vcp_double_vec invdr2_unmasked = one / c_dr2;
		const vcp_double_vec invdr2 = vcp_simd_applymask(invdr2_unmasked, forceMask);
		const vcp_double_vec invdr = vcp_simd_sqrt(invdr2);

		const vcp_double_vec qQ05per4pie0 = _05 * q * m;

		vcp_double_vec costj = vcp_simd_scalProd(ejj_x, ejj_y, ejj_z, c_dx, c_dy, c_dz);
		costj = costj * invdr;

		const vcp_double_vec qQinv4dr3 = qQ05per4pie0 * invdr * invdr2;
		const vcp_double_vec part1 = three * costj * costj;
		const vcp_double_vec upot = qQinv4dr3 * (part1 - one);

		/**********
		 * Force
		 **********/
		const vcp_double_vec minus_partialRijInvdr = three * upot * invdr2;
		const vcp_double_vec partialTjInvdr = six * costj * qQinv4dr3 * invdr;

		const vcp_double_vec fac = vcp_simd_fma(costj * partialTjInvdr, invdr, minus_partialRijInvdr);

		f_x = vcp_simd_fms(fac, c_dx, partialTjInvdr * ejj_x);
		f_y = vcp_simd_fms(fac, c_dy, partialTjInvdr * ejj_y);
		f_z = vcp_simd_fms(fac, c_dz, partialTjInvdr * ejj_z);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		V_x = m_dx * f_x;
		V_y = m_dy * f_y;
		V_z = m_dz * f_z;

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {
			sum_upotXpoles = sum_upotXpoles + upot;

			sum_virial = sum_virial + V_x + V_y + V_z;
		}

		/**********
		 * Torque
		 **********/
		const vcp_double_vec minuseXrij_x = vcp_simd_fms(ejj_z, c_dy, ejj_y * c_dz);
		const vcp_double_vec minuseXrij_y = vcp_simd_fms(ejj_x, c_dz, ejj_z * c_dx);
		const vcp_double_vec minuseXrij_z = vcp_simd_fms(ejj_y, c_dx, ejj_x * c_dy);

		M_x = partialTjInvdr * minuseXrij_x;
		M_y = partialTjInvdr * minuseXrij_y;
		M_z = partialTjInvdr * minuseXrij_z;
	}

	template<bool calculateMacroscopic>
	inline void VectorizedCellProcessor :: _loopBodyDipoleQuadrupole(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
			const vcp_double_vec& p,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
			const vcp_double_vec& m,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
			vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
			vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
			const vcp_mask_vec& forceMask) {


		const vcp_double_vec c_dx = r1_x - r2_x;
		const vcp_double_vec c_dy = r1_y - r2_y;
		const vcp_double_vec c_dz = r1_z - r2_z;

		const vcp_double_vec c_dr2 = vcp_simd_scalProd(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const vcp_double_vec invdr2_unmasked = one / c_dr2;
		const vcp_double_vec invdr2 = vcp_simd_applymask(invdr2_unmasked, forceMask);
		const vcp_double_vec invdr = vcp_simd_sqrt(invdr2);

		const vcp_double_vec myqfac = _1pt5 * p * m * invdr2 * invdr2;

		vcp_double_vec costi = vcp_simd_scalProd(eii_x, eii_y, eii_z, c_dx, c_dy, c_dz);
		costi = costi * invdr;

		vcp_double_vec costj = vcp_simd_scalProd(ejj_x, ejj_y, ejj_z, c_dx, c_dy, c_dz);
		costj = costj * invdr;

		const vcp_double_vec cos2tj = costj * costj;

		const vcp_double_vec cosgij = vcp_simd_scalProd(eii_x, eii_y, eii_z, ejj_x, ejj_y, ejj_z);

		/************
		 * Potential
		 ************/
		// TODO: Check if upot has to be multiplied by -1 according to DISS_STOLL S.178.
		// This affects also the implementation in potforce.h
		const vcp_double_vec _5cos2tjminus1 = vcp_simd_fms(five, cos2tj, one);
		const vcp_double_vec _2costj = two * costj;

		vcp_double_vec part1 = costi * _5cos2tjminus1;
		vcp_double_vec part2 = _2costj * cosgij;

		vcp_double_vec const upot = myqfac * (part2 - part1);

		const vcp_double_vec myqfacXinvdr = myqfac * invdr;
		const vcp_double_vec minus_partialRijInvdr = four * upot * invdr2;
		const vcp_double_vec minus_partialTiInvdr = myqfacXinvdr * _5cos2tjminus1;

		part1 = vcp_simd_fms(five, costi * costj, cosgij); // *-1!

		const vcp_double_vec minus_partialTjInvdr = myqfacXinvdr * two * part1;
		const vcp_double_vec partialGij = myqfac * _2costj;

		vcp_double_vec part3 = vcp_simd_fma(costi, minus_partialTiInvdr, costj * minus_partialTjInvdr);
		const vcp_double_vec fac = vcp_simd_fnma(part3, invdr, minus_partialRijInvdr);//minus_pRI - part3*infdr

		// Force components
		f_x = vcp_simd_scalProd(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dx, eii_x, ejj_x);

		f_y = vcp_simd_scalProd(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dy, eii_y, ejj_y);

		f_z = vcp_simd_scalProd(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dz, eii_z, ejj_z);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		V_x = m_dx * f_x;
		V_y = m_dy * f_y;
		V_z = m_dz * f_z;

		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {
			sum_upotXpoles = sum_upotXpoles + upot;

			sum_virial = sum_virial + V_x + V_y + V_z;
		}

		/**********
		 * Torque
		 **********/
		const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = vcp_simd_fms(eii_x, ejj_y, eii_y * ejj_x);
		const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = vcp_simd_fms(eii_y, ejj_z, eii_z * ejj_y);
		const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = vcp_simd_fms(eii_z, ejj_x, eii_x * ejj_z);

		const vcp_double_vec partialGij_eiXej_x = partialGij * eii_y_ejj_z_minus_eii_z_ejj_y;
		const vcp_double_vec partialGij_eiXej_y = partialGij * eii_z_ejj_x_minus_eii_x_ejj_z;
		const vcp_double_vec partialGij_eiXej_z = partialGij * eii_x_ejj_y_minus_eii_y_ejj_x;

		vcp_double_vec eXrij_x = vcp_simd_fms(eii_y, c_dz, eii_z * c_dy);
		vcp_double_vec eXrij_y = vcp_simd_fms(eii_z, c_dx, eii_x * c_dz);
		vcp_double_vec eXrij_z = vcp_simd_fms(eii_x, c_dy, eii_y * c_dx);

		M1_x = vcp_simd_fms(minus_partialTiInvdr, eXrij_x, partialGij_eiXej_x);
		M1_y = vcp_simd_fms(minus_partialTiInvdr, eXrij_y, partialGij_eiXej_y);
		M1_z = vcp_simd_fms(minus_partialTiInvdr, eXrij_z, partialGij_eiXej_z);

		eXrij_x = vcp_simd_fms(ejj_y, c_dz, ejj_z * c_dy);
		eXrij_y = vcp_simd_fms(ejj_z, c_dx, ejj_x * c_dz);
		eXrij_z = vcp_simd_fms(ejj_x, c_dy, ejj_y * c_dx);

		M2_x = vcp_simd_fma(minus_partialTjInvdr, eXrij_x, partialGij_eiXej_x);
		M2_y = vcp_simd_fma(minus_partialTjInvdr, eXrij_y, partialGij_eiXej_y);
		M2_z = vcp_simd_fma(minus_partialTjInvdr, eXrij_z, partialGij_eiXej_z);
	}

	template<bool calculateMacroscopic>
	inline void VectorizedCellProcessor :: _loopBodyQuadrupole(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
			const vcp_double_vec& mii,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
			const vcp_double_vec& mjj,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& Mii_x, vcp_double_vec& Mii_y, vcp_double_vec& Mii_z,
			vcp_double_vec& Mjj_x, vcp_double_vec& Mjj_y, vcp_double_vec& Mjj_z,
			vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
			const vcp_mask_vec& forceMask)
	{
		const vcp_double_vec c_dx = r1_x - r2_x;
		const vcp_double_vec c_dy = r1_y - r2_y;
		const vcp_double_vec c_dz = r1_z - r2_z;

		const vcp_double_vec c_dr2 = vcp_simd_scalProd(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const vcp_double_vec invdr2_unmasked = one / c_dr2;
		const vcp_double_vec invdr2 = vcp_simd_applymask(invdr2_unmasked, forceMask);
		const vcp_double_vec invdr = vcp_simd_sqrt(invdr2);

		vcp_double_vec qfac = _075 * invdr;
		qfac = qfac * (mii * mjj);
		qfac = qfac * (invdr2 * invdr2);

		vcp_double_vec costi = vcp_simd_scalProd(eii_x, eii_y, eii_z, c_dx, c_dy, c_dz);
		costi = costi * invdr;

		vcp_double_vec costj = vcp_simd_scalProd(ejj_x, ejj_y, ejj_z, c_dx, c_dy, c_dz);
		costj = costj * invdr;

		const vcp_double_vec cos2ti = costi * costi;
		const vcp_double_vec cos2tj = costj * costj;

		const vcp_double_vec cosgij = vcp_simd_scalProd(eii_x, eii_y, eii_z, ejj_x, ejj_y, ejj_z);

		vcp_double_vec term = five * (costi * costj);
		term = cosgij - term;

		/************
		 * Potential
		 ************/
		vcp_double_vec part2 = _15 * cos2ti * cos2tj;
		vcp_double_vec part3 = two * term * term;
		vcp_double_vec upot = vcp_simd_fma(five, (cos2ti + cos2tj), part2);
		upot = (one + part3) - upot;
		upot = qfac * upot;

		/**********
		 * Force
		 **********/
		const vcp_double_vec minus_partialRijInvdr = five * upot * invdr2;

		// partialTiInvdr & partialTjInvdr
		vcp_double_vec part1 = qfac * ten * invdr;
		part2 = two * term;

		// partialTiInvdr only
		part3 = three * costi * cos2tj;
		vcp_double_vec part4 = costi + vcp_simd_fma(part2, costj, part3);
		const vcp_double_vec minus_partialTiInvdr = part1 * part4;

		// partialTjInvdr only
		part3 = three * costj * cos2ti;
		part4 = costj + vcp_simd_fma(part2, costi, part3);
		const vcp_double_vec minus_partialTjInvdr = part1 * part4;

		const vcp_double_vec partialGij = qfac * four * term;

		// fac
		part1 = minus_partialTiInvdr * costi;
		part2 = minus_partialTjInvdr * costj;
		//part3 = (part1 + part2) * invdr;
		const vcp_double_vec fac = vcp_simd_fnma((part1 + part2), invdr, minus_partialRijInvdr); //min - part3 = -part3 + min

		// Force components


		f_x = vcp_simd_scalProd(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dx, eii_x, ejj_x);

		f_y = vcp_simd_scalProd(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dy, eii_y, ejj_y);

		f_z = vcp_simd_scalProd(fac, minus_partialTiInvdr, minus_partialTjInvdr, c_dz, eii_z, ejj_z);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		V_x = m_dx * f_x;
		V_y = m_dy * f_y;
		V_z = m_dz * f_z;

		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (calculateMacroscopic) {
			sum_upotXpoles = sum_upotXpoles + upot;


			sum_virial = sum_virial + V_x + V_y + V_z;
		}

		/**********
		 * Torque
		 **********/

		const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = vcp_simd_fms(eii_x, ejj_y, eii_y * ejj_x);
		const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = vcp_simd_fms(eii_y, ejj_z, eii_z * ejj_y);
		const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = vcp_simd_fms(eii_z, ejj_x, eii_x * ejj_z);

		const vcp_double_vec partialGij_eiXej_x = partialGij * eii_y_ejj_z_minus_eii_z_ejj_y;
		const vcp_double_vec partialGij_eiXej_y = partialGij * eii_z_ejj_x_minus_eii_x_ejj_z;
		const vcp_double_vec partialGij_eiXej_z = partialGij * eii_x_ejj_y_minus_eii_y_ejj_x;

		vcp_double_vec eXrij_x = vcp_simd_fms(eii_y, c_dz, eii_z * c_dy);
		vcp_double_vec eXrij_y = vcp_simd_fms(eii_z, c_dx, eii_x * c_dz);
		vcp_double_vec eXrij_z = vcp_simd_fms(eii_x, c_dy, eii_y * c_dx);

		Mii_x = vcp_simd_fms(minus_partialTiInvdr, eXrij_x, partialGij_eiXej_x);
		Mii_y = vcp_simd_fms(minus_partialTiInvdr, eXrij_y, partialGij_eiXej_y);
		Mii_z = vcp_simd_fms(minus_partialTiInvdr, eXrij_z, partialGij_eiXej_z);

		eXrij_x = vcp_simd_fms(ejj_y, c_dz, ejj_z * c_dy);
		eXrij_y = vcp_simd_fms(ejj_z, c_dx, ejj_x * c_dz);
		eXrij_z = vcp_simd_fms(ejj_x, c_dy, ejj_y * c_dx);

		Mjj_x = vcp_simd_fma(minus_partialTjInvdr, eXrij_x, partialGij_eiXej_x);
		Mjj_y = vcp_simd_fma(minus_partialTjInvdr, eXrij_y, partialGij_eiXej_y);
		Mjj_z = vcp_simd_fma(minus_partialTjInvdr, eXrij_z, partialGij_eiXej_z);
	}

template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
void VectorizedCellProcessor::_calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
	const int tid = omp_get_thread_num();
	VLJCPThreadData &my_threadData = *_threadData[tid];

	// initialize dist lookups
	if(my_threadData._centers_dist_lookup.get_size() < soa2._centers_size){
		soa2.resizeCentersZero(my_threadData._centers_dist_lookup, soa2._centers_size);
	}
	soa2.initDistLookupPointers(my_threadData._centers_dist_lookup, my_threadData._ljc_dist_lookup, my_threadData._charges_dist_lookup, my_threadData._dipoles_dist_lookup, my_threadData._quadrupoles_dist_lookup);

	// Pointer for molecules
	const double * const soa1_mol_pos_x = soa1._mol_pos.xBegin();
	const double * const soa1_mol_pos_y = soa1._mol_pos.yBegin();
	const double * const soa1_mol_pos_z = soa1._mol_pos.zBegin();

	// Pointer for LJ centers
	const double * const soa1_ljc_r_x = soa1.ljc_r_xBegin();
	const double * const soa1_ljc_r_y = soa1.ljc_r_yBegin();
	const double * const soa1_ljc_r_z = soa1.ljc_r_zBegin();
	      double * const soa1_ljc_f_x = soa1.ljc_f_xBegin();
	      double * const soa1_ljc_f_y = soa1.ljc_f_yBegin();
	      double * const soa1_ljc_f_z = soa1.ljc_f_zBegin();
	      double * const soa1_ljc_V_x = soa1.ljc_V_xBegin();
	      double * const soa1_ljc_V_y = soa1.ljc_V_yBegin();
	      double * const soa1_ljc_V_z = soa1.ljc_V_zBegin();
	const int * const soa1_mol_ljc_num = soa1._mol_ljc_num;
	const size_t * const soa1_ljc_id = soa1._ljc_id;

	const double * const soa2_ljc_m_r_x = soa2.ljc_m_r_xBegin();
	const double * const soa2_ljc_m_r_y = soa2.ljc_m_r_yBegin();
	const double * const soa2_ljc_m_r_z = soa2.ljc_m_r_zBegin();
	const double * const soa2_ljc_r_x = soa2.ljc_r_xBegin();
	const double * const soa2_ljc_r_y = soa2.ljc_r_yBegin();
	const double * const soa2_ljc_r_z = soa2.ljc_r_zBegin();
	      double * const soa2_ljc_f_x = soa2.ljc_f_xBegin();
	      double * const soa2_ljc_f_y = soa2.ljc_f_yBegin();
	      double * const soa2_ljc_f_z = soa2.ljc_f_zBegin();
	      double * const soa2_ljc_V_x = soa2.ljc_V_xBegin();
	      double * const soa2_ljc_V_y = soa2.ljc_V_yBegin();
	      double * const soa2_ljc_V_z = soa2.ljc_V_zBegin();
	const size_t * const soa2_ljc_id = soa2._ljc_id;

	vcp_lookupOrMask_single* const soa2_ljc_dist_lookup = my_threadData._ljc_dist_lookup;

	// Pointer for charges
	const double * const soa1_charges_r_x = soa1.charges_r_xBegin();
	const double * const soa1_charges_r_y = soa1.charges_r_yBegin();
	const double * const soa1_charges_r_z = soa1.charges_r_zBegin();
	      double * const soa1_charges_f_x = soa1.charges_f_xBegin();
	      double * const soa1_charges_f_y = soa1.charges_f_yBegin();
	      double * const soa1_charges_f_z = soa1.charges_f_zBegin();
	      double * const soa1_charges_V_x = soa1.charges_V_xBegin();
	      double * const soa1_charges_V_y = soa1.charges_V_yBegin();
	      double * const soa1_charges_V_z = soa1.charges_V_zBegin();
	const double * const soa1_charges_q = soa1._charges_q;
	const int * const soa1_mol_charges_num = soa1._mol_charges_num;

	const double * const soa2_charges_m_r_x = soa2.charges_m_r_xBegin();
	const double * const soa2_charges_m_r_y = soa2.charges_m_r_yBegin();
	const double * const soa2_charges_m_r_z = soa2.charges_m_r_zBegin();
	const double * const soa2_charges_r_x   = soa2.charges_r_xBegin();
	const double * const soa2_charges_r_y   = soa2.charges_r_yBegin();
	const double * const soa2_charges_r_z   = soa2.charges_r_zBegin();
	      double * const soa2_charges_f_x   = soa2.charges_f_xBegin();
	      double * const soa2_charges_f_y   = soa2.charges_f_yBegin();
	      double * const soa2_charges_f_z   = soa2.charges_f_zBegin();
	      double * const soa2_charges_V_x   = soa2.charges_V_xBegin();
	      double * const soa2_charges_V_y   = soa2.charges_V_yBegin();
	      double * const soa2_charges_V_z   = soa2.charges_V_zBegin();
	const double * const soa2_charges_q = soa2._charges_q;

	vcp_lookupOrMask_single* const soa2_charges_dist_lookup = my_threadData._charges_dist_lookup;

	// Pointer for dipoles
	const double * const soa1_dipoles_r_x = soa1.dipoles_r_xBegin();
	const double * const soa1_dipoles_r_y = soa1.dipoles_r_yBegin();
	const double * const soa1_dipoles_r_z = soa1.dipoles_r_zBegin();
	      double * const soa1_dipoles_f_x = soa1.dipoles_f_xBegin();
	      double * const soa1_dipoles_f_y = soa1.dipoles_f_yBegin();
	      double * const soa1_dipoles_f_z = soa1.dipoles_f_zBegin();
	      double * const soa1_dipoles_V_x = soa1.dipoles_V_xBegin();
	      double * const soa1_dipoles_V_y = soa1.dipoles_V_yBegin();
	      double * const soa1_dipoles_V_z = soa1.dipoles_V_zBegin();
	const double * const soa1_dipoles_p = soa1._dipoles_p;
	const double * const soa1_dipoles_e_x = soa1._dipoles_e.xBegin();
	const double * const soa1_dipoles_e_y = soa1._dipoles_e.yBegin();
	const double * const soa1_dipoles_e_z = soa1._dipoles_e.zBegin();
	double * const soa1_dipoles_M_x = soa1._dipoles_M.xBegin();
	double * const soa1_dipoles_M_y = soa1._dipoles_M.yBegin();
	double * const soa1_dipoles_M_z = soa1._dipoles_M.zBegin();
	const int * const soa1_mol_dipoles_num = soa1._mol_dipoles_num;

	const double * const soa2_dipoles_m_r_x = soa2.dipoles_m_r_xBegin();
	const double * const soa2_dipoles_m_r_y = soa2.dipoles_m_r_yBegin();
	const double * const soa2_dipoles_m_r_z = soa2.dipoles_m_r_zBegin();
	const double * const soa2_dipoles_r_x   = soa2.dipoles_r_xBegin();
	const double * const soa2_dipoles_r_y   = soa2.dipoles_r_yBegin();
	const double * const soa2_dipoles_r_z   = soa2.dipoles_r_zBegin();
	      double * const soa2_dipoles_f_x   = soa2.dipoles_f_xBegin();
	      double * const soa2_dipoles_f_y   = soa2.dipoles_f_yBegin();
	      double * const soa2_dipoles_f_z   = soa2.dipoles_f_zBegin();
	      double * const soa2_dipoles_V_x   = soa2.dipoles_V_xBegin();
	      double * const soa2_dipoles_V_y   = soa2.dipoles_V_yBegin();
	      double * const soa2_dipoles_V_z   = soa2.dipoles_V_zBegin();
	const double * const soa2_dipoles_p = soa2._dipoles_p;
	const double * const soa2_dipoles_e_x = soa2._dipoles_e.xBegin();
	const double * const soa2_dipoles_e_y = soa2._dipoles_e.yBegin();
	const double * const soa2_dipoles_e_z = soa2._dipoles_e.zBegin();
	double * const soa2_dipoles_M_x = soa2._dipoles_M.xBegin();
	double * const soa2_dipoles_M_y = soa2._dipoles_M.yBegin();
	double * const soa2_dipoles_M_z = soa2._dipoles_M.zBegin();

	vcp_lookupOrMask_single* const soa2_dipoles_dist_lookup = my_threadData._dipoles_dist_lookup;

	// Pointer for quadrupoles
	const double * const soa1_quadrupoles_r_x = soa1.quadrupoles_r_xBegin();
	const double * const soa1_quadrupoles_r_y = soa1.quadrupoles_r_yBegin();
	const double * const soa1_quadrupoles_r_z = soa1.quadrupoles_r_zBegin();
	      double * const soa1_quadrupoles_f_x = soa1.quadrupoles_f_xBegin();
	      double * const soa1_quadrupoles_f_y = soa1.quadrupoles_f_yBegin();
	      double * const soa1_quadrupoles_f_z = soa1.quadrupoles_f_zBegin();
	      double * const soa1_quadrupoles_V_x = soa1.quadrupoles_V_xBegin();
	      double * const soa1_quadrupoles_V_y = soa1.quadrupoles_V_yBegin();
	      double * const soa1_quadrupoles_V_z = soa1.quadrupoles_V_zBegin();
	const double * const soa1_quadrupoles_m = soa1._quadrupoles_m;
	const double * const soa1_quadrupoles_e_x = soa1._quadrupoles_e.xBegin();
	const double * const soa1_quadrupoles_e_y = soa1._quadrupoles_e.yBegin();
	const double * const soa1_quadrupoles_e_z = soa1._quadrupoles_e.zBegin();
	      double * const soa1_quadrupoles_M_x = soa1._quadrupoles_M.xBegin();
	      double * const soa1_quadrupoles_M_y = soa1._quadrupoles_M.yBegin();
	      double * const soa1_quadrupoles_M_z = soa1._quadrupoles_M.zBegin();
	const int * const soa1_mol_quadrupoles_num = soa1._mol_quadrupoles_num;

	const double * const soa2_quadrupoles_m_r_x = soa2.quadrupoles_m_r_xBegin();
	const double * const soa2_quadrupoles_m_r_y = soa2.quadrupoles_m_r_yBegin();
	const double * const soa2_quadrupoles_m_r_z = soa2.quadrupoles_m_r_zBegin();
	const double * const soa2_quadrupoles_r_x   = soa2.quadrupoles_r_xBegin();
	const double * const soa2_quadrupoles_r_y   = soa2.quadrupoles_r_yBegin();
	const double * const soa2_quadrupoles_r_z   = soa2.quadrupoles_r_zBegin();
	      double * const soa2_quadrupoles_f_x   = soa2.quadrupoles_f_xBegin();
	      double * const soa2_quadrupoles_f_y   = soa2.quadrupoles_f_yBegin();
	      double * const soa2_quadrupoles_f_z   = soa2.quadrupoles_f_zBegin();
	      double * const soa2_quadrupoles_V_x   = soa2.quadrupoles_V_xBegin();
	      double * const soa2_quadrupoles_V_y   = soa2.quadrupoles_V_yBegin();
	      double * const soa2_quadrupoles_V_z   = soa2.quadrupoles_V_zBegin();
	const double * const soa2_quadrupoles_m = soa2._quadrupoles_m;
	const double * const soa2_quadrupoles_e_x = soa2._quadrupoles_e.xBegin();
	const double * const soa2_quadrupoles_e_y = soa2._quadrupoles_e.yBegin();
	const double * const soa2_quadrupoles_e_z = soa2._quadrupoles_e.zBegin();
	      double * const soa2_quadrupoles_M_x = soa2._quadrupoles_M.xBegin();
	      double * const soa2_quadrupoles_M_y = soa2._quadrupoles_M.yBegin();
	      double * const soa2_quadrupoles_M_z = soa2._quadrupoles_M.zBegin();

	vcp_lookupOrMask_single* const soa2_quadrupoles_dist_lookup = my_threadData._quadrupoles_dist_lookup;




	vcp_double_vec sum_upot6lj = VCP_SIMD_ZEROV;
	vcp_double_vec sum_upotXpoles = VCP_SIMD_ZEROV;
	vcp_double_vec sum_virial = VCP_SIMD_ZEROV;
	vcp_double_vec sum_myRF = VCP_SIMD_ZEROV;

	const vcp_double_vec rc2 = vcp_simd_set1(_LJCutoffRadiusSquare);
	const vcp_double_vec cutoffRadiusSquare = vcp_simd_set1(_cutoffRadiusSquare);
	const vcp_double_vec epsRFInvrc3 = vcp_simd_broadcast(&_epsRFInvrc3);

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

	//if(soa1._mol_num < 8){
	//	printf("less than 8\n");
	//}

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._mol_num; ++i) {//over the molecules
		const vcp_double_vec m1_r_x = vcp_simd_broadcast(soa1_mol_pos_x + i);
		const vcp_double_vec m1_r_y = vcp_simd_broadcast(soa1_mol_pos_y + i);
		const vcp_double_vec m1_r_z = vcp_simd_broadcast(soa1_mol_pos_z + i);
		// Iterate over centers of second cell
		const countertype32 compute_molecule_ljc = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_ljc_idx, soa2._ljc_num, _LJCutoffRadiusSquare,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);
		const countertype32 compute_molecule_charges = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_charge_idx, soa2._charges_num, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);
		const countertype32 compute_molecule_dipoles = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_dipole_idx, soa2._dipoles_num, _cutoffRadiusSquare,
				soa2_dipoles_dist_lookup, soa2_dipoles_m_r_x, soa2_dipoles_m_r_y, soa2_dipoles_m_r_z,
				cutoffRadiusSquare,	end_dipoles_j, m1_r_x, m1_r_y, m1_r_z);
		const countertype32 compute_molecule_quadrupoles = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_quadrupole_idx, soa2._quadrupoles_num, _cutoffRadiusSquare,
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
				vcp_double_vec sum_fx1 = VCP_SIMD_ZEROV;
				vcp_double_vec sum_fy1 = VCP_SIMD_ZEROV;
				vcp_double_vec sum_fz1 = VCP_SIMD_ZEROV;

				vcp_double_vec sum_Vx1 = VCP_SIMD_ZEROV;
				vcp_double_vec sum_Vy1 = VCP_SIMD_ZEROV;
				vcp_double_vec sum_Vz1 = VCP_SIMD_ZEROV;

				const vcp_double_vec c_r_x1 = vcp_simd_broadcast(soa1_ljc_r_x + i_ljc_idx);
				const vcp_double_vec c_r_y1 = vcp_simd_broadcast(soa1_ljc_r_y + i_ljc_idx);
				const vcp_double_vec c_r_z1 = vcp_simd_broadcast(soa1_ljc_r_z + i_ljc_idx);

				// Iterate over each pair of centers in the second cell.
				size_t j = ForcePolicy::InitJ2(i_ljc_idx);
				for (; j < end_ljc_loop; j += VCP_VEC_SIZE) {//over (all/some) lj-centers -- amount depends on ForcePolicy::InitJ (two cells: all, within one cell: only do i,j not j,i)
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_ljc_dist_lookup, j);
					// Only go on if at least 1 of the forces has to be calculated.
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const vcp_double_vec c_r_x2 = MaskGatherChooser::load(soa2_ljc_r_x, j, lookupORforceMask);
						const vcp_double_vec c_r_y2 = MaskGatherChooser::load(soa2_ljc_r_y, j, lookupORforceMask);
						const vcp_double_vec c_r_z2 = MaskGatherChooser::load(soa2_ljc_r_z, j, lookupORforceMask);

						const vcp_double_vec m_r_x2 = MaskGatherChooser::load(soa2_ljc_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m_r_y2 = MaskGatherChooser::load(soa2_ljc_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m_r_z2 = MaskGatherChooser::load(soa2_ljc_m_r_z, j, lookupORforceMask);

						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						vcp_double_vec fx, fy, fz;
						vcp_double_vec Vx, Vy, Vz;

						vcp_double_vec eps_24;
						vcp_double_vec sig2;
						unpackEps24Sig2<MaskGatherChooser>(eps_24, sig2, _eps_sig[id_i], soa2_ljc_id, j, lookupORforceMask);

						vcp_double_vec shift6;
						unpackShift6<MaskGatherChooser>(shift6, _shift6[id_i], soa2_ljc_id, j, lookupORforceMask);

						_loopBodyLJ<CalculateMacroscopic>(
							m1_r_x, m1_r_y, m1_r_z, c_r_x1, c_r_y1, c_r_z1,
							m_r_x2, m_r_y2, m_r_z2, c_r_x2, c_r_y2, c_r_z2,
							fx, fy, fz,
							Vx, Vy, Vz,
							sum_upot6lj, sum_virial,
							MaskGatherChooser::getForceMask(lookupORforceMask),
							eps_24, sig2,
							shift6);

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_ljc_f_x, j, fx, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_ljc_f_y, j, fy, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_ljc_f_z, j, fz, lookupORforceMask);

						sum_fx1 = sum_fx1 + fx;
						sum_fy1 = sum_fy1 + fy;
						sum_fz1 = sum_fz1 + fz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_ljc_V_x, j, Vx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_ljc_V_y, j, Vy, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_ljc_V_z, j, Vz, lookupORforceMask);

						sum_Vx1 = sum_Vx1 + Vx;
						sum_Vy1 = sum_Vy1 + Vy;
						sum_Vz1 = sum_Vz1 + Vz;
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_ljc);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_ljc_dist_lookup, j, remainderM);

						const vcp_double_vec c_r_x2 = MaskGatherChooser::load(soa2_ljc_r_x, j, lookupORforceMask);
						const vcp_double_vec c_r_y2 = MaskGatherChooser::load(soa2_ljc_r_y, j, lookupORforceMask);
						const vcp_double_vec c_r_z2 = MaskGatherChooser::load(soa2_ljc_r_z, j, lookupORforceMask);

						const vcp_double_vec m_r_x2 = MaskGatherChooser::load(soa2_ljc_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m_r_y2 = MaskGatherChooser::load(soa2_ljc_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m_r_z2 = MaskGatherChooser::load(soa2_ljc_m_r_z, j, lookupORforceMask);

						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						vcp_double_vec fx, fy, fz;
						vcp_double_vec Vx, Vy, Vz;

						vcp_double_vec eps_24;
						vcp_double_vec sig2;
						unpackEps24Sig2<MaskGatherChooser>(eps_24, sig2, _eps_sig[id_i], soa2_ljc_id, j, lookupORforceMask);

						vcp_double_vec shift6;
						unpackShift6<MaskGatherChooser>(shift6, _shift6[id_i], soa2_ljc_id, j, lookupORforceMask);

						_loopBodyLJ<CalculateMacroscopic>(
							m1_r_x, m1_r_y, m1_r_z, c_r_x1, c_r_y1, c_r_z1,
							m_r_x2, m_r_y2, m_r_z2, c_r_x2, c_r_y2, c_r_z2,
							fx, fy, fz,
							Vx, Vy, Vz,
							sum_upot6lj, sum_virial,
							remainderM,//use remainder mask as forcemask
							eps_24, sig2,
							shift6);

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_ljc_f_x, j, fx, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_ljc_f_y, j, fy, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_ljc_f_z, j, fz, lookupORforceMask, remainderM);

						sum_fx1 = sum_fx1 + fx;
						sum_fy1 = sum_fy1 + fy;
						sum_fz1 = sum_fz1 + fz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_ljc_V_x, j, Vx, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_ljc_V_y, j, Vy, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_ljc_V_z, j, Vz, lookupORforceMask, remainderM);


						sum_Vx1 = sum_Vx1 + Vx;
						sum_Vy1 = sum_Vy1 + Vy;
						sum_Vz1 = sum_Vz1 +Vz;

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

				const vcp_double_vec q1 = vcp_simd_broadcast(soa1_charges_q + i_charge_idx + local_i);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_charges_r_x + i_charge_idx + local_i);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_charges_r_y + i_charge_idx + local_i);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_charges_r_z + i_charge_idx + local_i);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_charge_idx + local_i);

				for (; j < end_charges_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_charges_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const vcp_double_vec q2 = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));



						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask);

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask);
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_charges_dist_lookup, j, remainderM);

						const vcp_double_vec q2 = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);
						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								sum_upotXpoles, sum_virial,
								remainderM);



						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask, remainderM);

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
				const vcp_double_vec p = vcp_simd_broadcast(soa1_dipoles_p + i_dipole_charge_idx);
				const vcp_double_vec e_x = vcp_simd_broadcast(soa1_dipoles_e_x + i_dipole_charge_idx);
				const vcp_double_vec e_y = vcp_simd_broadcast(soa1_dipoles_e_y + i_dipole_charge_idx);
				const vcp_double_vec e_z = vcp_simd_broadcast(soa1_dipoles_e_z + i_dipole_charge_idx);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_dipoles_r_x + i_dipole_charge_idx);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_dipoles_r_y + i_dipole_charge_idx);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_dipoles_r_z + i_dipole_charge_idx);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_M_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M_z = VCP_SIMD_ZEROV;

				size_t j = ForcePolicy::InitJ2(i_charge_idx);
				for (; j < end_charges_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_charges_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const vcp_double_vec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						// Store forces

						sum_f1_x = sum_f1_x - f_x;//negative, since dipole-charge, not charge-dipole -> direction inversed
						sum_f1_y = sum_f1_y - f_y;
						sum_f1_z = sum_f1_z - f_z;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask);//newton 3

						//store virials
						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask);//newton 3

						// Store torque

						sum_M_x = vcp_simd_add(sum_M_x, M_x);
						sum_M_y = vcp_simd_add(sum_M_y, M_y);
						sum_M_z = vcp_simd_add(sum_M_z, M_z);
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_charges_dist_lookup, j, remainderM);
						const vcp_double_vec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						// Store forces

						sum_f1_x = sum_f1_x - f_x;
						sum_f1_y = sum_f1_y - f_y;
						sum_f1_z = sum_f1_z - f_z;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3

						// Store torque

						sum_M_x = vcp_simd_add(sum_M_x, M_x);
						sum_M_y = vcp_simd_add(sum_M_y, M_y);
						sum_M_z = vcp_simd_add(sum_M_z, M_z);
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
				const vcp_double_vec m = vcp_simd_broadcast(soa1_quadrupoles_m + i_quadrupole_charge_idx);
				const vcp_double_vec e_x = vcp_simd_broadcast(soa1_quadrupoles_e_x + i_quadrupole_charge_idx);
				const vcp_double_vec e_y = vcp_simd_broadcast(soa1_quadrupoles_e_y + i_quadrupole_charge_idx);
				const vcp_double_vec e_z = vcp_simd_broadcast(soa1_quadrupoles_e_z + i_quadrupole_charge_idx);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_quadrupoles_r_x + i_quadrupole_charge_idx);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_quadrupoles_r_y + i_quadrupole_charge_idx);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_quadrupoles_r_z + i_quadrupole_charge_idx);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_M1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_z = VCP_SIMD_ZEROV;

				size_t j = ForcePolicy::InitJ2(i_charge_idx);
				for (; j < end_charges_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_charges_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const vcp_double_vec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						// Store forces

						sum_f1_x = sum_f1_x - f_x;
						sum_f1_y = sum_f1_y - f_y;
						sum_f1_z = sum_f1_z - f_z;


						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask);//newton 3


						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M_z);
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_charges_dist_lookup, j, remainderM);
						const vcp_double_vec q = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						// Store forces

						sum_f1_x = sum_f1_x - f_x;
						sum_f1_y = sum_f1_y - f_y;
						sum_f1_z = sum_f1_z - f_z;


						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;


						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M_z);

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

				const vcp_double_vec p1 = vcp_simd_broadcast(soa1_dipoles_p + i_dipole_idx + local_i);
				const vcp_double_vec e1_x = vcp_simd_broadcast(soa1_dipoles_e_x + i_dipole_idx + local_i);
				const vcp_double_vec e1_y = vcp_simd_broadcast(soa1_dipoles_e_y + i_dipole_idx + local_i);
				const vcp_double_vec e1_z = vcp_simd_broadcast(soa1_dipoles_e_z + i_dipole_idx + local_i);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_dipoles_r_x + i_dipole_idx + local_i);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_dipoles_r_y + i_dipole_idx + local_i);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_dipoles_r_z + i_dipole_idx + local_i);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_M1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_z = VCP_SIMD_ZEROV;

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_dipole_idx + local_i);
				for (; j < end_dipoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_dipoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const vcp_double_vec p2 = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const vcp_double_vec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

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

						// Store forces

						sum_f1_x = vcp_simd_add(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_add(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_add(sum_f1_z, f_z);

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_z, j, f_z, lookupORforceMask);//newton 3

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
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_dipoles_dist_lookup, j, remainderM);
						const vcp_double_vec p2 = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const vcp_double_vec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

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

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

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

				const vcp_double_vec q = vcp_simd_broadcast(soa1_charges_q + i_charge_dipole_idx);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_charges_r_x + i_charge_dipole_idx);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_charges_r_y + i_charge_dipole_idx);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_charges_r_z + i_charge_dipole_idx);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				size_t j = ForcePolicy::InitJ2(i_dipole_idx);
				for (; j < end_dipoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_dipoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const vcp_double_vec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);

						const vcp_double_vec e_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);

						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_dipoles_f_z, j, f_z, lookupORforceMask);//newton 3

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
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_dipoles_dist_lookup, j, remainderM);
						const vcp_double_vec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);

						const vcp_double_vec e_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);

						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeDipole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_dipoles_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

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

				const vcp_double_vec m = vcp_simd_broadcast(soa1_quadrupoles_m + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_x = vcp_simd_broadcast(soa1_quadrupoles_e_x + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_y = vcp_simd_broadcast(soa1_quadrupoles_e_y + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_z = vcp_simd_broadcast(soa1_quadrupoles_e_z + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_quadrupoles_r_x + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_quadrupoles_r_y + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_quadrupoles_r_z + i_quadrupole_dipole_idx);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_M1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_z = VCP_SIMD_ZEROV;

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_dipole_idx);
				for (; j < end_dipoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_dipoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const vcp_double_vec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const vcp_double_vec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M2_x, M2_y, M2_z, M1_x, M1_y, M1_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						// Store forces

						sum_f1_x = sum_f1_x - f_x;
						sum_f1_y = sum_f1_y - f_y;
						sum_f1_z = sum_f1_z - f_z;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_f_z, j, f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask);//newton 3

						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);


						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_x, j, M2_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_y, j, M2_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_dipoles_M_z, j, M2_z, lookupORforceMask);//newton 3

					}
				}
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_dipoles);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_dipoles_dist_lookup, j, remainderM);
						const vcp_double_vec p = MaskGatherChooser::load(soa2_dipoles_p, j, lookupORforceMask);
						const vcp_double_vec e2_x = MaskGatherChooser::load(soa2_dipoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e2_y = MaskGatherChooser::load(soa2_dipoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e2_z = MaskGatherChooser::load(soa2_dipoles_e_z, j, lookupORforceMask);
						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_dipoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_dipoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_dipoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_dipoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_dipoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_dipoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M2_x, M2_y, M2_z, M1_x, M1_y, M1_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						// Store forces

						sum_f1_x = sum_f1_x - f_x;
						sum_f1_y = sum_f1_y - f_y;
						sum_f1_z = sum_f1_z - f_z;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_dipoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3



						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);


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
				const vcp_double_vec mii = vcp_simd_broadcast(soa1_quadrupoles_m + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_x = vcp_simd_broadcast(soa1_quadrupoles_e_x + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_y = vcp_simd_broadcast(soa1_quadrupoles_e_y + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_z = vcp_simd_broadcast(soa1_quadrupoles_e_z + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_x = vcp_simd_broadcast(soa1_quadrupoles_r_x + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_y = vcp_simd_broadcast(soa1_quadrupoles_r_y + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_z = vcp_simd_broadcast(soa1_quadrupoles_r_z + i_quadrupole_idx + local_i);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_M1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_z = VCP_SIMD_ZEROV;

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_quadrupole_idx + local_i);
				for (; j < end_quadrupoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_quadrupoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const vcp_double_vec mjj = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const vcp_double_vec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const vcp_double_vec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const vcp_double_vec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const vcp_double_vec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const vcp_double_vec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const vcp_double_vec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, mii,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, mjj,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_z, j, f_z, lookupORforceMask);//newton 3


						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask);//newton 3

						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M2_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M2_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M2_z, lookupORforceMask);//newton 3
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_quadrupoles_dist_lookup, j, remainderM);
						const vcp_double_vec mjj = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const vcp_double_vec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const vcp_double_vec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const vcp_double_vec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const vcp_double_vec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const vcp_double_vec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const vcp_double_vec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, mii,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, mjj,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

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
				const vcp_double_vec q = vcp_simd_broadcast(soa1_charges_q + i_charge_quadrupole_idx);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_charges_r_x + i_charge_quadrupole_idx);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_charges_r_y + i_charge_quadrupole_idx);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_charges_r_z + i_charge_quadrupole_idx);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				size_t j = ForcePolicy::InitJ2(i_quadrupole_idx);
				for (; j < end_quadrupoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_quadrupoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const vcp_double_vec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const vcp_double_vec e_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);


						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_z, j, f_z, lookupORforceMask);//newton 3

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
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_quadrupoles_dist_lookup, j, remainderM);
						const vcp_double_vec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const vcp_double_vec e_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const vcp_double_vec e_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const vcp_double_vec e_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const vcp_double_vec r2_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const vcp_double_vec r2_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const vcp_double_vec r2_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);


						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyChargeQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

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
				const vcp_double_vec p = vcp_simd_broadcast(soa1_dipoles_p + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_x = vcp_simd_broadcast(soa1_dipoles_e_x + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_y = vcp_simd_broadcast(soa1_dipoles_e_y + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_z = vcp_simd_broadcast(soa1_dipoles_e_z + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_x = vcp_simd_broadcast(soa1_dipoles_r_x + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_y = vcp_simd_broadcast(soa1_dipoles_r_y + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_z = vcp_simd_broadcast(soa1_dipoles_r_z + i_dipole_quadrupole_idx);

				vcp_double_vec sum_f1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_f1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_V1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_V1_z = VCP_SIMD_ZEROV;

				vcp_double_vec sum_M1_x = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_y = VCP_SIMD_ZEROV;
				vcp_double_vec sum_M1_z = VCP_SIMD_ZEROV;

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_quadrupole_idx);
				for (; j < end_quadrupoles_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_quadrupoles_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {

						const vcp_double_vec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const vcp_double_vec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const vcp_double_vec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const vcp_double_vec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const vcp_double_vec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const vcp_double_vec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const vcp_double_vec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, p,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_x, j, f_x, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_y, j, f_y, lookupORforceMask);//newton 3
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_quadrupoles_f_z, j, f_z, lookupORforceMask);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask);//newton 3


						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_x, j, M2_x, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_y, j, M2_y, lookupORforceMask);//newton 3
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_quadrupoles_M_z, j, M2_z, lookupORforceMask);//newton 3

					}
				}
#if VCP_VEC_TYPE == VCP_VEC_MIC_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_quadrupoles);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_quadrupoles_dist_lookup, j, remainderM);
						const vcp_double_vec m = MaskGatherChooser::load(soa2_quadrupoles_m, j, lookupORforceMask);
						const vcp_double_vec ejj_x = MaskGatherChooser::load(soa2_quadrupoles_e_x, j, lookupORforceMask);
						const vcp_double_vec ejj_y = MaskGatherChooser::load(soa2_quadrupoles_e_y, j, lookupORforceMask);
						const vcp_double_vec ejj_z = MaskGatherChooser::load(soa2_quadrupoles_e_z, j, lookupORforceMask);
						const vcp_double_vec rjj_x = MaskGatherChooser::load(soa2_quadrupoles_r_x, j, lookupORforceMask);
						const vcp_double_vec rjj_y = MaskGatherChooser::load(soa2_quadrupoles_r_y, j, lookupORforceMask);
						const vcp_double_vec rjj_z = MaskGatherChooser::load(soa2_quadrupoles_r_z, j, lookupORforceMask);

						const vcp_double_vec m2_r_x = MaskGatherChooser::load(soa2_quadrupoles_m_r_x, j, lookupORforceMask);
						const vcp_double_vec m2_r_y = MaskGatherChooser::load(soa2_quadrupoles_m_r_y, j, lookupORforceMask);
						const vcp_double_vec m2_r_z = MaskGatherChooser::load(soa2_quadrupoles_m_r_z, j, lookupORforceMask);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;
						vcp_double_vec Vx, Vy, Vz;

						_loopBodyDipoleQuadrupole<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, p,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, m,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								remainderM);

						// Store forces

						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_x, j, f_x, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_y, j, f_y, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_quadrupoles_f_z, j, f_z, lookupORforceMask, remainderM);//newton 3

						// Store virials

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_x, j, Vx, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_y, j, Vy, lookupORforceMask, remainderM);//newton 3
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_quadrupoles_V_z, j, Vz, lookupORforceMask, remainderM);//newton 3


						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

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

	hSum_Add_Store(my_threadData._upot6ljV, sum_upot6lj);
	hSum_Add_Store(my_threadData._upotXpolesV, sum_upotXpoles);
	hSum_Add_Store(my_threadData._virialV, sum_virial);
	hSum_Add_Store(my_threadData._myRFV, vcp_simd_sub(zero, sum_myRF));

} // void LennardJonesCellHandler::CalculatePairs_(LJSoA & soa1, LJSoA & soa2)

void VectorizedCellProcessor::processCell(ParticleCell & c) {
	CellDataSoA& soa = c.getCellDataSoA();
	if (c.isHaloCell() or soa._mol_num < 2) {
		return;
	}
	const bool CalculateMacroscopic = true;
	const bool ApplyCutoff = true;
	_calculatePairs<SingleCellPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa, soa);
}

// provisionally, the code from the legacy cell processor is used here
//
int VectorizedCellProcessor::countNeighbours(Molecule* m1, ParticleCell& cell2, double RR)
{
	int tn = 0;
	double distanceVector[3];

	int neighbourParticleCount = cell2.getMoleculeCount();

	for (int j = 0; j < neighbourParticleCount; j++) {
		Molecule& molecule2 = cell2.moleculesAt(j);
		if(m1->id() == molecule2.id()) continue;
		double dd = molecule2.dist2(*m1, distanceVector);
		if (dd < RR) tn++;
	}
	return tn;
}

void VectorizedCellProcessor::processCellPair(ParticleCell & c1, ParticleCell & c2) {
	assert(&c1 != &c2);
	const CellDataSoA& soa1 = c1.getCellDataSoA();
	const CellDataSoA& soa2 = c2.getCellDataSoA();
	const bool c1Halo = c1.isHaloCell();
	const bool c2Halo = c2.isHaloCell();

	// this variable determines whether
	// _calcPairs(soa1, soa2) or _calcPairs(soa2, soa1)
	// is more efficient
	const bool calc_soa1_soa2 = (soa1._mol_num <= soa2._mol_num);

	// if one cell is empty, or both cells are Halo, skip
	if (soa1._mol_num == 0 or soa2._mol_num == 0 or (c1Halo and c2Halo)) {
		return;
	}

	// Macroscopic conditions:
	// if none of the cells is halo, then compute
	// if one of them is halo:
	// 		if c1-index < c2-index, then compute
	// 		else, then don't compute
	// This saves the Molecule::isLessThan checks
	// and works similar to the "Half-Shell" scheme

	const bool ApplyCutoff = true;

	if ((not c1Halo and not c2Halo) or						// no cell is halo or
			(c1.getCellIndex() < c2.getCellIndex())) 		// one of them is halo, but c1.index < c2.index
	{
		const bool CalculateMacroscopic = true;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}

	} else {
		assert(c1Halo != c2Halo);							// one of them is halo and
		assert(not (c1.getCellIndex() < c2.getCellIndex()));// c1.index not < c2.index

		const bool CalculateMacroscopic = false;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}
	}
}
