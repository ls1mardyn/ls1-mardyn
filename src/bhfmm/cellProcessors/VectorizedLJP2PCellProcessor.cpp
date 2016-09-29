/**
 * \file
 * \brief VectorizedLJP2PCellProcessor.cpp
 */

#include "VectorizedLJP2PCellProcessor.h"

#include "particleContainer/adapter/CellDataSoA.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include <algorithm>
#include "particleContainer/adapter/vectorization/MaskGatherChooser.h"

using namespace Log;
using namespace std;
namespace bhfmm {
VectorizedLJP2PCellProcessor::VectorizedLJP2PCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		CellProcessor(cutoffRadius, LJcutoffRadius), _domain(domain),
		// maybe move the following to somewhere else:
		_epsRFInvrc3(2. * (domain.getepsilonRF() - 1.) / ((cutoffRadius * cutoffRadius * cutoffRadius) * (2. * domain.getepsilonRF() + 1.))), 
		_eps_sig(), _shift6(), _upot6lj(0.0), _virial(0.0), _ljc_dist_lookup(0) {

#if VCP_VEC_TYPE==VCP_NOVEC
	global_log->info() << "VectorizedLJP2PCellProcessor: using no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	global_log->info() << "VectorizedLJP2PCellProcessor: using SSE3 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	global_log->info() << "VectorizedLJP2PCellProcessor: using AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	global_log->info() << "VectorizedLJP2PCellProcessor: using AVX2 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_MIC
	global_log->info() << "VectorizedLJP2PCellProcessor: using MIC intrinsics." << std::endl;
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

#ifdef ENABLE_MPI
	_timer.set_sync(false);
#endif
}

VectorizedLJP2PCellProcessor :: ~VectorizedLJP2PCellProcessor () {
}

void VectorizedLJP2PCellProcessor::printTimers() {
	std::cout << "FMM: Time spent in LJ P2P " << _timer.get_etime() << std::endl;
}


void VectorizedLJP2PCellProcessor::initTraversal() {
	_timer.start();
	_virial = 0.0;
	_upot6lj = 0.0;

	global_log->debug() << "VectorizedLJCellProcessor::initTraversal()." << std::endl;
}


void VectorizedLJP2PCellProcessor::endTraversal() {
	_domain.setLocalVirial(_virial /*+ 3.0 * _myRF*/);
	_domain.setLocalUpot(_upot6lj / 6.0 /*+ _upotXpoles + _myRF*/);
	_timer.stop();
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
	void VectorizedLJP2PCellProcessor :: _loopBodyLJ(
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

template<class ForcePolicy, class MaskGatherChooser>
countertype32
inline VectorizedLJP2PCellProcessor::calcDistLookup (const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		vcp_lookupOrMask_single* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z,
		const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z) {

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	vcp_mask_vec initJ_mask = ForcePolicy :: InitJ_Mask(i_center_idx);

	MaskGatherChooser mgc(soa2_center_dist_lookup, j);

	for (; j < end_j; j += VCP_VEC_SIZE) {
		const vcp_double_vec m2_r_x = vcp_simd_load(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = vcp_simd_load(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = vcp_simd_load(soa2_m_r_z + j);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const vcp_mask_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);

		mgc.storeCalcDistLookup(j, forceMask);

	}
	const vcp_mask_vec remainderMask = vcp_simd_getRemainderMask(soa2_num_centers);
	if (vcp_simd_movemask(remainderMask)) {
		const vcp_double_vec m2_r_x = vcp_simd_maskload(soa2_m_r_x + j, remainderMask);
		const vcp_double_vec m2_r_y = vcp_simd_maskload(soa2_m_r_y + j, remainderMask);
		const vcp_double_vec m2_r_z = vcp_simd_maskload(soa2_m_r_z + j, remainderMask);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const vcp_mask_vec forceMask = vcp_simd_and(remainderMask, ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask));//AND remainderMask -> set unimportant ones to zero.
		mgc.storeCalcDistLookup(j, forceMask);
	}

	return mgc.getCount();	//do not compute stuff if nothing needs to be computed.

}

template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
void VectorizedLJP2PCellProcessor :: _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
	// initialize dist lookups
	if(_centers_dist_lookup.get_size() < soa2._ljc_size){
		soa2.resizeLastZero(_centers_dist_lookup, soa2._ljc_size, soa2._ljc_num);
	}
	_ljc_dist_lookup = _centers_dist_lookup;

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

	vcp_lookupOrMask_single* const soa2_ljc_dist_lookup = _ljc_dist_lookup;

	vcp_double_vec sum_upot6lj = VCP_SIMD_ZEROV;
	vcp_double_vec sum_virial = VCP_SIMD_ZEROV;

	const vcp_double_vec rc2 = vcp_simd_set1(_LJCutoffRadiusSquare);

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
	size_t i_ljc_idx = 0;

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

		size_t end_ljc_loop = MaskGatherChooser::getEndloop(end_ljc_j_longloop, compute_molecule_ljc);

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
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(end_ljc_j_cnt);
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

	}

	hSum_Add_Store(&_upot6lj, sum_upot6lj);
	hSum_Add_Store(&_virial, sum_virial);

}

void VectorizedLJP2PCellProcessor::processCell(ParticleCell & c) {
	CellDataSoA& soa = c.getCellDataSoA();
	if (c.isHaloCell() or soa._mol_num < 2) {
		return;
	}
	const bool CalculateMacroscopic = true;
	_calculatePairs<SingleCellPolicy_, CalculateMacroscopic, MaskGatherC>(soa, soa);
}

void VectorizedLJP2PCellProcessor::processCellPair(ParticleCell & c1, ParticleCell & c2) {
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

	if ((not c1Halo and not c2Halo) or						// no cell is halo or
			(c1.getCellIndex() < c2.getCellIndex())) 		// one of them is halo, but c1.index < c2.index
	{
		const bool CalculateMacroscopic = true;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}

	} else {
		assert(c1Halo != c2Halo);							// one of them is halo and
		assert(not (c1.getCellIndex() < c2.getCellIndex()));// c1.index not < c2.index

		const bool CalculateMacroscopic = false;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}
	}
}

} // namespace bhfmm

