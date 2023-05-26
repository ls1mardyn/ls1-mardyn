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
		_eps_sig(), _shift6(), _upot6lj(0.0), _virial(0.0){

#if VCP_VEC_TYPE==VCP_NOVEC
	global_log->info() << "VectorizedLJP2PCellProcessor: using no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	global_log->info() << "VectorizedLJP2PCellProcessor: using SSE3 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	global_log->info() << "VectorizedLJP2PCellProcessor: using AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	global_log->info() << "VectorizedLJP2PCellProcessor: using AVX2 intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_KNL) || (VCP_VEC_TYPE==VCP_VEC_KNL_GATHER)
	global_log->info() << "VectorizedLJP2PCellProcessor: using KNL intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_AVX512F) || (VCP_VEC_TYPE==VCP_VEC_AVX512F_GATHER)
	global_log->info() << "VectorizedLJP2PCellProcessor: using SKX intrinsics." << std::endl;
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
	global_log->info() << "VectorizedLJP2PCellProcessor: allocate data for " << _numThreads << " threads." << std::endl;
	_threadData.resize(_numThreads);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		VLJP2PCPThreadData * myown = new VLJP2PCPThreadData();
		const int myid = mardyn_get_thread_num();
		_threadData[myid] = myown;
	} // end pragma omp parallel

#ifdef ENABLE_MPI
	global_simulation->timers()->setOutputString("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P", "FMM: Time spent in LJ P2P ");
	//global_simulation->setSyncTimer("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P", false); //it is per default false
#endif
}

VectorizedLJP2PCellProcessor :: ~VectorizedLJP2PCellProcessor () {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int myid = mardyn_get_thread_num();
		delete _threadData[myid];
	}
}

void VectorizedLJP2PCellProcessor::printTimers() {
	std::cout << "FMM: Time spent in LJ P2P " << global_simulation->timers()->getTime("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P") << std::endl;
	global_simulation->timers()->print("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P");
}

void VectorizedLJP2PCellProcessor::initTraversal() {
	global_simulation->timers()->start("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P");

	#if defined(_OPENMP)
	#pragma omp master
	#endif
	{
		_upot6lj = 0.0;
		_virial = 0.0;
	} // end pragma omp master

	global_log->debug() << "VectorizedLJP2PCellProcessor::initTraversal()." << std::endl;
}

void VectorizedLJP2PCellProcessor::endTraversal() {
	vcp_real_accum glob_upot6lj = 0.0;
	vcp_real_accum glob_virial = 0.0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:glob_upot6lj, glob_virial)
	#endif
	{
		const int tid = mardyn_get_thread_num();

		// reduce vectors and clear local variable
		vcp_real_accum thread_upot = 0.0, thread_virial = 0.0;

		load_hSum_Store_Clear(&thread_upot, _threadData[tid]->_upot6ljV);
		load_hSum_Store_Clear(&thread_virial, _threadData[tid]->_virialV);

		// add to global sum
		glob_upot6lj += thread_upot;
		glob_virial += thread_virial;
	} // end pragma omp parallel reduction

	_upot6lj = glob_upot6lj;
	_virial = glob_virial;
    _domain.setLocalVirial(_virial /*+ 3.0 * _myRF*/ + _domain.getLocalVirial());
    _domain.setLocalUpot(_upot6lj / 6.0 /*+ _upotXpoles + _myRF*/ + _domain.getLocalUpot());
	global_simulation->timers()->stop("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P");
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
inline
void VectorizedLJP2PCellProcessor :: _loopBodyLJ(
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
	const RealCalcVec r2_inv_unmasked = one / c_r2;
	const RealCalcVec r2_inv = RealCalcVec::apply_mask(r2_inv_unmasked, forceMask);


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

		const RealCalcVec upot_sh = RealCalcVec::fmadd(eps_24, lj12m6, shift6); //2 FP upot	//shift6 is not masked -> we have to mask upot_shifted
		const RealCalcVec upot_masked = RealCalcVec::apply_mask(upot_sh, forceMask); //mask it
		const RealAccumVec upot_masked_accum = RealAccumVec::convertCalcToAccum(upot_masked);

		sum_upot6lj = sum_upot6lj + upot_masked_accum;//1FP (sum macro)

		sum_virial = sum_virial +  V_x + V_y + V_z;//1 FP (sum macro) + 2 FP (virial)
	}
}

template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
void VectorizedLJP2PCellProcessor::_calculatePairs(CellDataSoA & soa1, CellDataSoA & soa2) {
	const int tid = mardyn_get_thread_num();
	VLJP2PCPThreadData &my_threadData = *_threadData[tid];

	// initialize dist lookups
	soa2.initDistLookupPointersSingle(my_threadData._centers_dist_lookup,
				my_threadData._ljc_dist_lookup, soa2._ljc_num);

	// Pointer for molecules
	const vcp_real_calc * const soa1_mol_pos_x = soa1._mol_pos.xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = soa1._mol_pos.yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = soa1._mol_pos.zBegin();

	//for better readability:
	constexpr ConcSites::SiteType LJC = ConcSites::SiteType::LJC;
	typedef ConcSites::CoordinateType	Coordinate;
	typedef CellDataSoA::QuantityType	QuantityType;

	// Pointer for LJ centers
	const vcp_real_calc * const soa1_ljc_r_x = soa1.getBeginCalc(QuantityType::CENTER_POSITION, LJC, Coordinate::X);
	const vcp_real_calc * const soa1_ljc_r_y = soa1.getBeginCalc(QuantityType::CENTER_POSITION, LJC, Coordinate::Y);
	const vcp_real_calc * const soa1_ljc_r_z = soa1.getBeginCalc(QuantityType::CENTER_POSITION, LJC, Coordinate::Z);
	     vcp_real_accum * const soa1_ljc_f_x = soa1.getBeginAccum(QuantityType::FORCE, LJC, Coordinate::X);
	     vcp_real_accum * const soa1_ljc_f_y = soa1.getBeginAccum(QuantityType::FORCE, LJC, Coordinate::Y);
	     vcp_real_accum * const soa1_ljc_f_z = soa1.getBeginAccum(QuantityType::FORCE, LJC, Coordinate::Z);
	     vcp_real_accum * const soa1_ljc_V_x = soa1.getBeginAccum(QuantityType::VIRIAL, LJC, Coordinate::X);
	     vcp_real_accum * const soa1_ljc_V_y = soa1.getBeginAccum(QuantityType::VIRIAL, LJC, Coordinate::Y);
	     vcp_real_accum * const soa1_ljc_V_z = soa1.getBeginAccum(QuantityType::VIRIAL, LJC, Coordinate::Z);
	const int * const soa1_mol_ljc_num = soa1._mol_ljc_num;
	const vcp_ljc_id_t * const soa1_ljc_id = soa1._ljc_id;

	const vcp_real_calc * const soa2_ljc_m_r_x = soa2.getBeginCalc(QuantityType::MOL_POSITION, LJC, Coordinate::X);
	const vcp_real_calc * const soa2_ljc_m_r_y = soa2.getBeginCalc(QuantityType::MOL_POSITION, LJC, Coordinate::Y);
	const vcp_real_calc * const soa2_ljc_m_r_z = soa2.getBeginCalc(QuantityType::MOL_POSITION, LJC, Coordinate::Z);
	const vcp_real_calc * const soa2_ljc_r_x = soa2.getBeginCalc(QuantityType::CENTER_POSITION, LJC, Coordinate::X);
	const vcp_real_calc * const soa2_ljc_r_y = soa2.getBeginCalc(QuantityType::CENTER_POSITION, LJC, Coordinate::Y);
	const vcp_real_calc * const soa2_ljc_r_z = soa2.getBeginCalc(QuantityType::CENTER_POSITION, LJC, Coordinate::Z);
	     vcp_real_accum * const soa2_ljc_f_x = soa2.getBeginAccum(QuantityType::FORCE, LJC, Coordinate::X);
	     vcp_real_accum * const soa2_ljc_f_y = soa2.getBeginAccum(QuantityType::FORCE, LJC, Coordinate::Y);
	     vcp_real_accum * const soa2_ljc_f_z = soa2.getBeginAccum(QuantityType::FORCE, LJC, Coordinate::Z);
	     vcp_real_accum * const soa2_ljc_V_x = soa2.getBeginAccum(QuantityType::VIRIAL, LJC, Coordinate::X);
	     vcp_real_accum * const soa2_ljc_V_y = soa2.getBeginAccum(QuantityType::VIRIAL, LJC, Coordinate::Y);
	     vcp_real_accum * const soa2_ljc_V_z = soa2.getBeginAccum(QuantityType::VIRIAL, LJC, Coordinate::Z);
	const vcp_ljc_id_t * const soa2_ljc_id = soa2._ljc_id;

	vcp_lookupOrMask_single* const soa2_ljc_dist_lookup = my_threadData._ljc_dist_lookup;

	RealAccumVec sum_upot6lj = RealAccumVec::zero();
	RealAccumVec sum_virial = RealAccumVec::zero();

	const RealCalcVec rc2 = RealCalcVec::set1(_LJCutoffRadiusSquare);

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
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);

		size_t end_ljc_loop = MaskGatherChooser::getEndloop(end_ljc_j_longloop, compute_molecule_ljc);

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

						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						RealCalcVec fx, fy, fz;
						RealAccumVec Vx, Vy, Vz;

						RealCalcVec eps_24;
						RealCalcVec sig2;
						unpackEps24Sig2<MaskGatherChooser>(eps_24, sig2, _eps_sig[id_i], soa2_ljc_id, j, lookupORforceMask);

						RealCalcVec shift6;
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
						unpackEps24Sig2<MaskGatherChooser>(eps_24, sig2, _eps_sig[id_i], soa2_ljc_id, j, lookupORforceMask);

						RealCalcVec shift6;
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
	}

	sum_upot6lj.aligned_load_add_store(&my_threadData._upot6ljV[0]);
	sum_virial.aligned_load_add_store(&my_threadData._virialV[0]);
}

void VectorizedLJP2PCellProcessor::processCell(ParticleCell & c) {
	FullParticleCell & full_c = downcastCellReferenceFull(c);
	CellDataSoA& soa = full_c.getCellDataSoA();
	if (full_c.isHaloCell() or soa.getMolNum() < 2) {
		return;
	}
	const bool CalculateMacroscopic = true;
	const bool ApplyCutoff = true;
	_calculatePairs<SingleCellPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa, soa);
}

void VectorizedLJP2PCellProcessor::processCellPair(ParticleCell & c1, ParticleCell & c2, bool sumAll) {
	mardyn_assert(&c1 != &c2);
	FullParticleCell & full_c1 = downcastCellReferenceFull(c1);
	FullParticleCell & full_c2 = downcastCellReferenceFull(c2);

	CellDataSoA& soa1 = full_c1.getCellDataSoA();
	CellDataSoA& soa2 = full_c2.getCellDataSoA();
	const bool c1Halo = full_c1.isHaloCell();
	const bool c2Halo = full_c2.isHaloCell();

	// this variable determines whether
	// _calcPairs(soa1, soa2) or _calcPairs(soa2, soa1)
	// is more efficient
	const bool calc_soa1_soa2 = (soa1.getMolNum() <= soa2.getMolNum());

	
	if(sumAll) { // sumAll
		// if one cell is empty, skip
		if (soa1.getMolNum() == 0 or soa2.getMolNum() == 0) {
			return;
		}

		// Macroscopic conditions: Sum all

		const bool ApplyCutoff = true;

		const bool CalculateMacroscopic = true;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}
	} else { // sumHalf
		// if one cell is empty, or both cells are Halo, skip
		if (soa1.getMolNum() == 0 or soa2.getMolNum() == 0 or (c1Halo and c2Halo)) {
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

} // namespace bhfmm

