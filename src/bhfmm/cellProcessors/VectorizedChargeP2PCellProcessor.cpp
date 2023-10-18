/**
 * \file
 * \brief VectorizedChargeP2PCellProcessor.cpp
 */

#include "particleContainer/adapter/CellDataSoA.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include <algorithm>
#include "particleContainer/adapter/vectorization/MaskGatherChooser.h"
#include "VectorizedChargeP2PCellProcessor.h"

namespace bhfmm {

VectorizedChargeP2PCellProcessor::VectorizedChargeP2PCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		_domain(domain),
		// maybe move the following to somewhere else:
		_upotXpoles(0.0), _virial(0.0){

#if VCP_VEC_TYPE==VCP_NOVEC
	Log::global_log->info() << "VectorizedChargeP2PCellProcessor: using no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	Log::global_log->info() << "VectorizedChargeP2PCellProcessor: using SSE3 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	Log::global_log->info() << "VectorizedChargeP2PCellProcessor: using AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	Log::global_log->info() << "VectorizedChargeP2PCellProcessor: using AVX2 intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_KNL) || (VCP_VEC_TYPE==VCP_VEC_KNL_GATHER)
	Log::global_log->info() << "VectorizedChargeP2PCellProcessor: using KNL intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_AVX512F) || (VCP_VEC_TYPE==VCP_VEC_AVX512F_GATHER)
	Log::global_log->info() << "VectorizedChargeP2PCellProcessor: using SKX intrinsics." << std::endl;
#endif

	// initialize thread data
	_numThreads = mardyn_get_max_threads();
	Log::global_log->info() << "VectorizedChargeP2PCellProcessor: allocate data for " << _numThreads << " threads." << std::endl;
	_threadData.resize(_numThreads);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		VCP2PCPThreadData * myown = new VCP2PCPThreadData();
		const int myid = mardyn_get_thread_num();
		_threadData[myid] = myown;
	} // end pragma omp parallel

#ifdef ENABLE_MPI
	global_simulation->timers()->setOutputString("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P", "FMM: Time spent in Charge P2P ");
	//global_simulation->setSyncTimer("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P", false); //it is per default false
#endif
}

VectorizedChargeP2PCellProcessor :: ~VectorizedChargeP2PCellProcessor () {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int myid = mardyn_get_thread_num();
		delete _threadData[myid];
	}
}

void VectorizedChargeP2PCellProcessor::printTimers() {
	std::cout << "FMM: Time spent in Charge P2P " << global_simulation->timers()->getTime("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P") << std::endl;
	global_simulation->timers()->print("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P");
}

void VectorizedChargeP2PCellProcessor::initTraversal() {
	global_simulation->timers()->start("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P");

	#if defined(_OPENMP)
	#pragma omp master
	#endif
	{
		_upotXpoles = 0.0;
		_virial = 0.0;
	} // end pragma omp master

}

void VectorizedChargeP2PCellProcessor::endTraversal() {
	double currentVirial = _domain.getLocalVirial();
	double currentUpot = _domain.getLocalUpot();
	double glob_upotXpoles = 0.0;
	double glob_virial = 0.0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:glob_upotXpoles, glob_virial)
	#endif
	{
		const int tid = mardyn_get_thread_num();

		// reduce vectors and clear local variable
		vcp_real_accum thread_upotXpoles = 0.0, thread_virial = 0.0;

		load_hSum_Store_Clear(&thread_upotXpoles, _threadData[tid]->_upotXpolesV);
		load_hSum_Store_Clear(&thread_virial, _threadData[tid]->_virialV);

		// add to global sum
		glob_upotXpoles += thread_upotXpoles;
		glob_virial += thread_virial;
	} // end pragma omp parallel reduction

	_upotXpoles = glob_upotXpoles;
	_virial = glob_virial;
	_domain.setLocalVirial(currentVirial + _virial);
	_domain.setLocalUpot(currentUpot + _upotXpoles);
	global_simulation->timers()->stop("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P");
}

void VectorizedChargeP2PCellProcessor::preprocessCell(ParticleCellPointers & c) {
	// as pre new integration of Caches in SoAs,
	// this function work as before, as it builds secondary SoAs

	// Determine the total number of centers.
	size_t numMolecules = c.getMoleculeCount();
	if(numMolecules == 0)
		return;
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;

	for (size_t i = 0; i < numMolecules; ++i) {
		nCharges += c.moleculesAt(i).numCharges();
	}

	// Construct the SoA.
	CellDataSoA & soa = c.getCellDataSoA();
	soa.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);

	ComponentList components = *(_simulation.getEnsemble()->getComponents());

	//for better readability:
	constexpr ConcSites::SiteType CHARGE = ConcSites::SiteType::CHARGE;
	typedef ConcSites::CoordinateType Coordinate;
	typedef CellDataSoA::QuantityType QuantityType;

	vcp_real_calc* const soa_charges_m_r_x = soa.getBeginCalc(QuantityType::MOL_POSITION, CHARGE, Coordinate::X);
	vcp_real_calc* const soa_charges_m_r_y = soa.getBeginCalc(QuantityType::MOL_POSITION, CHARGE, Coordinate::Y);
	vcp_real_calc* const soa_charges_m_r_z = soa.getBeginCalc(QuantityType::MOL_POSITION, CHARGE, Coordinate::Z);
	vcp_real_calc* const soa_charges_r_x = soa.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::X);
	vcp_real_calc* const soa_charges_r_y = soa.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::Y);
	vcp_real_calc* const soa_charges_r_z = soa.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::Z);
	vcp_real_accum* const soa_charges_f_x = soa.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_f_y = soa.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_f_z = soa.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Z);
	vcp_real_accum* const soa_charges_V_xx = soa.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_V_yy = soa.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_V_zz = soa.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Z);
	vcp_real_accum* const soa_charges_V_xy = soa.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_V_xz = soa.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_V_yz = soa.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Z);
	vcp_real_accum* const soa_charges_V_yx = soa.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_V_zx = soa.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_V_zy = soa.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Z);

	size_t iCharges = 0;
	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < numMolecules; ++i) {
		Molecule& mol = c.moleculesAt(i);

		const size_t mol_charges_num = mol.numCharges();
		const vcp_real_calc mol_pos_x = mol.r(0);
		const vcp_real_calc mol_pos_y = mol.r(1);
		const vcp_real_calc mol_pos_z = mol.r(2);

		soa._mol_pos.x(i) = mol_pos_x;
		soa._mol_pos.y(i) = mol_pos_y;
		soa._mol_pos.z(i) = mol_pos_z;
		soa._mol_charges_num[i] = mol_charges_num;

		for (size_t j = 0; j < mol_charges_num; ++j, ++iCharges)
		{
			soa_charges_m_r_x[iCharges] = mol_pos_x;
			soa_charges_m_r_y[iCharges] = mol_pos_y;
			soa_charges_m_r_z[iCharges] = mol_pos_z;
			soa_charges_r_x[iCharges] = mol.charge_d(j)[0] + mol_pos_x;
			soa_charges_r_y[iCharges] = mol.charge_d(j)[1] + mol_pos_y;
			soa_charges_r_z[iCharges] = mol.charge_d(j)[2] + mol_pos_z;
			soa_charges_f_x[iCharges] = 0.0;
			soa_charges_f_y[iCharges] = 0.0;
			soa_charges_f_z[iCharges] = 0.0;
			soa_charges_V_xx[iCharges] = 0.0;
			soa_charges_V_yy[iCharges] = 0.0;
			soa_charges_V_zz[iCharges] = 0.0;
			soa_charges_V_xy[iCharges] = 0.0;
			soa_charges_V_xz[iCharges] = 0.0;
			soa_charges_V_yz[iCharges] = 0.0;
			soa_charges_V_yx[iCharges] = 0.0;
			soa_charges_V_zx[iCharges] = 0.0;
			soa_charges_V_zy[iCharges] = 0.0;
			//soa._charges_dist_lookup[iCharges] = 0.0;
			// Get the charge
			soa._charges_q[iCharges] = components[mol.componentid()].charge(j).q();
		}
	}
}

void VectorizedChargeP2PCellProcessor::postprocessCell(ParticleCellPointers & c) {
	// as pre new integration of Caches in SoAs,
	// this function work as before, as it builds secondary SoAs
	using std::isnan; // C++11 required

	size_t numMolecules = c.getMoleculeCount();
	if(numMolecules == 0)
		return;

	CellDataSoA& soa = c.getCellDataSoA();

	//for better readability:
	constexpr ConcSites::SiteType CHARGE = ConcSites::SiteType::CHARGE;
	typedef ConcSites::CoordinateType Coordinate;
	typedef CellDataSoA::QuantityType QuantityType;

	vcp_real_accum* const soa_charges_f_x = soa.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_f_y = soa.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_f_z = soa.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Z);
	vcp_real_accum* const soa_charges_V_xx = soa.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_V_yy = soa.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_V_zz = soa.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Z);
	vcp_real_accum* const soa_charges_V_xy = soa.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_V_xz = soa.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_V_yz = soa.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Z);
	vcp_real_accum* const soa_charges_V_yx = soa.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::X);
	vcp_real_accum* const soa_charges_V_zx = soa.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Y);
	vcp_real_accum* const soa_charges_V_zy = soa.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Z);

	// For each molecule iterate over all its centers.
	size_t iCharges = 0;
	for (size_t i = 0; i < numMolecules; ++i) {
		Molecule& m = c.moleculesAt(i);

		const size_t mol_charges_num = m.numCharges();

		for (size_t j = 0; j < mol_charges_num; ++j, ++iCharges) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = static_cast<double>(soa_charges_f_x[iCharges]);
			f[1] = static_cast<double>(soa_charges_f_y[iCharges]);
			f[2] = static_cast<double>(soa_charges_f_z[iCharges]);
			mardyn_assert(!std::isnan(f[0]));
			mardyn_assert(!std::isnan(f[1]));
			mardyn_assert(!std::isnan(f[2]));
			m.Fchargeadd(j, f);

			// Store the resulting virial in the molecule.
			double V[9];
			V[0] = static_cast<double>(soa_charges_V_xx[iCharges]*0.5);
			V[1] = static_cast<double>(soa_charges_V_yy[iCharges]*0.5);
			V[2] = static_cast<double>(soa_charges_V_zz[iCharges]*0.5);
			V[3] = static_cast<double>(soa_charges_V_xy[iCharges]*0.5);
			V[4] = static_cast<double>(soa_charges_V_xz[iCharges]*0.5);
			V[5] = static_cast<double>(soa_charges_V_yz[iCharges]*0.5);
			V[6] = static_cast<double>(soa_charges_V_yx[iCharges]*0.5);
			V[7] = static_cast<double>(soa_charges_V_zx[iCharges]*0.5);
			V[8] = static_cast<double>(soa_charges_V_zy[iCharges]*0.5);
			mardyn_assert(!std::isnan(V[0]));
			mardyn_assert(!std::isnan(V[1]));
			mardyn_assert(!std::isnan(V[2]));
			mardyn_assert(!std::isnan(V[3]));
			mardyn_assert(!std::isnan(V[4]));
			mardyn_assert(!std::isnan(V[5]));
			mardyn_assert(!std::isnan(V[6]));
			mardyn_assert(!std::isnan(V[7]));
			mardyn_assert(!std::isnan(V[8]));
			m.Viadd(V);
		}
	}
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
inline void VectorizedChargeP2PCellProcessor :: _loopBodyCharge(
		const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
		const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
		const RealCalcVec& qii,
		const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
		const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
		const RealCalcVec& qjj,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealAccumVec& V_xx, RealAccumVec& V_yy, RealAccumVec& V_zz,
		RealAccumVec& V_xy, RealAccumVec& V_xz, RealAccumVec& V_yz,
		RealAccumVec& V_yx, RealAccumVec& V_zx, RealAccumVec& V_zy,
		RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
		const MaskCalcVec& forceMask)
{
	const RealCalcVec c_dx = r1_x - r2_x;
	const RealCalcVec c_dy = r1_y - r2_y;
	const RealCalcVec c_dz = r1_z - r2_z;//fma not possible since they will be reused...

	const RealCalcVec c_dr2 = RealCalcVec::scal_prod(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

	const RealCalcVec c_dr2_inv_unmasked = one / c_dr2;
	const RealCalcVec c_dr2_inv = RealCalcVec::apply_mask(c_dr2_inv_unmasked, forceMask);//masked
	const RealCalcVec c_dr_inv = RealCalcVec::sqrt(c_dr2_inv);//masked

	const RealCalcVec q1q2per4pie0 = qii * qjj;
	const RealCalcVec upot = q1q2per4pie0 * c_dr_inv;//masked
	const RealCalcVec fac = upot * c_dr2_inv;//masked

	f_x = c_dx * fac;
	f_y = c_dy * fac;
	f_z = c_dz * fac;
	const RealCalcVec m_dx = m1_r_x - m2_r_x;
	const RealCalcVec m_dy = m1_r_y - m2_r_y;
	const RealCalcVec m_dz = m1_r_z - m2_r_z;

	V_xx = RealAccumVec::convertCalcToAccum(m_dx * f_x);
	V_yy = RealAccumVec::convertCalcToAccum(m_dy * f_y);
	V_zz = RealAccumVec::convertCalcToAccum(m_dz * f_z);
	V_xy = RealAccumVec::convertCalcToAccum(m_dx * f_y);
	V_xz = RealAccumVec::convertCalcToAccum(m_dx * f_z);
	V_yz = RealAccumVec::convertCalcToAccum(m_dy * f_z);
	V_yx = RealAccumVec::convertCalcToAccum(m_dy * f_x);
	V_zx = RealAccumVec::convertCalcToAccum(m_dz * f_x);
	V_zy = RealAccumVec::convertCalcToAccum(m_dz * f_y);
	// Check if we have to add the macroscopic values up
	if (calculateMacroscopic) {
		const RealAccumVec upot_accum = RealAccumVec::convertCalcToAccum(upot);
		sum_upotXpoles = sum_upotXpoles + upot_accum;
		sum_virial = sum_virial + V_xx + V_yy + V_zz;//DoubleVec::scal_prod(m_dx, m_dy, m_dz, f_x, f_y, f_z);
	}
}

template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
void VectorizedChargeP2PCellProcessor::_calculatePairs(CellDataSoA & soa1, CellDataSoA & soa2) {
	const int tid = mardyn_get_thread_num();
	VCP2PCPThreadData &my_threadData = *_threadData[tid];

	// initialize dist lookups
	soa2.initDistLookupPointersSingle(my_threadData._centers_dist_lookup,
			my_threadData._charges_dist_lookup, soa2._charges_num);

	// Pointer for molecules
	const vcp_real_calc * const soa1_mol_pos_x = soa1._mol_pos.xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = soa1._mol_pos.yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = soa1._mol_pos.zBegin();

	//for better readability:
	constexpr ConcSites::SiteType CHARGE = ConcSites::SiteType::CHARGE;
	typedef ConcSites::CoordinateType	Coordinate;
	typedef CellDataSoA::QuantityType							QuantityType;

	// Pointer for charges
	const vcp_real_calc * const soa1_charges_r_x = soa1.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::X);
	const vcp_real_calc * const soa1_charges_r_y = soa1.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::Y);
	const vcp_real_calc * const soa1_charges_r_z = soa1.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa1_charges_f_x = soa1.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa1_charges_f_y = soa1.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa1_charges_f_z = soa1.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa1_charges_V_xx = soa1.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa1_charges_V_yy = soa1.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa1_charges_V_zz = soa1.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa1_charges_V_xy = soa1.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa1_charges_V_xz = soa1.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa1_charges_V_yz = soa1.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa1_charges_V_yx = soa1.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa1_charges_V_zx = soa1.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa1_charges_V_zy = soa1.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Z);
	const vcp_real_calc * const soa1_charges_q = soa1._charges_q;
	const int * const soa1_mol_charges_num = soa1._mol_charges_num;

	const vcp_real_calc * const soa2_charges_m_r_x = soa2.getBeginCalc(QuantityType::MOL_POSITION, CHARGE, Coordinate::X);
	const vcp_real_calc * const soa2_charges_m_r_y = soa2.getBeginCalc(QuantityType::MOL_POSITION, CHARGE, Coordinate::Y);
	const vcp_real_calc * const soa2_charges_m_r_z = soa2.getBeginCalc(QuantityType::MOL_POSITION, CHARGE, Coordinate::Z);
	const vcp_real_calc * const soa2_charges_r_x = soa2.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::X);
	const vcp_real_calc * const soa2_charges_r_y = soa2.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::Y);
	const vcp_real_calc * const soa2_charges_r_z = soa2.getBeginCalc(QuantityType::CENTER_POSITION, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa2_charges_f_x = soa2.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa2_charges_f_y = soa2.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa2_charges_f_z = soa2.getBeginAccum(QuantityType::FORCE, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa2_charges_V_xx = soa2.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa2_charges_V_yy = soa2.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa2_charges_V_zz = soa2.getBeginAccum(QuantityType::VIRIAL, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa2_charges_V_xy = soa2.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa2_charges_V_xz = soa2.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa2_charges_V_yz = soa2.getBeginAccum(QuantityType::VIRIALND1, CHARGE, Coordinate::Z);
		 vcp_real_accum * const soa2_charges_V_yx = soa2.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::X);
		 vcp_real_accum * const soa2_charges_V_zx = soa2.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Y);
		 vcp_real_accum * const soa2_charges_V_zy = soa2.getBeginAccum(QuantityType::VIRIALND2, CHARGE, Coordinate::Z);
	const vcp_real_calc * const soa2_charges_q = soa2._charges_q;

	vcp_lookupOrMask_single* const soa2_charges_dist_lookup = my_threadData._charges_dist_lookup;


	RealAccumVec sum_upotXpoles = RealAccumVec::zero();
	RealAccumVec sum_virial = RealAccumVec::zero();

	const RealCalcVec cutoffRadiusSquare = RealCalcVec::set1(_cutoffRadiusSquare);

	/*
	 *  Here different end values for the loops are defined. For loops, which do not vectorize over the last (possibly "uneven") amount of indices, the normal values are computed. These mark the end of the vectorized part.
	 *  The longloop values mark the end of the vectorized part, if vectorization is performed for all elements. For these, various conditions have to be fulfilled, to be sure, that no NaN values are stored and no segfaults arise:
	 *  * arrays have to be long enough
	 *  * arrays have to be filled with something
	 *  * _ljc_id has to be set to existing values for each index
	 *  All of these conditions are fulfilled by setting the non existing values within CellDataSoA.h and AlignedArray.h to zero.
	 */
	const size_t end_charges_j = vcp_floor_to_vec_size(soa2._charges_num);
	const size_t end_charges_j_longloop = vcp_ceil_to_vec_size(soa2._charges_num);//this is ceil _charges_num, VCP_VEC_SIZE
	size_t i_charge_idx = 0;

	// Iterate over each center in the first cell.
	const size_t soa1_mol_num = soa1.getMolNum();
	for (size_t i = 0; i < soa1_mol_num; ++i) {//over the molecules
		const RealCalcVec m1_r_x = RealCalcVec::broadcast(soa1_mol_pos_x + i);
		const RealCalcVec m1_r_y = RealCalcVec::broadcast(soa1_mol_pos_y + i);
		const RealCalcVec m1_r_z = RealCalcVec::broadcast(soa1_mol_pos_z + i);
		// Iterate over centers of second cell
		const countertype32 compute_molecule_charges = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_charge_idx, soa2._charges_num,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);

		size_t end_charges_loop = MaskGatherChooser::getEndloop(end_charges_j_longloop, compute_molecule_charges);

		// Computation of site interactions with charges

		if (compute_molecule_charges==0) {
			i_charge_idx += soa1_mol_charges_num[i];
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

				RealAccumVec sum_V1_xx = RealAccumVec::zero();
				RealAccumVec sum_V1_yy = RealAccumVec::zero();
				RealAccumVec sum_V1_zz = RealAccumVec::zero();
				RealAccumVec sum_V1_xy = RealAccumVec::zero();
				RealAccumVec sum_V1_xz = RealAccumVec::zero();
				RealAccumVec sum_V1_yz = RealAccumVec::zero();
				RealAccumVec sum_V1_yx = RealAccumVec::zero();
				RealAccumVec sum_V1_zx = RealAccumVec::zero();
				RealAccumVec sum_V1_zy = RealAccumVec::zero();

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
						RealAccumVec Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, Vyx, Vzx, Vzy;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vxx, Vyy, Vzz,
								Vxy, Vxz, Vyz,
								Vyx, Vzx, Vzy,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));

						RealAccumVec a_fx = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_fy = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_fz = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_fx;
						sum_f1_y = sum_f1_y + a_fy;
						sum_f1_z = sum_f1_z + a_fz;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_x, j, a_fx, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_y, j, a_fy, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_z, j, a_fz, lookupORforceMask);

						sum_V1_xx = sum_V1_xx + Vxx;
						sum_V1_yy = sum_V1_yy + Vyy;
						sum_V1_zz = sum_V1_zz + Vzz;
						sum_V1_xy = sum_V1_xy + Vxy;
						sum_V1_xz = sum_V1_xz + Vxz;
						sum_V1_yz = sum_V1_yz + Vyz;
						sum_V1_yx = sum_V1_yx + Vyx;
						sum_V1_zx = sum_V1_zx + Vzx;
						sum_V1_zy = sum_V1_zy + Vzy;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_xx, j, Vxx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_yy, j, Vyy, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_zz, j, Vzz, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_xy, j, Vxy, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_xz, j, Vxz, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_yz, j, Vyz, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_yx, j, Vyx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_zx, j, Vzx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_zy, j, Vzy, lookupORforceMask);
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER  or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
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
						RealAccumVec Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, Vyx, Vzx, Vzy;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vxx, Vyy, Vzz,
								Vxy, Vxz, Vyz,
								Vyx, Vzx, Vzy,
								sum_upotXpoles, sum_virial,
								remainderM);

						RealAccumVec a_fx = RealAccumVec::convertCalcToAccum(f_x);
						RealAccumVec a_fy = RealAccumVec::convertCalcToAccum(f_y);
						RealAccumVec a_fz = RealAccumVec::convertCalcToAccum(f_z);

						sum_f1_x = sum_f1_x + a_fx;
						sum_f1_y = sum_f1_y + a_fy;
						sum_f1_z = sum_f1_z + a_fz;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, a_fx, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, a_fy, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, a_fz, lookupORforceMask, remainderM);

						sum_V1_xx = sum_V1_xx + Vxx;
						sum_V1_yy = sum_V1_yy + Vyy;
						sum_V1_zz = sum_V1_zz + Vzz;
						sum_V1_xy = sum_V1_xy + Vxy;
						sum_V1_xz = sum_V1_xz + Vxz;
						sum_V1_yz = sum_V1_yz + Vyz;
						sum_V1_yx = sum_V1_yx + Vyx;
						sum_V1_zx = sum_V1_zx + Vzx;
						sum_V1_zy = sum_V1_zy + Vzy;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_xx, j, Vxx, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_yy, j, Vyy, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_zz, j, Vzz, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_xy, j, Vxy, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_xz, j, Vxz, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_yz, j, Vyz, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_yx, j, Vyx, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_zx, j, Vzx, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_zy, j, Vzy, lookupORforceMask, remainderM);
					}
				}
#endif

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx + local_i, sum_f1_z);
				// Add old virial and summed calculated virials for center 1
				hSum_Add_Store(soa1_charges_V_xx + i_charge_idx + local_i, sum_V1_xx);
				hSum_Add_Store(soa1_charges_V_yy + i_charge_idx + local_i, sum_V1_yy);
				hSum_Add_Store(soa1_charges_V_zz + i_charge_idx + local_i, sum_V1_zz);
				hSum_Add_Store(soa1_charges_V_xy + i_charge_idx + local_i, sum_V1_xy);
				hSum_Add_Store(soa1_charges_V_xz + i_charge_idx + local_i, sum_V1_xz);
				hSum_Add_Store(soa1_charges_V_yz + i_charge_idx + local_i, sum_V1_yz);
				hSum_Add_Store(soa1_charges_V_yx + i_charge_idx + local_i, sum_V1_yx);
				hSum_Add_Store(soa1_charges_V_zx + i_charge_idx + local_i, sum_V1_zx);
				hSum_Add_Store(soa1_charges_V_zy + i_charge_idx + local_i, sum_V1_zy);
			}
			i_charge_idx += soa1_mol_charges_num[i];
		}
	}

	sum_upotXpoles.aligned_load_add_store(&my_threadData._upotXpolesV[0]);
	sum_virial.aligned_load_add_store(&my_threadData._virialV[0]);
} // void VectorizedChargeP2PCellProcessor::_calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2)

void VectorizedChargeP2PCellProcessor::processCell(ParticleCellPointers & c) {
	CellDataSoA& soa = c.getCellDataSoA();
	if (c.isHaloCell() or soa.getMolNum() < 2) {
		return;
	}
	const bool CalculateMacroscopic = true;
	const bool ApplyCutoff = false;
	_calculatePairs<SingleCellPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa, soa);
}

void VectorizedChargeP2PCellProcessor::processCellPair(ParticleCellPointers & c1, ParticleCellPointers & c2) {
	mardyn_assert(&c1 != &c2);
	CellDataSoA& soa1 = c1.getCellDataSoA();
	CellDataSoA& soa2 = c2.getCellDataSoA();
	const bool c1Halo = c1.isHaloCell();
	const bool c2Halo = c2.isHaloCell();

	// this variable determines whether
	// _calcPairs(soa1, soa2) or _calcPairs(soa2, soa1)
	// is more efficient
	const bool calc_soa1_soa2 = (soa1.getMolNum() <= soa2.getMolNum());

	// if one cell is empty, or both cells are Halo, skip
	if (soa1.getMolNum() == 0 or soa2.getMolNum() == 0 or (c1Halo and c2Halo)) {
		return;
	}

	// Macroscopic conditions:
	// if none of the cells is halo, then compute
	// if one of them is halo:
	// 		if c1-index < c2-index, then compute
	// 		else, then don't compute
	// This saves the Molecule::isLessThan checks
	// and works similar to the "Half-Shell" scheme

	const bool ApplyCutoff = false;

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
		mardyn_assert(c1Halo != c2Halo);							// one of them is halo and
		mardyn_assert(not (c1.getCellIndex() < c2.getCellIndex()));// c1.index not < c2.index

		const bool CalculateMacroscopic = false;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}
	}
}

} // namespace bhfmm

