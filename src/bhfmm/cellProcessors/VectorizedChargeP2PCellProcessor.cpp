/**
 * \file
 * \brief VectorizedChargeP2PCellProcessor.cpp
 */

#include "particleContainer/adapter/CellDataSoA.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include <algorithm>
#include "particleContainer/adapter/vectorization/MaskGatherChooser.h"
#include "VectorizedChargeP2PCellProcessor.h"

using namespace Log;
namespace bhfmm {

VectorizedChargeP2PCellProcessor::VectorizedChargeP2PCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		CellProcessor(cutoffRadius, LJcutoffRadius), _domain(domain),
		// maybe move the following to somewhere else:
		_compIDs(), _upotXpoles(0.0), _virial(0.0), _charges_dist_lookup(0) {

#if VCP_VEC_TYPE==VCP_NOVEC
	global_log->info() << "VectorizedChargeP2PCellProcessor: using no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	global_log->info() << "VectorizedChargeP2PCellProcessor: using SSE3 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	global_log->info() << "VectorizedChargeP2PCellProcessor: using AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	global_log->info() << "VectorizedChargeP2PCellProcessor: using AVX2 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_MIC
	global_log->info() << "VectorizedChargeP2PCellProcessor: using MIC intrinsics." << std::endl;
#endif

	ComponentList components = *(_simulation.getEnsemble()->components());
	// Get the maximum Component ID.
	size_t maxID = 0;
	const ComponentList::const_iterator end = components.end();
	for (ComponentList::const_iterator c = components.begin(); c != end; ++c)
		maxID = std::max(maxID, static_cast<size_t>(c->ID()));

	// Assign a center list start index for each component.
	_compIDs.resize(maxID + 1, 0);
	size_t centers = 0;
	for (ComponentList::const_iterator c = components.begin(); c != end; ++c) {
		_compIDs[c->ID()] = centers;
		centers += c->numLJcenters();
	}

#ifdef ENABLE_MPI
	_timer.set_sync(false);
#endif
}

VectorizedChargeP2PCellProcessor :: ~VectorizedChargeP2PCellProcessor () {
}

void VectorizedChargeP2PCellProcessor::printTimers() {
	std::cout << "FMM: Time spent in Charge P2P " << _timer.get_etime() << std::endl;
}


void VectorizedChargeP2PCellProcessor::initTraversal(const size_t numCells) {
	_timer.start();
	_virial = 0.0;
	_upotXpoles = 0.0;
}


void VectorizedChargeP2PCellProcessor::endTraversal() {
	double currentVirial = _domain.getLocalVirial();
	double currentUpot = _domain.getLocalUpot();

	_domain.setLocalVirial(currentVirial + _virial);
	_domain.setLocalUpot(currentUpot + _upotXpoles);
	_timer.stop();
}


void VectorizedChargeP2PCellProcessor::preprocessCell(ParticleCell & c) {
	const MoleculeList & molecules = c.getParticlePointers();

	// Determine the total number of centers.
	size_t numMolecules = molecules.size();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;
	
	for (size_t m = 0;  m < numMolecules; ++m) {
		nCharges += molecules[m]->numCharges();
	}

	// Construct the SoA.
	CellDataSoA & soa = c.getCellDataSoA();
	soa.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);

	ComponentList components = *(_simulation.getEnsemble()->components());

	double* const soa_charges_m_r_x = soa.charges_m_r_xBegin();
	double* const soa_charges_m_r_y = soa.charges_m_r_yBegin();
	double* const soa_charges_m_r_z = soa.charges_m_r_zBegin();
	double* const soa_charges_r_x = soa.charges_r_xBegin();
	double* const soa_charges_r_y = soa.charges_r_yBegin();
	double* const soa_charges_r_z = soa.charges_r_zBegin();
	double* const soa_charges_f_x = soa.charges_f_xBegin();
	double* const soa_charges_f_y = soa.charges_f_yBegin();
	double* const soa_charges_f_z = soa.charges_f_zBegin();
	double* const soa_charges_V_x = soa.charges_V_xBegin();
	double* const soa_charges_V_y = soa.charges_V_yBegin();
	double* const soa_charges_V_z = soa.charges_V_zBegin();

	size_t iCharges = 0;
	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < molecules.size(); ++i) {
		const size_t mol_charges_num = molecules[i]->numCharges();
		const double mol_pos_x = molecules[i]->r(0);
		const double mol_pos_y = molecules[i]->r(1);
		const double mol_pos_z = molecules[i]->r(2);

		soa._mol_pos.x(i) = mol_pos_x;
		soa._mol_pos.y(i) = mol_pos_y;
		soa._mol_pos.z(i) = mol_pos_z;
		soa._mol_charges_num[i] = mol_charges_num;

		for (size_t j = 0; j < mol_charges_num; ++j, ++iCharges)
		{
			soa_charges_m_r_x[iCharges] = mol_pos_x;
			soa_charges_m_r_y[iCharges] = mol_pos_y;
			soa_charges_m_r_z[iCharges] = mol_pos_z;
			soa_charges_r_x[iCharges] = molecules[i]->charge_d(j)[0] + mol_pos_x;
			soa_charges_r_y[iCharges] = molecules[i]->charge_d(j)[1] + mol_pos_y;
			soa_charges_r_z[iCharges] = molecules[i]->charge_d(j)[2] + mol_pos_z;
			soa_charges_f_x[iCharges] = 0.0;
			soa_charges_f_y[iCharges] = 0.0;
			soa_charges_f_z[iCharges] = 0.0;
			soa_charges_V_x[iCharges] = 0.0;
			soa_charges_V_y[iCharges] = 0.0;
			soa_charges_V_z[iCharges] = 0.0;
			//soa._charges_dist_lookup[iCharges] = 0.0;
			// Get the charge
			soa._charges_q[iCharges] = components[molecules[i]->componentid()].charge(j).q();
		}
	}
}


void VectorizedChargeP2PCellProcessor::postprocessCell(ParticleCell & c) {
	CellDataSoA& soa = c.getCellDataSoA();

	MoleculeList & molecules = c.getParticlePointers();

	double* const soa_charges_f_x = soa.charges_f_xBegin();
	double* const soa_charges_f_y = soa.charges_f_yBegin();
	double* const soa_charges_f_z = soa.charges_f_zBegin();
	double* const soa_charges_V_x = soa.charges_V_xBegin();
	double* const soa_charges_V_y = soa.charges_V_yBegin();
	double* const soa_charges_V_z = soa.charges_V_zBegin();

	// For each molecule iterate over all its centers.
	size_t iCharges = 0;
	size_t numMols = molecules.size();
	for (size_t m = 0; m < numMols; ++m) {
		const size_t mol_ljc_num = molecules[m]->numLJcenters();
		const size_t mol_charges_num = molecules[m]->numCharges();
		const size_t mol_dipoles_num = molecules[m]->numDipoles();
		const size_t mol_quadrupoles_num = molecules[m]->numQuadrupoles();

		for (size_t i = 0; i < mol_charges_num; ++i, ++iCharges) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa_charges_f_x[iCharges];
			f[1] = soa_charges_f_y[iCharges];
			f[2] = soa_charges_f_z[iCharges];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			molecules[m]->Fchargeadd(i, f);

			// Store the resulting virial in the molecule.
			double V[3];
			V[0] = soa_charges_V_x[iCharges]*0.5;
			V[1] = soa_charges_V_y[iCharges]*0.5;
			V[2] = soa_charges_V_z[iCharges]*0.5;
			assert(!isnan(V[0]));
			assert(!isnan(V[1]));
			assert(!isnan(V[2]));
			molecules[m]->Viadd(V);
		}
	}
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
	inline void VectorizedChargeP2PCellProcessor :: _loopBodyCharge(
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


template<class ForcePolicy>
vcp_mask_vec
inline VectorizedChargeP2PCellProcessor::calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		vcp_lookupOrMask_single* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z,
		const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z,
		countertype32& counter
		) {

#if VCP_VEC_TYPE==VCP_NOVEC

	bool compute_molecule = false;

	for (size_t j = ForcePolicy :: InitJ(i_center_idx); j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos.x(i) - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos.y(i) - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos.z(i) - soa2_m_r_z[j];
		const double m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const bool forceMask = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? true : false;
		*(soa2_center_dist_lookup + j) = forceMask;
		compute_molecule |= forceMask;

	}

	return compute_molecule;

#elif VCP_VEC_TYPE==VCP_VEC_SSE3

	vcp_mask_vec compute_molecule = VCP_SIMD_ZEROVM;

	// Iterate over centers of second cell
	size_t j = ForcePolicy :: InitJ(i_center_idx);
	for (; j < end_j; j+=VCP_VEC_SIZE) {//end_j is chosen accordingly when function is called. (teilbar durch VCP_VEC_SIZE)
		const vcp_double_vec m2_r_x = vcp_simd_load(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = vcp_simd_load(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = vcp_simd_load(soa2_m_r_z + j);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const vcp_mask_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD);
		vcp_simd_store(soa2_center_dist_lookup + j, forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask);
	}

	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos.x(i) - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos.y(i) - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos.z(i) - soa2_m_r_z[j];

		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		const unsigned long forceMask = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;

		*(soa2_center_dist_lookup + j) = forceMask;
		const vcp_mask_vec forceMask_vec = vcp_simd_set1(forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask_vec);
	}
	memset(soa2_center_dist_lookup + j, 0, ((VCP_VEC_SIZE-(soa2_num_centers-end_j)) % VCP_VEC_SIZE) * sizeof(vcp_mask_single));//set the remaining values to zero.
	//This is needed to allow vectorization even of the last elements, their count does not necessarily divide by VCP_VEC_SIZE.
	//The array size is however long enough to vectorize over the last few entries. This sets the entries, that do not make sense in that vectorization to zero.

	return compute_molecule;

#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2

	vcp_mask_vec compute_molecule = VCP_SIMD_ZEROVM;

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	vcp_mask_vec initJ_mask = ForcePolicy::InitJ_Mask(i_center_idx);
	// Iterate over centers of second cell
	for (; j < end_j; j+=VCP_VEC_SIZE) {//end_j is chosen accordingly when function is called. (teilbar durch VCP_VEC_SIZE)
		const vcp_double_vec m2_r_x = vcp_simd_load(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = vcp_simd_load(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = vcp_simd_load(soa2_m_r_z + j);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const vcp_mask_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);
		vcp_simd_store(soa2_center_dist_lookup + j, forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask);
	}

	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos.x(i) - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos.y(i) - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos.z(i) - soa2_m_r_z[j];
		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		// can we do this nicer?
		vcp_mask_single forceMask;
			// DetectSingleCell() = false for SingleCellDistinctPolicy and CellPairPolicy, true for SingleCellPolicy
			// we need this, since in contrast to sse3 we can no longer guarantee, that j>=i by default (j==i is handled by ForcePolicy::Condition).
			// however only one of the branches should be chosen by the compiler, since the class is known at compile time.
		if (ForcePolicy::DetectSingleCell()) {
			forceMask = (ForcePolicy::Condition(m_r2, cutoffRadiusSquare) && j > i_center_idx) ? ~0l : 0l;
		} else {
			forceMask = ForcePolicy::Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;
		}

		*(soa2_center_dist_lookup + j) = forceMask;
		const vcp_mask_vec forceMask_vec = vcp_simd_set1(forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask_vec);
	}
	memset(soa2_center_dist_lookup + j, 0, ((VCP_VEC_SIZE-(soa2_num_centers-end_j)) % VCP_VEC_SIZE) * sizeof(vcp_mask_single));//set the remaining values to zero.
	//This is needed to allow vectorization even of the last elements, their count does not necessarily divide by VCP_VEC_SIZE.
	//The array size is however long enough to vectorize over the last few entries. This sets the entries, that do not make sense in that vectorization to zero.

	return compute_molecule;
#elif VCP_VEC_TYPE==VCP_VEC_MIC
	vcp_mask_vec compute_molecule = VCP_SIMD_ZEROVM;

	size_t j = ForcePolicy :: InitJ(i_center_idx);//j=0 or multiple of vec_size
	vcp_mask_vec initJ_mask = ForcePolicy::InitJ_Mask(i_center_idx);
	// Iterate over centers of second cell
	for (; j < end_j; j+=VCP_VEC_SIZE) {//end_j is chosen accordingly when function is called. (teilbar durch VCP_VEC_SIZE)
		const vcp_double_vec m2_r_x = vcp_simd_load(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = vcp_simd_load(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = vcp_simd_load(soa2_m_r_z + j);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const vcp_mask_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);
		vcp_simd_store(soa2_center_dist_lookup + j/VCP_INDICES_PER_LOOKUP_SINGLE, forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask);
	}//end_j is floor_to_vec_size
	//j is now floor_to_vec_size(soa2_num_centers)=end_j

	//last indices are more complicated for MIC, since masks look different
	size_t k = (end_j+VCP_INDICES_PER_LOOKUP_SINGLE_M1)/VCP_INDICES_PER_LOOKUP_SINGLE;
	vcp_mask_vec forceMask = VCP_SIMD_ZEROVM;
	unsigned char bitmultiplier = 1;
	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j, bitmultiplier *= 2) {
		const double m_dx = soa1._mol_pos.x(i) - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos.y(i) - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos.z(i) - soa2_m_r_z[j];

		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		// can we do this nicer?
		unsigned char forceMask_local;
			// DetectSingleCell() = false for SingleCellDistinctPolicy and CellPairPolicy, true for SingleCellPolicy
			// we need this, since in contrast to sse3 we can no longer guarantee, that j>=i by default (j==i is handled by ForcePolicy::Condition).
			// however only one of the branches should be chosen by the compiler, since the class is known at compile time.
		if (ForcePolicy::DetectSingleCell()) {
			forceMask_local = (ForcePolicy::Condition(m_r2, cutoffRadiusSquare) && j > i_center_idx) ? 1 : 0;
		} else {
			forceMask_local = ForcePolicy::Condition(m_r2, cutoffRadiusSquare) ? 1 : 0;
		}
		forceMask += forceMask_local * bitmultiplier;

		//*(soa2_center_dist_lookup + j) = forceMask;
	}
	compute_molecule = vcp_simd_or(compute_molecule, forceMask);//from last iteration
	vcp_simd_store(soa2_center_dist_lookup + k, forceMask);//only one store, since there is only one additional mask.
	//no memset needed, since each element of soa2_center_dist_lookup corresponds to 8 masked elements.
	//every element of soa2_center_dist_lookup is therefore set.

	return compute_molecule;
#elif VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
	counter=0;
	static const __m512i eight = _mm512_set_epi32(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08
	);
	static const __m512i first_indices = _mm512_set_epi32(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00
	);
	static const __mmask8 possibleRemainderMasks[8] = { 0x00, 0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3F, 0x7F  };

	// distance and force mask computation
	//size_t i_center_idx = soa1._mol_num_ljc_acc[i];

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	__mmask8 initJ_mask = ForcePolicy :: InitJ_Mask(i_center_idx);

	__m512i indices = _mm512_mask_add_epi32(_mm512_setzero_epi32(), static_cast<__mmask16>(0x00FF), first_indices, _mm512_set1_epi32(j));
	const __mmask8 remainderMask = possibleRemainderMasks[soa2_num_centers & static_cast<size_t>(7)];

	for (; j < end_j; j += 8) {
		const __m512d m_r_x2 = _mm512_load_pd(soa2_m_r_x + j);
		const __m512d m_r_y2 = _mm512_load_pd(soa2_m_r_y + j);
		const __m512d m_r_z2 = _mm512_load_pd(soa2_m_r_z + j);

		const __m512d m_dx = _mm512_sub_pd(m1_r_x, m_r_x2);
		const __m512d m_dy = _mm512_sub_pd(m1_r_y, m_r_y2);
		const __m512d m_dz = _mm512_sub_pd(m1_r_z, m_r_z2);

		const __m512d m_dxdx = _mm512_mul_pd(m_dx, m_dx);
		const __m512d m_dxdx_dydy = _mm512_fmadd_pd(m_dy, m_dy, m_dxdx);
		const __m512d m_r2 = _mm512_fmadd_pd(m_dz, m_dz, m_dxdx_dydy);

		const __mmask8 forceMask = ForcePolicy :: GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);

		_mm512_mask_packstorelo_epi32(soa2_center_dist_lookup + counter, static_cast<__mmask16>(forceMask), indices);//these two lines are an unaligned store
		_mm512_mask_packstorehi_epi32(soa2_center_dist_lookup + counter + (64 / sizeof(vcp_lookupOrMask_single)), static_cast<__mmask16>(forceMask), indices);//these two lines are an unaligned store

		indices = _mm512_add_epi32(indices, eight);
		counter += _popcnt32(forceMask);
	}
	if (remainderMask != 0x00) {
		const __m512d m_r_x2 = _mm512_mask_load_pd(zero, remainderMask, soa2_m_r_x + j);
		const __m512d m_r_y2 = _mm512_mask_load_pd(zero, remainderMask, soa2_m_r_y + j);
		const __m512d m_r_z2 = _mm512_mask_load_pd(zero, remainderMask, soa2_m_r_z + j);

		const __m512d m_dx = _mm512_sub_pd(m1_r_x, m_r_x2);
		const __m512d m_dy = _mm512_sub_pd(m1_r_y, m_r_y2);
		const __m512d m_dz = _mm512_sub_pd(m1_r_z, m_r_z2);

		const __m512d m_dxdx = _mm512_mul_pd(m_dx, m_dx);
		const __m512d m_dxdx_dydy = _mm512_fmadd_pd(m_dy, m_dy, m_dxdx);
		const __m512d m_r2 = _mm512_fmadd_pd(m_dz, m_dz, m_dxdx_dydy);

		const __mmask8 forceMask = remainderMask & ForcePolicy :: GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);//AND remainderMask -> set unimportant ones to zero.

		_mm512_mask_packstorelo_epi32(soa2_center_dist_lookup + counter, static_cast<__mmask16>(forceMask), indices);//these two lines are an unaligned store
		_mm512_mask_packstorehi_epi32(soa2_center_dist_lookup + counter + (64 / sizeof(vcp_lookupOrMask_single)), static_cast<__mmask16>(forceMask), indices);//these two lines are an unaligned store

		indices = _mm512_add_epi32(indices, eight);
		counter += _popcnt32(forceMask);
	}

	return counter>0?VCP_SIMD_ONESVM:VCP_SIMD_ZEROVM;//do not compute stuff if nothing needs to be computed.
#endif
}

template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
void VectorizedChargeP2PCellProcessor :: _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
	// initialize dist lookups
	if(_centers_dist_lookup.get_size() < soa2._charges_size){
		soa2.resizeLastZero(_centers_dist_lookup, soa2._charges_size, soa2._charges_num);
	}
	_charges_dist_lookup = _centers_dist_lookup;

	// Pointer for molecules
	const double * const soa1_mol_pos_x = soa1._mol_pos.xBegin();
	const double * const soa1_mol_pos_y = soa1._mol_pos.yBegin();
	const double * const soa1_mol_pos_z = soa1._mol_pos.zBegin();

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

	vcp_lookupOrMask_single* const soa2_charges_dist_lookup = _charges_dist_lookup;


	vcp_double_vec sum_upotXpoles = VCP_SIMD_ZEROV;
	vcp_double_vec sum_virial = VCP_SIMD_ZEROV;

	const vcp_double_vec rc2 = vcp_simd_set1(_LJCutoffRadiusSquare);
	const vcp_double_vec cutoffRadiusSquare = vcp_simd_set1(_cutoffRadiusSquare);

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
	countertype32 end_charges_j_cnt = 0;//count for gather
	size_t i_charge_idx = 0;

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._mol_num; ++i) {//over the molecules
		const vcp_double_vec m1_r_x = vcp_simd_broadcast(soa1_mol_pos_x + i);
		const vcp_double_vec m1_r_y = vcp_simd_broadcast(soa1_mol_pos_y + i);
		const vcp_double_vec m1_r_z = vcp_simd_broadcast(soa1_mol_pos_z + i);
		// Iterate over centers of second cell
		const vcp_mask_vec compute_molecule_charges = calcDistLookup<ForcePolicy>(soa1, i, i_charge_idx, soa2._charges_num, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z, end_charges_j_cnt);

		size_t end_charges_loop = MaskGatherChooser::getEndloop(end_charges_j_longloop, end_charges_j_cnt);

		// Computation of site interactions with charges

		if (!vcp_simd_movemask(compute_molecule_charges)) {
			i_charge_idx += soa1_mol_charges_num[i];
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
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(end_charges_j_cnt);
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
			i_charge_idx += soa1_mol_charges_num[i];
		}


	}

	hSum_Add_Store(&_upotXpoles, sum_upotXpoles);
	hSum_Add_Store(&_virial, sum_virial);

} // void LennardJonesCellHandler::CalculatePairs_(LJSoA & soa1, LJSoA & soa2)

void VectorizedChargeP2PCellProcessor::processCell(ParticleCell & c) {
	if (c.isHaloCell() || (c.getCellDataSoA()._mol_num < 2)) {
		return;
	}
	_calculatePairs<SingleCellPolicy_, true, MaskGatherC>(c.getCellDataSoA(), c.getCellDataSoA());
}

void VectorizedChargeP2PCellProcessor::processCellPair(ParticleCell & c1, ParticleCell & c2) {
	assert(&c1 != &c2);

	if ((c1.getCellDataSoA()._mol_num == 0) || (c2.getCellDataSoA()._mol_num == 0)) {
		return;
	}
	if (!(c1.isHaloCell() || c2.isHaloCell())) {//no cell is halo
		_calculatePairs<CellPairPolicy_, true, MaskGatherC>(c1.getCellDataSoA(), c2.getCellDataSoA());
	} else if (c1.isHaloCell() == (!c2.isHaloCell())) {//exactly one cell is halo, therefore we only calculate some of the interactions.
		if (c1.getCellIndex() < c2.getCellIndex()){//using this method one can neglect the macroscopic boundary condition.
			_calculatePairs<CellPairPolicy_, true, MaskGatherC>(c1.getCellDataSoA(), c2.getCellDataSoA());
		}
		else {
			_calculatePairs<CellPairPolicy_, false, MaskGatherC>(c1.getCellDataSoA(), c2.getCellDataSoA());
		}
	} else {//both cells halo -> do nothing
		return;
	}
}

} // namespace bhfmm

