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
namespace bhfmm {
VectorizedLJP2PCellProcessor::VectorizedLJP2PCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		CellProcessor(cutoffRadius, LJcutoffRadius), _domain(domain),
		// maybe move the following to somewhere else:
		_epsRFInvrc3(2. * (domain.getepsilonRF() - 1.) / ((cutoffRadius * cutoffRadius * cutoffRadius) * (2. * domain.getepsilonRF() + 1.))), 
		_compIDs(), _eps_sig(), _shift6(), _upot6lj(0.0), _virial(0.0) {

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
					if ((components[comp_i].ID() == components[comp_j].ID())
							&& (components[comp_i].numTersoff() > 0
									|| components[comp_j].numTersoff() > 0)) {
						// No LJ interaction between solid atoms of the same component.
						_eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j)] = 0.0;
						_eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j) + 1] = 0.0;
						_shift6[_compIDs[comp_i] + center_i][_compIDs[comp_j] + center_j] = 0.0;
					} else {
						// Extract epsilon*24.0, sigma^2 and shift*6.0 from paramStreams.
						p >> _eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j)];
						p >> _eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j) + 1];
						p >> _shift6[_compIDs[comp_i] + center_i][_compIDs[comp_j] + center_j];
					}
				}
			}
		}
	}
}

VectorizedLJP2PCellProcessor :: ~VectorizedLJP2PCellProcessor () {
	for (size_t i = 0; i < _particleCellDataVector.size(); ++i) {
		delete _particleCellDataVector[i];
	}
	_particleCellDataVector.clear();
}


void VectorizedLJP2PCellProcessor::initTraversal(const size_t numCells) {
	_virial = 0.0;
	_upot6lj = 0.0;

	global_log->debug() << "VectorizedLJCellProcessor::initTraversal() to " << numCells << " cells." << std::endl;

	if (numCells > _particleCellDataVector.size()) {
		for (size_t i = _particleCellDataVector.size(); i < numCells; i++) {
			_particleCellDataVector.push_back(new CellDataSoA(64,64,64,64,64));
		}
		global_log->debug() << "resize CellDataSoA to " << numCells << " cells." << std::endl;
	}
}


void VectorizedLJP2PCellProcessor::endTraversal() {
	_domain.setLocalVirial(_virial /*+ 3.0 * _myRF*/);
	_domain.setLocalUpot(_upot6lj / 6.0 /*+ _upotXpoles + _myRF*/);
}


void VectorizedLJP2PCellProcessor::preprocessCell(ParticleCell & c) {
	assert(!c.getCellDataSoA());

	const MoleculeList & molecules = c.getParticlePointers();

	// Determine the total number of centers.
	size_t numMolecules = molecules.size();
	size_t nLJCenters = 0;
	
	for (size_t m = 0;  m < numMolecules; ++m) {
		nLJCenters += molecules[m]->numLJcenters();
	}

	// Construct the SoA.
	assert(!_particleCellDataVector.empty()); 
	CellDataSoA* soaPtr = _particleCellDataVector.back();
	CellDataSoA & soa = *soaPtr;
	soa.resize(numMolecules,nLJCenters,0,0,0);
	c.setCellDataSoA(soaPtr);
	_particleCellDataVector.pop_back();

	ComponentList components = *(_simulation.getEnsemble()->components());

	size_t iLJCenters = 0;
	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < molecules.size(); ++i) {
		const size_t mol_ljc_num = molecules[i]->numLJcenters();
		const double mol_pos_x = molecules[i]->r(0);
		const double mol_pos_y = molecules[i]->r(1);
		const double mol_pos_z = molecules[i]->r(2);

		soa._mol_pos_x[i] = mol_pos_x;
		soa._mol_pos_y[i] = mol_pos_y;
		soa._mol_pos_z[i] = mol_pos_z;
		soa._mol_ljc_num[i] = mol_ljc_num;

		for (size_t j = 0; j < mol_ljc_num; ++j, ++iLJCenters) {
			// Store a copy of the molecule position for each center, and the position of
			// each center. Assign each LJ center its ID and set the force to 0.0.
			soa._ljc_m_r_x[iLJCenters] = mol_pos_x;
			soa._ljc_m_r_y[iLJCenters] = mol_pos_y;
			soa._ljc_m_r_z[iLJCenters] = mol_pos_z;
			soa._ljc_r_x[iLJCenters] = molecules[i]->ljcenter_d(j)[0] + mol_pos_x;
			soa._ljc_r_y[iLJCenters] = molecules[i]->ljcenter_d(j)[1] + mol_pos_y;
			soa._ljc_r_z[iLJCenters] = molecules[i]->ljcenter_d(j)[2] + mol_pos_z;
			soa._ljc_f_x[iLJCenters] = 0.0;
			soa._ljc_f_y[iLJCenters] = 0.0;
			soa._ljc_f_z[iLJCenters] = 0.0;
			soa._ljc_V_x[iLJCenters] = 0.0;
			soa._ljc_V_y[iLJCenters] = 0.0;
			soa._ljc_V_z[iLJCenters] = 0.0;
			soa._ljc_id[iLJCenters] = _compIDs[molecules[i]->componentid()] + j;
			soa._ljc_dist_lookup[iLJCenters] = 0.0;
		}
	}
}


void VectorizedLJP2PCellProcessor::postprocessCell(ParticleCell & c) {
	assert(c.getCellDataSoA());
	CellDataSoA& soa = *c.getCellDataSoA();

	MoleculeList & molecules = c.getParticlePointers();

	// For each molecule iterate over all its centers.
	size_t iLJCenters = 0;
	size_t numMols = molecules.size();
	for (size_t m = 0; m < numMols; ++m) {
		const size_t mol_ljc_num = molecules[m]->numLJcenters();

		for (size_t i = 0; i < mol_ljc_num; ++i, ++iLJCenters) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa._ljc_f_x[iLJCenters];
			f[1] = soa._ljc_f_y[iLJCenters];
			f[2] = soa._ljc_f_z[iLJCenters];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			molecules[m]->Fljcenteradd(i, f);

			// Store the resulting virial in the molecule.
			double V[3];
			V[0] = soa._ljc_V_x[iLJCenters]*0.5;
			V[1] = soa._ljc_V_y[iLJCenters]*0.5;
			V[2] = soa._ljc_V_z[iLJCenters]*0.5;
			assert(!isnan(V[0]));
			assert(!isnan(V[1]));
			assert(!isnan(V[2]));
			molecules[m]->Viadd(V);
		}
	}
	// Delete the SoA.
	_particleCellDataVector.push_back(&soa);
	c.setCellDataSoA(0);
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


template<class ForcePolicy>
vcp_mask_vec
inline VectorizedLJP2PCellProcessor::calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		vcp_lookupOrMask_single* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z,
		const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z,
		countertype32& counter
		) {

#if VCP_VEC_TYPE==VCP_NOVEC

	bool compute_molecule = false;

	for (size_t j = ForcePolicy :: InitJ(i_center_idx); j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];
		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		const bool forceMask = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? (~0l) : 0l;
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

		const vcp_double_vec m_r2 = vcp_simd_fma(m_dx, m_dx, vcp_simd_fma(m_dy, m_dy, m_dz * m_dz));

		const vcp_mask_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD);
		vcp_simd_store(soa2_center_dist_lookup + j, forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask);
	}

	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];

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
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];

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
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];

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
void VectorizedLJP2PCellProcessor :: _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
	// Pointer for molecules
	const double * const soa1_mol_pos_x = soa1._mol_pos_x;
	const double * const soa1_mol_pos_y = soa1._mol_pos_y;
	const double * const soa1_mol_pos_z = soa1._mol_pos_z;

	// Pointer for LJ centers
	const double * const soa1_ljc_r_x = soa1._ljc_r_x;
	const double * const soa1_ljc_r_y = soa1._ljc_r_y;
	const double * const soa1_ljc_r_z = soa1._ljc_r_z;
	double * const soa1_ljc_f_x = soa1._ljc_f_x;
	double * const soa1_ljc_f_y = soa1._ljc_f_y;
	double * const soa1_ljc_f_z = soa1._ljc_f_z;
	double * const soa1_ljc_V_x = soa1._ljc_V_x;
	double * const soa1_ljc_V_y = soa1._ljc_V_y;
	double * const soa1_ljc_V_z = soa1._ljc_V_z;
	const int * const soa1_mol_ljc_num = soa1._mol_ljc_num;
	const size_t * const soa1_ljc_id = soa1._ljc_id;

	const double * const soa2_ljc_m_r_x = soa2._ljc_m_r_x;
	const double * const soa2_ljc_m_r_y = soa2._ljc_m_r_y;
	const double * const soa2_ljc_m_r_z = soa2._ljc_m_r_z;
	const double * const soa2_ljc_r_x = soa2._ljc_r_x;
	const double * const soa2_ljc_r_y = soa2._ljc_r_y;
	const double * const soa2_ljc_r_z = soa2._ljc_r_z;
	double * const soa2_ljc_f_x = soa2._ljc_f_x;
	double * const soa2_ljc_f_y = soa2._ljc_f_y;
	double * const soa2_ljc_f_z = soa2._ljc_f_z;
	double * const soa2_ljc_V_x = soa2._ljc_V_x;
	double * const soa2_ljc_V_y = soa2._ljc_V_y;
	double * const soa2_ljc_V_z = soa2._ljc_V_z;
	const size_t * const soa2_ljc_id = soa2._ljc_id;

	vcp_lookupOrMask_single* const soa2_ljc_dist_lookup = soa2._ljc_dist_lookup;

	vcp_double_vec sum_upot6lj = VCP_SIMD_ZEROV;
	vcp_double_vec sum_virial = VCP_SIMD_ZEROV;

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
	countertype32 end_ljc_j_cnt = 0;//count for gather
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
		const vcp_mask_vec compute_molecule_ljc = calcDistLookup<ForcePolicy>(soa1, i, i_ljc_idx, soa2._ljc_num, _LJCutoffRadiusSquare,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z, end_ljc_j_cnt);

		size_t end_ljc_loop = MaskGatherChooser::getEndloop(end_ljc_j_longloop, end_ljc_j_cnt);

		if (!vcp_simd_movemask(compute_molecule_ljc)) {
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
	assert(c.getCellDataSoA());
	if (c.isHaloCell() || (c.getCellDataSoA()->_mol_num < 2)) {
		return;
	}
//printf("---------------singleCell-------------\n");
	_calculatePairs<SingleCellPolicy_, true, MaskGatherC>(*(c.getCellDataSoA()), *(c.getCellDataSoA()));
}

void VectorizedLJP2PCellProcessor::processCellPair(ParticleCell & c1, ParticleCell & c2) {
	assert(&c1 != &c2);
	assert(c1.getCellDataSoA());
	assert(c2.getCellDataSoA());

	if ((c1.getCellDataSoA()->_mol_num == 0) || (c2.getCellDataSoA()->_mol_num == 0)) {
		return;
	}
	//printf("---------------cellPair-------------\n");
	if (!(c1.isHaloCell() || c2.isHaloCell())) {//no cell is halo
		_calculatePairs<CellPairPolicy_, true, MaskGatherC>(*(c1.getCellDataSoA()), *(c2.getCellDataSoA()));
	} else if (c1.isHaloCell() == (!c2.isHaloCell())) {//exactly one cell is halo, therefore we only calculate some of the interactions.
		if (c1.getCellIndex() < c2.getCellIndex()){//using this method one can neglect the macroscopic boundary condition.
			_calculatePairs<CellPairPolicy_, true, MaskGatherC>(*(c1.getCellDataSoA()), *(c2.getCellDataSoA()));
		}
		else {
			_calculatePairs<CellPairPolicy_, false, MaskGatherC>(*(c1.getCellDataSoA()), *(c2.getCellDataSoA()));
		}
	} else {//both cells halo -> do nothing
		return;
	}
}

} // namespace bhfmm

