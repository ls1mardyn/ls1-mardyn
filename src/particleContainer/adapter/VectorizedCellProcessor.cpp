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

using namespace Log;

VectorizedCellProcessor::VectorizedCellProcessor(Domain & domain,
		double cutoffRadius) :
		_domain(domain), _rc2(cutoffRadius * cutoffRadius), _compIDs(), _eps_sig(), _shift6(), _upot6lj(
				0.0), _virial(0.0), _center_dist_lookup(128) {
#if VLJCP_VEC_TYPE==VLJCP_NOVEC
	Log::global_log->info() << "VectorizedLJCellProcessor: no vectorization."
	<< std::endl;
#elif VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
	Log::global_log->info() << "VectorizedLJCellProcessor: using SSE3."
			<< std::endl;
#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX
	Log::global_log->info() << "VectorizedLJCellProcessor: using AVX."
			<< std::endl;
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

VectorizedCellProcessor :: ~VectorizedCellProcessor () {
	for (size_t i = 0; i < _particleCellDataVector.size(); ++i) {
		delete _particleCellDataVector[i];
	}
	_particleCellDataVector.clear();
}


void VectorizedCellProcessor::initTraversal(const size_t numCells) {
	_virial = 0.0;
	_upot6lj = 0.0;

	global_log->debug() << "VectorizedLJCellProcessor::initTraversal() to " << numCells << " cells." << std::endl;

	if (numCells > _particleCellDataVector.size()) {
//		_particleCellDataVector.resize(numCells);
		for (size_t i = _particleCellDataVector.size(); i < numCells; i++) {
			_particleCellDataVector.push_back(new CellDataSoA(64,64));
		}
		global_log->debug() << "resize CellDataSoA to " << numCells << " cells." << std::endl;
	}

}


void VectorizedCellProcessor::endTraversal() {
	_domain.setLocalVirial(_virial);
	_domain.setLocalUpot(_upot6lj / 6.0);
}


void VectorizedCellProcessor::preprocessCell(ParticleCell & c) {
	assert(!c.getCellDataSoA());

	const MoleculeList & molecules = c.getParticlePointers();

	// Determine the total number of LJ centers.
	size_t numMolecules = molecules.size();
	size_t nLJCenters = 0;
	for (size_t m = 0;  m < numMolecules; ++m) {
		nLJCenters += molecules[m]->numLJcenters();
	}

	// Construct the SoA.
	assert(!_particleCellDataVector.empty());
	CellDataSoA* soaPtr = _particleCellDataVector.back();
	global_log->debug() << " _particleCellDataVector.size()=" << _particleCellDataVector.size() << " soaPtr=" << soaPtr << " nLJCenters=" << nLJCenters << std::endl;
	soaPtr->resize(numMolecules,nLJCenters);
	soaPtr->_num_ljcenters = nLJCenters;
	soaPtr->_num_molecules = numMolecules;
	c.setCellDataSoA(soaPtr);
	_particleCellDataVector.pop_back();
	CellDataSoA & soa = *soaPtr;

	size_t n = 0;
	// For each molecule iterate over all its LJ centers.
	for (size_t i = 0; i < molecules.size(); ++i) {
		const size_t nLJC = molecules[i]->numLJcenters();
		const double mol_pos_x = molecules[i]->r(0);
		const double mol_pos_y = molecules[i]->r(1);
		const double mol_pos_z = molecules[i]->r(2);

		soa._mol_pos_x[i] = mol_pos_x;
		soa._mol_pos_y[i] = mol_pos_y;
		soa._mol_pos_z[i] = mol_pos_z;
		soa._mol_num_ljc[i] = nLJC;

		for (size_t j = 0; j < nLJC; ++j, ++n) {
			// Store a copy of the molecule position for each center, and the position of
			// each LJ center. Assign each LJ center its ID and set the force to 0.0.
			soa._m_r_x[n] = mol_pos_x;
			soa._m_r_y[n] = mol_pos_y;
			soa._m_r_z[n] = mol_pos_z;
			soa._ljc_r_x[n] = molecules[i]->ljcenter_d(j)[0] + mol_pos_x;
			soa._ljc_r_y[n] = molecules[i]->ljcenter_d(j)[1] + mol_pos_y;
			soa._ljc_r_z[n] = molecules[i]->ljcenter_d(j)[2] + mol_pos_z;
			soa._ljc_f_x[n] = 0.0;
			soa._ljc_f_y[n] = 0.0;
			soa._ljc_f_z[n] = 0.0;
			soa._ljc_id[n] = _compIDs[molecules[i]->componentid()] + j;
		}
	}

	if (_center_dist_lookup.get_size() < nLJCenters) {
		_center_dist_lookup.resize(nLJCenters);
	}
}


void VectorizedCellProcessor::postprocessCell(ParticleCell & c) {
	assert(c.getCellDataSoA());
	CellDataSoA& soa = *c.getCellDataSoA();

	MoleculeList & molecules = c.getParticlePointers();

	// For each molecule iterate over all its centers.
	size_t n = 0;
	size_t numMols = molecules.size();
	for (size_t m = 0; m < numMols; ++m) {
		const size_t end = molecules[m]->numLJcenters();
		for (size_t i = 0; i < end; ++i) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa._ljc_f_x[n];
			f[1] = soa._ljc_f_y[n];
			f[2] = soa._ljc_f_z[n];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			molecules[m]->Fljcenteradd(i, f);
			++n;
		}
		molecules[m]->calcFM();
	}
	// Delete the SoA.
	_particleCellDataVector.push_back(&soa);
	c.setCellDataSoA(0);
}

template<class ForcePolicy, class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovec (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {
		// Distance vector from center 2 to center 1. Together with the rest of the
		// calculation, we get the proper force from center 1 to center 2.
		// This is done because the calculation of the potential fits in nicely that way.
		const double c_dx = soa1._ljc_r_x[i] - soa2._ljc_r_x[j];
		const double c_dy = soa1._ljc_r_y[i] - soa2._ljc_r_y[j];
		const double c_dz = soa1._ljc_r_z[i] - soa2._ljc_r_z[j];

		const double c_r2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;
		const double r2_inv = 1.0 / c_r2;

		const double eps_24 = _eps_sig[soa1._ljc_id[i]][2 * soa2._ljc_id[j]];
		const double sig2 = _eps_sig[soa1._ljc_id[i]][2 * soa2._ljc_id[j] + 1];

		const double lj2 = sig2 * r2_inv;
		const double lj6 = lj2 * lj2 * lj2;
		const double lj12 = lj6 * lj6;
		const double lj12m6 = lj12 - lj6;

		const double scale = eps_24 * r2_inv * (lj12 + lj12m6);

		const double fx = c_dx * scale;
		const double fy = c_dy * scale;
		const double fz = c_dz * scale;

		const double m_dx = soa1._m_r_x[i] - soa2._m_r_x[j];
		const double m_dy = soa1._m_r_y[i] - soa2._m_r_y[j];
		const double m_dz = soa1._m_r_z[i] - soa2._m_r_z[j];

		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueCondition(m_dx, m_dy, m_dz)) {
			_upot6lj += eps_24 * lj12m6 + _shift6[soa1._ljc_id[i]][soa2._ljc_id[j]];
			_virial += m_dx * fx + m_dy * fy + m_dz * fz;
		}
		// Add the force to center 1, and subtract it from center 2.
		soa1._ljc_f_x[i] += fx;
		soa1._ljc_f_y[i] += fy;
		soa1._ljc_f_z[i] += fz;

		soa2._ljc_f_x[j] -= fx;
		soa2._ljc_f_y[j] -= fy;
		soa2._ljc_f_z[j] -= fz;
	}
}  /* end of method VectorizedLJCellProcessor :: _loopBodyNovec */

template<class ForcePolicy, class MacroPolicy>
void VectorizedCellProcessor::_calculatePairs(const CellDataSoA & soa1,
		const CellDataSoA & soa2) {
#if VLJCP_VEC_TYPE==VLJCP_NOVEC
	// For the unvectorized version, we only have to iterate over all pairs of
	// LJ centers and apply the unvectorized loop body.
	size_t i_center_idx = 0;
	assert(_center_dist_lookup.get_size() >= soa2._ljcenters_size);
	for (size_t i = 0; i < soa1._num_molecules; ++i) {

		unsigned long compute_molecule = 0;
		for (size_t j = ForcePolicy :: InitJ(i_center_idx); j < soa2._num_ljcenters; ++j) {
			const double m_dx = soa1._mol_pos_x[i] - soa2._m_r_x[j];
			const double m_dy = soa1._mol_pos_y[i] - soa2._m_r_y[j];
			const double m_dz = soa1._mol_pos_z[i] - soa2._m_r_z[j];
			const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

			const signed long forceMask = ForcePolicy :: Condition(m_r2, _rc2) ? (~0l) : 0l;
			compute_molecule |= forceMask;
			*(_center_dist_lookup + j) = forceMask;
		}

		if (!compute_molecule) {
			i_center_idx += soa1._mol_num_ljc[i];
			continue;
		}

		for (int local_i = 0; local_i < soa1._mol_num_ljc[i]; local_i++ ) {
			for (size_t j = ForcePolicy :: InitJ(i_center_idx); j < soa2._num_ljcenters; ++j) {
				_loopBodyNovec<CellPairPolicy_, MacroPolicy>(soa1, i_center_idx, soa2, j, _center_dist_lookup + j);
			}
			i_center_idx++;
		}
	}


#elif VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
	const double * const p_mol_rx1 = soa1._mol_pos_x;
	const double * const p_mol_ry1 = soa1._mol_pos_y;
	const double * const p_mol_rz1 = soa1._mol_pos_z;
	const double * const p_crx1 = soa1._ljc_r_x;
	const double * const p_cry1 = soa1._ljc_r_y;
	const double * const p_crz1 = soa1._ljc_r_z;
	double * const p_cfx1 = soa1._ljc_f_x;
	double * const p_cfy1 = soa1._ljc_f_y;
	double * const p_cfz1 = soa1._ljc_f_z;
	const size_t * const p_cid1 = soa1._ljc_id;

	const double * const p_mrx2 = soa2._m_r_x;
	const double * const p_mry2 = soa2._m_r_y;
	const double * const p_mrz2 = soa2._m_r_z;
	const double * const p_crx2 = soa2._ljc_r_x;
	const double * const p_cry2 = soa2._ljc_r_y;
	const double * const p_crz2 = soa2._ljc_r_z;
	double * const p_cfx2 = soa2._ljc_f_x;
	double * const p_cfy2 = soa2._ljc_f_y;
	double * const p_cfz2 = soa2._ljc_f_z;
	const size_t * const p_cid2 = soa2._ljc_id;

	double* const p_center_dist_lookup = _center_dist_lookup;
	const size_t end_j = soa2._num_ljcenters & (~1);
	const __m128d one = _mm_set1_pd(1.0);
	const __m128d rc2 = _mm_set1_pd(_rc2);
	__m128d sum_upot = _mm_setzero_pd();
	__m128d sum_virial = _mm_setzero_pd();

	size_t i_center_idx = 0;
	assert(_center_dist_lookup.get_size() >= soa2._ljcenters_size);
	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._num_molecules; ++i) {
		const __m128d m_r_x1 = _mm_loaddup_pd(p_mol_rx1 + i);
		const __m128d m_r_y1 = _mm_loaddup_pd(p_mol_ry1 + i);
		const __m128d m_r_z1 = _mm_loaddup_pd(p_mol_rz1 + i);

		// distance and forca mask computation
		__m128d compute_molecule = _mm_setzero_pd();
		size_t j = ForcePolicy :: InitJ(i_center_idx);
		for (; j < end_j; j+=2) {
			const __m128d m_r_x2 = _mm_load_pd(p_mrx2 + j);
			const __m128d m_dx = _mm_sub_pd(m_r_x1, m_r_x2);
			const __m128d m_r_y2 = _mm_load_pd(p_mry2 + j);
			const __m128d m_dy = _mm_sub_pd(m_r_y1, m_r_y2);
			const __m128d m_r_z2 = _mm_load_pd(p_mrz2 + j);
			const __m128d m_dz = _mm_sub_pd(m_r_z1, m_r_z2);
			const __m128d m_dxdx = _mm_mul_pd(m_dx, m_dx);
			const __m128d m_dydy = _mm_mul_pd(m_dy, m_dy);
			const __m128d m_dzdz = _mm_mul_pd(m_dz, m_dz);
			const __m128d m_dxdx_dydy = _mm_add_pd(m_dxdx, m_dydy);
			const __m128d m_r2 = _mm_add_pd(m_dxdx_dydy, m_dzdz);
			const __m128d forceMask = ForcePolicy::GetForceMask(m_r2, rc2);
			_mm_store_pd(p_center_dist_lookup + j, forceMask);
			compute_molecule = _mm_or_pd(compute_molecule, forceMask);
		}
		for (; j < soa2._num_ljcenters; ++j) {
			const double m_dx = soa1._mol_pos_x[i] - soa2._m_r_x[j];
			const double m_dy = soa1._mol_pos_y[i] - soa2._m_r_y[j];
			const double m_dz = soa1._mol_pos_z[i] - soa2._m_r_z[j];

			const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;
			const signed long forceMask_l = ForcePolicy :: Condition(m_r2, _rc2) ? ~0l : 0l;
			// this casting via void* is required for gcc
			const void* forceMask_tmp = reinterpret_cast<const void*>(&forceMask_l);
			double forceMask = *reinterpret_cast<double const* const>(forceMask_tmp);

			*(p_center_dist_lookup + j) = forceMask;
			const __m128d forceMask_128 = _mm_set1_pd(forceMask);
			compute_molecule = _mm_or_pd(compute_molecule, forceMask_128);
		}

		if (!_mm_movemask_pd(compute_molecule)) {
			i_center_idx += soa1._mol_num_ljc[i];
			continue;
		}



		// actual force computation
		for (int local_i = 0; local_i < soa1._mol_num_ljc[i]; local_i++ ) {
			__m128d sum_fx1 = _mm_setzero_pd();
			__m128d sum_fy1 = _mm_setzero_pd();
			__m128d sum_fz1 = _mm_setzero_pd();
			const __m128d c_r_x1 = _mm_loaddup_pd(p_crx1 + i_center_idx);
			const __m128d c_r_y1 = _mm_loaddup_pd(p_cry1 + i_center_idx);
			const __m128d c_r_z1 = _mm_loaddup_pd(p_crz1 + i_center_idx);
			// Iterate over each pair of centers in the second cell.
			size_t j = ForcePolicy::InitJ(i_center_idx);
			for (; j < end_j; j += 2) {
				const __m128d forceMask = _mm_load_pd(p_center_dist_lookup + j);
				// Only go on if at least 1 of the forces has to be calculated.
				if (_mm_movemask_pd(forceMask) > 0) {
					const __m128d c_r_x2 = _mm_load_pd(p_crx2 + j);
					const __m128d c_dx = _mm_sub_pd(c_r_x1, c_r_x2);
					const __m128d c_r_y2 = _mm_load_pd(p_cry2 + j);
					const __m128d c_dy = _mm_sub_pd(c_r_y1, c_r_y2);
					const __m128d c_r_z2 = _mm_load_pd(p_crz2 + j);
					const __m128d c_dz = _mm_sub_pd(c_r_z1, c_r_z2);
					const __m128d c_dxdx = _mm_mul_pd(c_dx, c_dx);
					const __m128d c_dydy = _mm_mul_pd(c_dy, c_dy);
					const __m128d c_dzdz = _mm_mul_pd(c_dz, c_dz);
					const __m128d c_dxdx_dydy = _mm_add_pd(c_dxdx, c_dydy);
					const __m128d c_r2 = _mm_add_pd(c_dxdx_dydy, c_dzdz);
					const __m128d r2_inv_unmasked = _mm_div_pd(one, c_r2);
					const __m128d r2_inv = _mm_and_pd(r2_inv_unmasked, forceMask);
					const size_t id_i = p_cid1[i_center_idx];
					const size_t id_j0 = p_cid2[j];
					const size_t id_j1 = p_cid2[j + 1];
					const __m128d e1s1 = _mm_load_pd(_eps_sig[id_i] + 2 * id_j0);
					const __m128d e2s2 = _mm_load_pd(_eps_sig[id_i] + 2 * id_j1);
					const __m128d eps_24 = _mm_unpacklo_pd(e1s1, e2s2);
					const __m128d sig2 = _mm_unpackhi_pd(e1s1, e2s2);
					const __m128d lj2 = _mm_mul_pd(sig2, r2_inv);
					const __m128d lj4 = _mm_mul_pd(lj2, lj2);
					const __m128d lj6 = _mm_mul_pd(lj4, lj2);
					const __m128d lj12 = _mm_mul_pd(lj6, lj6);
					const __m128d lj12m6 = _mm_sub_pd(lj12, lj6);
					const __m128d eps24r2inv = _mm_mul_pd(eps_24, r2_inv);
					const __m128d lj12lj12m6 = _mm_add_pd(lj12, lj12m6);
					const __m128d scale = _mm_mul_pd(eps24r2inv, lj12lj12m6);
					const __m128d fx = _mm_mul_pd(c_dx, scale);
					const __m128d fy = _mm_mul_pd(c_dy, scale);
					const __m128d fz = _mm_mul_pd(c_dz, scale);

					const __m128d m_r_x2 = _mm_load_pd(p_mrx2 + j);
					const __m128d m_dx = _mm_sub_pd(m_r_x1, m_r_x2);
					const __m128d m_r_y2 = _mm_load_pd(p_mry2 + j);
					const __m128d m_dy = _mm_sub_pd(m_r_y1, m_r_y2);
					const __m128d m_r_z2 = _mm_load_pd(p_mrz2 + j);
					const __m128d m_dz = _mm_sub_pd(m_r_z1, m_r_z2);

					const __m128d macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
					// Only go on if at least 1 macroscopic value has to be calculated.
					if (_mm_movemask_pd(macroMask) > 0) {
						const __m128d sh1 = _mm_load_sd(_shift6[id_i] + id_j0);
						const __m128d sh2 = _mm_load_sd(_shift6[id_i] + id_j1);
						const __m128d shift6 = _mm_unpacklo_pd(sh1, sh2);
						const __m128d upot = _mm_mul_pd(eps_24, lj12m6);
						const __m128d upot_sh = _mm_add_pd(shift6, upot);
						const __m128d upot_masked = _mm_and_pd(upot_sh, macroMask);
						sum_upot = _mm_add_pd(sum_upot, upot_masked);
						const __m128d vir_x = _mm_mul_pd(m_dx, fx);
						const __m128d vir_y = _mm_mul_pd(m_dy, fy);
						const __m128d vir_z = _mm_mul_pd(m_dz, fz);
						const __m128d vir_xy = _mm_add_pd(vir_x, vir_y);
						const __m128d virial = _mm_add_pd(vir_xy, vir_z);
						const __m128d vir_masked = _mm_and_pd(virial, macroMask);
						sum_virial = _mm_add_pd(sum_virial, vir_masked);
					}
					const __m128d old_fx2 = _mm_load_pd(p_cfx2 + j);
					const __m128d new_fx2 = _mm_sub_pd(old_fx2, fx);
					_mm_store_pd(p_cfx2 + j, new_fx2);
					const __m128d old_fy2 = _mm_load_pd(p_cfy2 + j);
					const __m128d new_fy2 = _mm_sub_pd(old_fy2, fy);
					_mm_store_pd(p_cfy2 + j, new_fy2);
					const __m128d old_fz2 = _mm_load_pd(p_cfz2 + j);
					const __m128d new_fz2 = _mm_sub_pd(old_fz2, fz);
					_mm_store_pd(p_cfz2 + j, new_fz2);
					sum_fx1 = _mm_add_pd(sum_fx1, fx);
					sum_fy1 = _mm_add_pd(sum_fy1, fy);
					sum_fz1 = _mm_add_pd(sum_fz1, fz);
				}
			}
			_mm_store_sd(
					p_cfx1 + i_center_idx,
					_mm_add_sd(_mm_hadd_pd(sum_fx1, sum_fx1),
							_mm_load_sd(p_cfx1 + i_center_idx)));
			_mm_store_sd(
					p_cfy1 + i_center_idx,
					_mm_add_sd(_mm_hadd_pd(sum_fy1, sum_fy1),
							_mm_load_sd(p_cfy1 + i_center_idx)));
			_mm_store_sd(
					p_cfz1 + i_center_idx,
					_mm_add_sd(_mm_hadd_pd(sum_fz1, sum_fz1),
							_mm_load_sd(p_cfz1 + i_center_idx)));

			// Unvectorized calculation for leftover pairs.
			switch (soa2._num_ljcenters & 1) {
				case 1: {
					_loopBodyNovec<ForcePolicy, MacroPolicy>(soa1, i_center_idx, soa2, end_j, p_center_dist_lookup + j);
				}
				break;
			}

			i_center_idx++;
		}
	}
	_mm_store_sd(
			&_upot6lj,
			_mm_add_sd(_mm_hadd_pd(sum_upot, sum_upot),
					_mm_load_sd(&_upot6lj)));
	_mm_store_sd(
			&_virial,
			_mm_add_sd(_mm_hadd_pd(sum_virial, sum_virial),
					_mm_load_sd(&_virial)));

#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX

	const double * const p_mol_rx1 = soa1._mol_pos_x;
	const double * const p_mol_ry1 = soa1._mol_pos_y;
	const double * const p_mol_rz1 = soa1._mol_pos_z;
	const double * const p_crx1 = soa1._ljc_r_x;
	const double * const p_cry1 = soa1._ljc_r_y;
	const double * const p_crz1 = soa1._ljc_r_z;
	double * const p_cfx1 = soa1._ljc_f_x;
	double * const p_cfy1 = soa1._ljc_f_y;
	double * const p_cfz1 = soa1._ljc_f_z;
	const size_t * const p_cid1 = soa1._ljc_id;

	const double * const p_mrx2 = soa2._m_r_x;
	const double * const p_mry2 = soa2._m_r_y;
	const double * const p_mrz2 = soa2._m_r_z;
	const double * const p_crx2 = soa2._ljc_r_x;
	const double * const p_cry2 = soa2._ljc_r_y;
	const double * const p_crz2 = soa2._ljc_r_z;
	double * const p_cfx2 = soa2._ljc_f_x;
	double * const p_cfy2 = soa2._ljc_f_y;
	double * const p_cfz2 = soa2._ljc_f_z;
	const size_t * const p_cid2 = soa2._ljc_id;

	double* const p_center_dist_lookup = _center_dist_lookup;
	const size_t end_j = soa2._num_ljcenters & ~static_cast<size_t>(3);
	const __m256d one = _mm256_set1_pd(1.0);
	const __m256d rc2 = _mm256_set1_pd(_rc2);

	__m256d sum_upot = _mm256_setzero_pd();
	__m256d sum_virial = _mm256_setzero_pd();

	static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
	static const __m256i memoryMask_first_second = _mm256_set_epi32(0, 0, 0, 0, 1<<31, 0, 1<<31, 0);

	size_t i_center_idx = 0;
	assert(_center_dist_lookup.get_size() >= soa2._num_ljcenters);

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._num_molecules; ++i) {
		const __m256d m_r_x1 = _mm256_broadcast_sd(p_mol_rx1 + i);
		const __m256d m_r_y1 = _mm256_broadcast_sd(p_mol_ry1 + i);
		const __m256d m_r_z1 = _mm256_broadcast_sd(p_mol_rz1 + i);

		// distance and forca mask computation
		__m256d compute_molecule = _mm256_setzero_pd();
		size_t j = ForcePolicy :: InitJ(i_center_idx);
		__m256d initJ_mask = ForcePolicy::InitJ_Mask(i_center_idx);
		for (; j < end_j; j+=4) {
			const __m256d m_r_x2 = _mm256_load_pd(p_mrx2 + j);
			const __m256d m_dx = _mm256_sub_pd(m_r_x1, m_r_x2);
			const __m256d m_r_y2 = _mm256_load_pd(p_mry2 + j);
			const __m256d m_dy = _mm256_sub_pd(m_r_y1, m_r_y2);
			const __m256d m_r_z2 = _mm256_load_pd(p_mrz2 + j);
			const __m256d m_dz = _mm256_sub_pd(m_r_z1, m_r_z2);
			const __m256d m_dxdx = _mm256_mul_pd(m_dx, m_dx);
			const __m256d m_dydy = _mm256_mul_pd(m_dy, m_dy);
			const __m256d m_dzdz = _mm256_mul_pd(m_dz, m_dz);
			const __m256d m_dxdx_dydy = _mm256_add_pd(m_dxdx, m_dydy);
			const __m256d m_r2 = _mm256_add_pd(m_dxdx_dydy, m_dzdz);
			const __m256d forceMask = ForcePolicy::GetForceMask(m_r2, rc2, initJ_mask);
			_mm256_store_pd(p_center_dist_lookup + j, forceMask);
			compute_molecule = _mm256_or_pd(compute_molecule, forceMask);
		}
		for (; j < soa2._num_ljcenters; ++j) {
			const double m_dx = soa1._mol_pos_x[i] - soa2._m_r_x[j];
			const double m_dy = soa1._mol_pos_y[i] - soa2._m_r_y[j];
			const double m_dz = soa1._mol_pos_z[i] - soa2._m_r_z[j];

			const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;
			signed long forceMask_l;
			if (ForcePolicy::DetectSingleCell()) {
				forceMask_l = (ForcePolicy::Condition(m_r2, _rc2) && j > i_center_idx) ? ~0l : 0l;
			} else {
				forceMask_l = ForcePolicy::Condition(m_r2, _rc2) ? ~0l : 0l;
			}

//			this casting via void* is required for gcc
			void* forceMask_tmp = reinterpret_cast<void*>(&forceMask_l);
			double forceMask = *reinterpret_cast<double const* const>(forceMask_tmp);

			*(p_center_dist_lookup + j) = forceMask;
			const __m256d forceMask_256 = _mm256_set1_pd(forceMask);
			compute_molecule = _mm256_or_pd(compute_molecule, forceMask_256);
		}

		if (!_mm256_movemask_pd(compute_molecule)) {
			i_center_idx += soa1._mol_num_ljc[i];
			continue;
		}

		// actual force computation
		for (int local_i = 0; local_i < soa1._mol_num_ljc[i]; local_i++ ) {
			__m256d sum_fx1 = _mm256_setzero_pd();
			__m256d sum_fy1 = _mm256_setzero_pd();
			__m256d sum_fz1 = _mm256_setzero_pd();
			const __m256d c_r_x1 = _mm256_broadcast_sd(p_crx1 + i_center_idx);
			const __m256d c_r_y1 = _mm256_broadcast_sd(p_cry1 + i_center_idx);
			const __m256d c_r_z1 = _mm256_broadcast_sd(p_crz1 + i_center_idx);
			// Iterate over each pair of centers in the second cell.
			size_t j = ForcePolicy::InitJ(i_center_idx);
			for (; j < end_j; j += 4) {
				const __m256d forceMask = _mm256_load_pd(p_center_dist_lookup + j);
				// Only go on if at least 1 of the forces has to be calculated.
				if (_mm256_movemask_pd(forceMask) > 0) {
					const __m256d c_r_x2 = _mm256_load_pd(p_crx2 + j);
					const __m256d c_dx = _mm256_sub_pd(c_r_x1, c_r_x2);
					const __m256d c_r_y2 = _mm256_load_pd(p_cry2 + j);
					const __m256d c_dy = _mm256_sub_pd(c_r_y1, c_r_y2);
					const __m256d c_r_z2 = _mm256_load_pd(p_crz2 + j);
					const __m256d c_dz = _mm256_sub_pd(c_r_z1, c_r_z2);
					const __m256d c_dxdx = _mm256_mul_pd(c_dx, c_dx);
					const __m256d c_dydy = _mm256_mul_pd(c_dy, c_dy);
					const __m256d c_dzdz = _mm256_mul_pd(c_dz, c_dz);
					const __m256d c_dxdx_dydy = _mm256_add_pd(c_dxdx, c_dydy);
					const __m256d c_r2 = _mm256_add_pd(c_dxdx_dydy, c_dzdz);
					const __m256d r2_inv_unmasked = _mm256_div_pd(one, c_r2);
					const __m256d r2_inv = _mm256_and_pd(r2_inv_unmasked, forceMask);

					const size_t id_i = p_cid1[i_center_idx];
					const size_t id_j0 = p_cid2[j];
					const size_t id_j1 = p_cid2[j + 1];
					const size_t id_j2 = p_cid2[j + 2];
					const size_t id_j3 = p_cid2[j + 3];

					const __m256d e0s0 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j0, memoryMask_first_second);
					const __m256d e1s1 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j1, memoryMask_first_second);
					const __m256d e2s2 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j2, memoryMask_first_second);
					const __m256d e3s3 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j3, memoryMask_first_second);

					const __m256d e0e1 = _mm256_unpacklo_pd(e0s0, e1s1);
					const __m256d s0s1 = _mm256_unpackhi_pd(e0s0, e1s1);
					const __m256d e2e3 = _mm256_unpacklo_pd(e2s2, e3s3);
					const __m256d s2s3 = _mm256_unpackhi_pd(e2s2, e3s3);

					const __m256d eps_24 = _mm256_permute2f128_pd(e0e1, e2e3, 1<<5);
					const __m256d sig2 = _mm256_permute2f128_pd(s0s1, s2s3, 1<<5);

					const __m256d lj2 = _mm256_mul_pd(sig2, r2_inv);
					const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
					const __m256d lj6 = _mm256_mul_pd(lj4, lj2);
					const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
					const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);

					const __m256d eps24r2inv = _mm256_mul_pd(eps_24, r2_inv);
					const __m256d lj12lj12m6 = _mm256_add_pd(lj12, lj12m6);
					const __m256d scale = _mm256_mul_pd(eps24r2inv, lj12lj12m6);

					const __m256d fx = _mm256_mul_pd(c_dx, scale);
					const __m256d fy = _mm256_mul_pd(c_dy, scale);
					const __m256d fz = _mm256_mul_pd(c_dz, scale);

					const __m256d m_r_x2 = _mm256_load_pd(p_mrx2 + j);
					const __m256d m_dx = _mm256_sub_pd(m_r_x1, m_r_x2);
					const __m256d m_r_y2 = _mm256_load_pd(p_mry2 + j);
					const __m256d m_dy = _mm256_sub_pd(m_r_y1, m_r_y2);
					const __m256d m_r_z2 = _mm256_load_pd(p_mrz2 + j);
					const __m256d m_dz = _mm256_sub_pd(m_r_z1, m_r_z2);

					const __m256d macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);

					// Only go on if at least 1 macroscopic value has to be calculated.
					if (_mm256_movemask_pd(macroMask) > 0) {
						const __m256d sh0 = _mm256_maskload_pd(_shift6[id_i] + id_j0, memoryMask_first);
						const __m256d sh1 = _mm256_maskload_pd(_shift6[id_i] + id_j1, memoryMask_first);
						const __m256d sh2 = _mm256_maskload_pd(_shift6[id_i] + id_j2, memoryMask_first);
						const __m256d sh3 = _mm256_maskload_pd(_shift6[id_i] + id_j3, memoryMask_first);

						const __m256d sh0sh1 = _mm256_unpacklo_pd(sh0, sh1);
						const __m256d sh2sh3 = _mm256_unpacklo_pd(sh2, sh3);

						const __m256d shift6 = _mm256_permute2f128_pd(sh0sh1, sh2sh3, 1<<5);

						const __m256d upot = _mm256_mul_pd(eps_24, lj12m6);
						const __m256d upot_sh = _mm256_add_pd(shift6, upot);
						const __m256d upot_masked = _mm256_and_pd(upot_sh, macroMask);

						sum_upot = _mm256_add_pd(sum_upot, upot_masked);

						const __m256d vir_x = _mm256_mul_pd(m_dx, fx);
						const __m256d vir_y = _mm256_mul_pd(m_dy, fy);
						const __m256d vir_z = _mm256_mul_pd(m_dz, fz);

						const __m256d vir_xy = _mm256_add_pd(vir_x, vir_y);
						const __m256d virial = _mm256_add_pd(vir_xy, vir_z);
						const __m256d vir_masked = _mm256_and_pd(virial, macroMask);

						sum_virial = _mm256_add_pd(sum_virial, vir_masked);
					}
					const __m256d old_fx2 = _mm256_load_pd(p_cfx2 + j);
					const __m256d new_fx2 = _mm256_sub_pd(old_fx2, fx);
					_mm256_store_pd(p_cfx2 + j, new_fx2);
					const __m256d old_fy2 = _mm256_load_pd(p_cfy2 + j);
					const __m256d new_fy2 = _mm256_sub_pd(old_fy2, fy);
					_mm256_store_pd(p_cfy2 + j, new_fy2);
					const __m256d old_fz2 = _mm256_load_pd(p_cfz2 + j);
					const __m256d new_fz2 = _mm256_sub_pd(old_fz2, fz);
					_mm256_store_pd(p_cfz2 + j, new_fz2);
					sum_fx1 = _mm256_add_pd(sum_fx1, fx);
					sum_fy1 = _mm256_add_pd(sum_fy1, fy);
					sum_fz1 = _mm256_add_pd(sum_fz1, fz);
				}
			}

			const __m256d sum_fx1_t1 = _mm256_permute2f128_pd(sum_fx1, sum_fx1, 0x1);
			const __m256d sum_fx1_t2 = _mm256_hadd_pd(sum_fx1, sum_fx1_t1);
			const __m256d sum_fx1_t3 = _mm256_hadd_pd(sum_fx1_t2, sum_fx1_t2);

			_mm256_maskstore_pd(
					p_cfx1 + i_center_idx,
					memoryMask_first,
					_mm256_add_pd(
							sum_fx1_t3,
							_mm256_maskload_pd(p_cfx1 + i_center_idx, memoryMask_first)
					)
			);

			const __m256d sum_fy1_t1 = _mm256_permute2f128_pd(sum_fy1, sum_fy1, 0x1);
			const __m256d sum_fy1_t2 = _mm256_hadd_pd(sum_fy1, sum_fy1_t1);
			const __m256d sum_fy1_t3 = _mm256_hadd_pd(sum_fy1_t2, sum_fy1_t2);
			_mm256_maskstore_pd(
					p_cfy1 + i_center_idx,
					memoryMask_first,
					_mm256_add_pd(
							sum_fy1_t3,
							_mm256_maskload_pd(p_cfy1 + i_center_idx, memoryMask_first)
					)
			);

			const __m256d sum_fz1_t1 = _mm256_permute2f128_pd(sum_fz1, sum_fz1, 0x1);
			const __m256d sum_fz1_t2 = _mm256_hadd_pd(sum_fz1, sum_fz1_t1);
			const __m256d sum_fz1_t3 = _mm256_hadd_pd(sum_fz1_t2, sum_fz1_t2);
			_mm256_maskstore_pd(
					p_cfz1 + i_center_idx,
					memoryMask_first,
					_mm256_add_pd(
							sum_fz1_t3,
							_mm256_maskload_pd(p_cfz1 + i_center_idx, memoryMask_first)
					)
			);

			// Unvectorized calculation for leftover pairs.
			for (; j < soa2._num_ljcenters; ++j) {
				_loopBodyNovec<ForcePolicy, MacroPolicy>(soa1, i_center_idx, soa2, j, p_center_dist_lookup + j);
			}
			
			i_center_idx++;
		}
	}
	const __m256d sum_upot_t1 = _mm256_permute2f128_pd(sum_upot, sum_upot, 0x1);
	const __m256d sum_upot_t2 = _mm256_hadd_pd(sum_upot, sum_upot_t1);
	const __m256d sum_upot_t3 = _mm256_hadd_pd(sum_upot_t2, sum_upot_t2);
	_mm256_maskstore_pd(
			&_upot6lj,
			memoryMask_first,
			_mm256_add_pd(
				sum_upot_t3,
				_mm256_maskload_pd(&_upot6lj, memoryMask_first)
			)

	);

	const __m256d sum_virial_t1 = _mm256_permute2f128_pd(sum_virial, sum_virial, 0x1);
	const __m256d sum_virial_t2 = _mm256_hadd_pd(sum_virial, sum_virial_t1);
	const __m256d sum_virial_t3 = _mm256_hadd_pd(sum_virial_t2, sum_virial_t2);
	_mm256_maskstore_pd(
			&_virial,
			memoryMask_first,
			_mm256_add_pd(
				sum_virial_t3,
				_mm256_maskload_pd(&_virial, memoryMask_first)
			)

	);

#endif
} // void LennardJonesCellHandler::CalculatePairs_(LJSoA & soa1, LJSoA & soa2)

void VectorizedCellProcessor::processCell(ParticleCell & c) {
	assert(c.getCellDataSoA());
	if (c.isHaloCell() || (c.getCellDataSoA()->_num_ljcenters < 2))
		return;

	_calculatePairs<SingleCellPolicy_, AllMacroPolicy_>(*(c.getCellDataSoA()), *(c.getCellDataSoA()));
}

double VectorizedCellProcessor::processSingleMolecule(Molecule* m1, ParticleCell& cell2)
{
	exit(666);
	return 666;
}

void VectorizedCellProcessor::processCellPair(ParticleCell & c1,
		ParticleCell & c2) {
	assert(&c1 != &c2);
	assert(c1.getCellDataSoA());
	assert(c2.getCellDataSoA());

	if ((c1.getCellDataSoA()->_num_ljcenters == 0) || (c2.getCellDataSoA()->_num_ljcenters == 0)) {
		return;
	}

	if (!(c1.isHaloCell() || c2.isHaloCell())) {
		_calculatePairs<CellPairPolicy_, AllMacroPolicy_>(*(c1.getCellDataSoA()), *(c2.getCellDataSoA()));
	} else if (c1.isHaloCell() == (!c2.isHaloCell())) {
		_calculatePairs<CellPairPolicy_, SomeMacroPolicy_>(*(c1.getCellDataSoA()), *(c2.getCellDataSoA()));
	} else {
		return;
	}
}
