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

VectorizedCellProcessor::VectorizedCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		_domain(domain), _cutoffRadiusSquare(cutoffRadius * cutoffRadius), _LJcutoffRadiusSquare(LJcutoffRadius * LJcutoffRadius),
		// maybe move the following to somewhere else:
		_epsRFInvrc3(2. * (domain.getepsilonRF() - 1.) / ((cutoffRadius * cutoffRadius * cutoffRadius) * (2. * domain.getepsilonRF() + 1.))), 
		_compIDs(), _eps_sig(), _shift6(), _upot6lj(0.0), _upotXpoles(0.0), _virial(0.0), _myRF(0.0), _centers_dist_lookup(128) {

#if VCP_VEC_TYPE==VCP_NOVEC
	global_log->info() << "VectorizedCellProcessor: using no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	global_log->info() << "VectorizedCellProcessor: using SSE3 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	global_log->info() << "VectorizedCellProcessor: using AVX intrinsics." << std::endl;
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
	_upotXpoles = 0.0;
	_myRF = 0.0;

	global_log->debug() << "VectorizedLJCellProcessor::initTraversal() to " << numCells << " cells." << std::endl;

	if (numCells > _particleCellDataVector.size()) {
		for (size_t i = _particleCellDataVector.size(); i < numCells; i++) {
			_particleCellDataVector.push_back(new CellDataSoA(64,64,64,64,64));
		}
		global_log->debug() << "resize CellDataSoA to " << numCells << " cells." << std::endl;
	}

}


void VectorizedCellProcessor::endTraversal() {
	_domain.setLocalVirial(_virial + 3.0 * _myRF);
	_domain.setLocalUpot(_upot6lj / 6.0 + _upotXpoles + _myRF);
}


void VectorizedCellProcessor::preprocessCell(ParticleCell & c) {
	assert(!c.getCellDataSoA());

	const MoleculeList & molecules = c.getParticlePointers();

	// Determine the total number of centers.
	size_t numMolecules = molecules.size();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;
	
	for (size_t m = 0;  m < numMolecules; ++m) {
		nLJCenters += molecules[m]->numLJcenters();
		nCharges += molecules[m]->numCharges();
		nDipoles += molecules[m]->numDipoles();
		nQuadrupoles += molecules[m]->numQuadrupoles();
	}

	// Construct the SoA.
	assert(!_particleCellDataVector.empty()); 
	CellDataSoA* soaPtr = _particleCellDataVector.back();
//	global_log->debug() << " _particleCellDataVector.size()=" << _particleCellDataVector.size() << " soaPtr=" << soaPtr
//						<< " nCenters=" << nLJCenters + nCharges + nDipoles + nQuadrupoles << std::endl;
	CellDataSoA & soa = *soaPtr;
	soa.setDistLookup(_centers_dist_lookup);
	soa.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);
	c.setCellDataSoA(soaPtr);
	_particleCellDataVector.pop_back();

	ComponentList components = *(_simulation.getEnsemble()->components());

	size_t iLJCenters = 0;
	size_t iCharges = 0;
	size_t iDipoles = 0;
	size_t iQuadrupoles = 0;
	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < molecules.size(); ++i) {
		const size_t mol_ljc_num = molecules[i]->numLJcenters();
		const size_t mol_charges_num = molecules[i]->numCharges();
		const size_t mol_dipoles_num = molecules[i]->numDipoles();
		const size_t mol_quadrupoles_num = molecules[i]->numQuadrupoles();
		const double mol_pos_x = molecules[i]->r(0);
		const double mol_pos_y = molecules[i]->r(1);
		const double mol_pos_z = molecules[i]->r(2);

		soa._mol_pos_x[i] = mol_pos_x;
		soa._mol_pos_y[i] = mol_pos_y;
		soa._mol_pos_z[i] = mol_pos_z;
		soa._mol_ljc_num[i] = mol_ljc_num;
		soa._mol_charges_num[i] = mol_charges_num;
		soa._mol_dipoles_num[i] = mol_dipoles_num;
		soa._mol_quadrupoles_num[i] = mol_quadrupoles_num;

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
			soa._ljc_id[iLJCenters] = _compIDs[molecules[i]->componentid()] + j;
			soa._ljc_dist_lookup[iLJCenters] = 0.0;
		}

		for (size_t j = 0; j < mol_charges_num; ++j, ++iCharges)
		{
			soa._charges_m_r_x[iCharges] = mol_pos_x;
			soa._charges_m_r_y[iCharges] = mol_pos_y;
			soa._charges_m_r_z[iCharges] = mol_pos_z;
			soa._charges_r_x[iCharges] = molecules[i]->charge_d(j)[0] + mol_pos_x;
			soa._charges_r_y[iCharges] = molecules[i]->charge_d(j)[1] + mol_pos_y;
			soa._charges_r_z[iCharges] = molecules[i]->charge_d(j)[2] + mol_pos_z;
			soa._charges_f_x[iCharges] = 0.0;
			soa._charges_f_y[iCharges] = 0.0;
			soa._charges_f_z[iCharges] = 0.0;
			soa._charges_dist_lookup[iCharges] = 0.0;
			// Get the charge
			soa._charges_q[iCharges] = components[molecules[i]->componentid()].charge(j).q();
		}

		for (size_t j = 0; j < mol_dipoles_num; ++j, ++iDipoles)
		{
			soa._dipoles_m_r_x[iDipoles] = mol_pos_x;
			soa._dipoles_m_r_y[iDipoles] = mol_pos_y;
			soa._dipoles_m_r_z[iDipoles] = mol_pos_z;
			soa._dipoles_r_x[iDipoles] = molecules[i]->dipole_d(j)[0] + mol_pos_x;
			soa._dipoles_r_y[iDipoles] = molecules[i]->dipole_d(j)[1] + mol_pos_y;
			soa._dipoles_r_z[iDipoles] = molecules[i]->dipole_d(j)[2] + mol_pos_z;
			soa._dipoles_f_x[iDipoles] = 0.0;
			soa._dipoles_f_y[iDipoles] = 0.0;
			soa._dipoles_f_z[iDipoles] = 0.0;
			soa._dipoles_dist_lookup[iDipoles] = 0.0;
			// Get the dipole moment
			soa._dipoles_p[iDipoles] = components[molecules[i]->componentid()].dipole(j).absMy();
			soa._dipoles_e_x[iDipoles] = molecules[i]->dipole_e(j)[0];
			soa._dipoles_e_y[iDipoles] = molecules[i]->dipole_e(j)[1];
			soa._dipoles_e_z[iDipoles] = molecules[i]->dipole_e(j)[2];
			soa._dipoles_M_x[iDipoles] = 0.0;
			soa._dipoles_M_y[iDipoles] = 0.0;
			soa._dipoles_M_z[iDipoles] = 0.0;
		}

		for (size_t j = 0; j < mol_quadrupoles_num; ++j, ++iQuadrupoles)
		{
			soa._quadrupoles_m_r_x[iQuadrupoles] = mol_pos_x;
			soa._quadrupoles_m_r_y[iQuadrupoles] = mol_pos_y;
			soa._quadrupoles_m_r_z[iQuadrupoles] = mol_pos_z;
			soa._quadrupoles_r_x[iQuadrupoles] = molecules[i]->quadrupole_d(j)[0] + mol_pos_x;
			soa._quadrupoles_r_y[iQuadrupoles] = molecules[i]->quadrupole_d(j)[1] + mol_pos_y;
			soa._quadrupoles_r_z[iQuadrupoles] = molecules[i]->quadrupole_d(j)[2] + mol_pos_z;
			soa._quadrupoles_f_x[iQuadrupoles] = 0.0;
			soa._quadrupoles_f_y[iQuadrupoles] = 0.0;
			soa._quadrupoles_f_z[iQuadrupoles] = 0.0;
			soa._quadrupoles_dist_lookup[iQuadrupoles] = 0.0;
			// Get the quadrupole moment
			soa._quadrupoles_m[iQuadrupoles] = components[molecules[i]->componentid()].quadrupole(j).absQ();
			soa._quadrupoles_e_x[iQuadrupoles] = molecules[i]->quadrupole_e(j)[0];
			soa._quadrupoles_e_y[iQuadrupoles] = molecules[i]->quadrupole_e(j)[1];
			soa._quadrupoles_e_z[iQuadrupoles] = molecules[i]->quadrupole_e(j)[2];
			soa._quadrupoles_M_x[iQuadrupoles] = 0.0;
			soa._quadrupoles_M_y[iQuadrupoles] = 0.0;
			soa._quadrupoles_M_z[iQuadrupoles] = 0.0;
		}
	}
}


void VectorizedCellProcessor::postprocessCell(ParticleCell & c) {
	assert(c.getCellDataSoA());
	CellDataSoA& soa = *c.getCellDataSoA();

	MoleculeList & molecules = c.getParticlePointers();

	// For each molecule iterate over all its centers.
	size_t iLJCenters = 0;
	size_t iCharges = 0;
	size_t iDipoles = 0;
	size_t iQuadrupoles = 0;
	size_t numMols = molecules.size();
	for (size_t m = 0; m < numMols; ++m) {
		const size_t mol_ljc_num = molecules[m]->numLJcenters();
		const size_t mol_charges_num = molecules[m]->numCharges();
		const size_t mol_dipoles_num = molecules[m]->numDipoles();
		const size_t mol_quadrupoles_num = molecules[m]->numQuadrupoles();

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
		}

		for (size_t i = 0; i < mol_charges_num; ++i, ++iCharges) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa._charges_f_x[iCharges];
			f[1] = soa._charges_f_y[iCharges];
			f[2] = soa._charges_f_z[iCharges];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			molecules[m]->Fchargeadd(i, f);
		}

		for (size_t i = 0; i < mol_dipoles_num; ++i, ++iDipoles) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa._dipoles_f_x[iDipoles];
			f[1] = soa._dipoles_f_y[iDipoles];
			f[2] = soa._dipoles_f_z[iDipoles];
			double M[3];
			M[0] = soa._dipoles_M_x[iDipoles];
			M[1] = soa._dipoles_M_y[iDipoles];
			M[2] = soa._dipoles_M_z[iDipoles];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			assert(!isnan(M[0]));
			assert(!isnan(M[1]));
			assert(!isnan(M[2]));
			molecules[m]->Fdipoleadd(i, f);
			molecules[m]->Madd(M);
		}

		for (size_t i = 0; i < mol_quadrupoles_num; ++i, ++iQuadrupoles) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa._quadrupoles_f_x[iQuadrupoles];
			f[1] = soa._quadrupoles_f_y[iQuadrupoles];
			f[2] = soa._quadrupoles_f_z[iQuadrupoles];
			double M[3];
			M[0] = soa._quadrupoles_M_x[iQuadrupoles];
			M[1] = soa._quadrupoles_M_y[iQuadrupoles];
			M[2] = soa._quadrupoles_M_z[iQuadrupoles];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			assert(!isnan(M[0]));
			assert(!isnan(M[1]));
			assert(!isnan(M[2]));
			molecules[m]->Fquadrupoleadd(i, f);
			molecules[m]->Madd(M);
		}
	}
	// Delete the SoA.
	_particleCellDataVector.push_back(&soa);
	c.setCellDataSoA(0);
}


template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecLJ (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
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

		const double m_dx = soa1._ljc_m_r_x[i] - soa2._ljc_m_r_x[j];
		const double m_dy = soa1._ljc_m_r_y[i] - soa2._ljc_m_r_y[j];
		const double m_dz = soa1._ljc_m_r_z[i] - soa2._ljc_m_r_z[j];

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
}  /* end of method VectorizedCellProcessor :: _loopBodyNovec */


template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecCharges (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
{
	if (*forceMask) {
		const double q1q2per4pie0 = soa1._charges_q[i] * soa2._charges_q[j];

		const double c_dx = soa1._charges_r_x[i] - soa2._charges_r_x[j];
		const double c_dy = soa1._charges_r_y[i] - soa2._charges_r_y[j];
		const double c_dz = soa1._charges_r_z[i] - soa2._charges_r_z[j];

		const double c_dr2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;
		const double c_dr2_inv = 1.0 / c_dr2;
		const double c_dr_inv = sqrt(c_dr2_inv);

		const double upot = q1q2per4pie0 * c_dr_inv;
		const double fac = upot * c_dr2_inv;

		const double f_x = c_dx * fac;
		const double f_y = c_dy * fac;
		const double f_z = c_dz * fac;

		const double m_dx = soa1._charges_m_r_x[i] - soa2._charges_m_r_x[j];
		const double m_dy = soa1._charges_m_r_y[i] - soa2._charges_m_r_y[j];
		const double m_dz = soa1._charges_m_r_z[i] - soa2._charges_m_r_z[j];

		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueCondition(m_dx, m_dy, m_dz)) {
			_upotXpoles += upot;
			_virial += m_dx * f_x + m_dy * f_y + m_dz * f_z;
		}

		// Add the force to center 1, and subtract it from center 2.
		soa1._charges_f_x[i] += f_x;
		soa1._charges_f_y[i] += f_y;
		soa1._charges_f_z[i] += f_z;

		soa2._charges_f_x[j] -= f_x;
		soa2._charges_f_y[j] -= f_y;
		soa2._charges_f_z[j] -= f_z;
	}
}

template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecChargesDipoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask, const bool& switched)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {
		const double e_x = soa2._dipoles_e_x[j];
		const double e_y = soa2._dipoles_e_y[j];
		const double e_z = soa2._dipoles_e_z[j];

		const double c_dx = soa1._charges_r_x[i] - soa2._dipoles_r_x[j];
		const double c_dy = soa1._charges_r_y[i] - soa2._dipoles_r_y[j];
		const double c_dz = soa1._charges_r_z[i] - soa2._dipoles_r_z[j];

		const double c_dr2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;
		const double c_dr2_inv = 1.0 / c_dr2;
		const double c_dr_inv = sqrt(c_dr2_inv);

		const double re = c_dx * e_x + c_dy * e_y + c_dz * e_z;
		const double qpper4pie0 = soa1._charges_q[i] * soa2._dipoles_p[j];

		const double fac = 3.0 * c_dr2_inv * re;
		const double qpper4pie0dr3 = qpper4pie0 * c_dr2_inv * c_dr_inv;

		// TODO: SHOULD BE NEGATED!!! In potforce.h (non-vectorized version) there is also a potential sign error.
		// If so, both implementations have to be fixed. For now, they are consistent. (Uwe Ehmann)
		const double f_x = qpper4pie0dr3 * (e_x - c_dx * fac);
		const double f_y = qpper4pie0dr3 * (e_y - c_dy * fac );
		const double f_z = qpper4pie0dr3 * (e_z - c_dz * fac);

		const double m_dx = soa1._charges_m_r_x[i] - soa2._dipoles_m_r_x[j];
		const double m_dy = soa1._charges_m_r_y[i] - soa2._dipoles_m_r_y[j];
		const double m_dz = soa1._charges_m_r_z[i] - soa2._dipoles_m_r_z[j];

		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueConditionSwitched(m_dx, m_dy, m_dz, switched)) {
			const double minusUpot = qpper4pie0dr3 * re;
			_upotXpoles -= minusUpot;
			_virial += m_dx * f_x + m_dy * f_y + m_dz * f_z;
		}

		const double M_x = qpper4pie0dr3 * (e_y * c_dz - e_z * c_dy);
		const double M_y = qpper4pie0dr3 * (e_z * c_dx - e_x * c_dz);
		const double M_z = qpper4pie0dr3 * (e_x * c_dy - e_y * c_dx);

		soa1._charges_f_x[i] += f_x;
		soa1._charges_f_y[i] += f_y;
		soa1._charges_f_z[i] += f_z;

		soa2._dipoles_f_x[j] -= f_x;
		soa2._dipoles_f_y[j] -= f_y;
		soa2._dipoles_f_z[j] -= f_z;

		soa2._dipoles_M_x[j] += M_x;
		soa2._dipoles_M_y[j] += M_y;
		soa2._dipoles_M_z[j] += M_z;
	}
}

template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecDipoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {
		const double e1_x = soa1._dipoles_e_x[i];
		const double e1_y = soa1._dipoles_e_y[i];
		const double e1_z = soa1._dipoles_e_z[i];

		const double e2_x = soa2._dipoles_e_x[j];
		const double e2_y = soa2._dipoles_e_y[j];
		const double e2_z = soa2._dipoles_e_z[j];

		const double c_dx = soa1._dipoles_r_x[i] - soa2._dipoles_r_x[j];
		const double c_dy = soa1._dipoles_r_y[i] - soa2._dipoles_r_y[j];
		const double c_dz = soa1._dipoles_r_z[i] - soa2._dipoles_r_z[j];

		const double c_dr2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;
		const double c_dr2_inv = 1. / c_dr2;
		const double c_dr_inv = sqrt(c_dr2_inv);
		const double c_dr2three_inv = 3. * c_dr2_inv;

		const double p1p2 = soa1._dipoles_p[i] * soa2._dipoles_p[j];
		const double p1p2per4pie0 = p1p2;
		const double rffac = p1p2 * _epsRFInvrc3;

		const double p1p2per4pie0r3 = p1p2per4pie0 * c_dr_inv * c_dr2_inv;
		const double p1p2threeper4pie0r5 = p1p2per4pie0r3 * c_dr2three_inv;

		const double e1e2 = e1_x * e2_x + e1_y * e2_y + e1_z * e2_z;
		const double re1 = c_dx * e1_x + c_dy * e1_y + c_dz * e1_z;
		const double re2 = c_dx * e2_x + c_dy * e2_y + c_dz * e2_z;

		const double re1threeperr2 = re1 * c_dr2three_inv;
		const double re2threeperr2 = re2 * c_dr2three_inv;
		const double re1re2perr2 = re1 * re2 * c_dr2_inv;

		const double e1e2minus5re1re2perr2 = e1e2 - 5. * re1re2perr2;

		const double f_x = p1p2threeper4pie0r5 * ( c_dx * e1e2minus5re1re2perr2 + e1_x * re2 + e2_x * re1 );
		const double f_y = p1p2threeper4pie0r5 * ( c_dy * e1e2minus5re1re2perr2 + e1_y * re2 + e2_y * re1 );
		const double f_z = p1p2threeper4pie0r5 * ( c_dz * e1e2minus5re1re2perr2 + e1_z * re2 + e2_z * re1 );

		const double m_dx = soa1._dipoles_m_r_x[i] - soa2._dipoles_m_r_x[j];
		const double m_dy = soa1._dipoles_m_r_y[i] - soa2._dipoles_m_r_y[j];
		const double m_dz = soa1._dipoles_m_r_z[i] - soa2._dipoles_m_r_z[j];

		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueCondition(m_dx, m_dy, m_dz)) {
			const double upot = p1p2per4pie0r3 * (e1e2 - 3.0 * re1re2perr2);
			_upotXpoles += upot;
			_virial += m_dx * f_x + m_dy * f_y + m_dz * f_z;
			_myRF -= rffac * e1e2;
		}

		const double e1_x_e2_y = e1_x * e2_y;
		const double e1_x_e2_z = e1_x * e2_z;
		const double e1_y_e2_x = e1_y * e2_x;
		const double e1_y_e2_z = e1_y * e2_z;
		const double e1_z_e2_x = e1_z * e2_x;
		const double e1_z_e2_y = e1_z * e2_y;

		const double e1_x_e2_y_minus_e1_y_e2_x = e1_x_e2_y - e1_y_e2_x;
		const double e1_y_e2_z_minus_e1_z_e2_y = e1_y_e2_z - e1_z_e2_y;
		const double e1_z_e2_x_minus_e1_x_e2_z = e1_z_e2_x - e1_x_e2_z;

		const double M1_x = p1p2per4pie0r3 * (re2threeperr2 * (e1_y * c_dz - e1_z * c_dy) - e1_y_e2_z_minus_e1_z_e2_y) + rffac * e1_y_e2_z_minus_e1_z_e2_y;
		const double M1_y = p1p2per4pie0r3 * (re2threeperr2 * (e1_z * c_dx - e1_x * c_dz) - e1_z_e2_x_minus_e1_x_e2_z) + rffac * e1_z_e2_x_minus_e1_x_e2_z;
		const double M1_z = p1p2per4pie0r3 * (re2threeperr2 * (e1_x * c_dy - e1_y * c_dx) - e1_x_e2_y_minus_e1_y_e2_x) + rffac * e1_x_e2_y_minus_e1_y_e2_x;

		const double M2_x = p1p2per4pie0r3 * (re1threeperr2 * (e2_y * c_dz - e2_z * c_dy) + e1_y_e2_z_minus_e1_z_e2_y) - rffac * e1_y_e2_z_minus_e1_z_e2_y;
		const double M2_y = p1p2per4pie0r3 * (re1threeperr2 * (e2_z * c_dx - e2_x * c_dz) + e1_z_e2_x_minus_e1_x_e2_z) - rffac * e1_z_e2_x_minus_e1_x_e2_z;
		const double M2_z = p1p2per4pie0r3 * (re1threeperr2 * (e2_x * c_dy - e2_y * c_dx) + e1_x_e2_y_minus_e1_y_e2_x) - rffac * e1_x_e2_y_minus_e1_y_e2_x;

		soa1._dipoles_f_x[i] += f_x;
		soa1._dipoles_f_y[i] += f_y;
		soa1._dipoles_f_z[i] += f_z;

		soa2._dipoles_f_x[j] -= f_x;
		soa2._dipoles_f_y[j] -= f_y;
		soa2._dipoles_f_z[j] -= f_z;

		soa1._dipoles_M_x[i] += M1_x;
		soa1._dipoles_M_y[i] += M1_y;
		soa1._dipoles_M_z[i] += M1_z;

		soa2._dipoles_M_x[j] += M2_x;
		soa2._dipoles_M_y[j] += M2_y;
		soa2._dipoles_M_z[j] += M2_z;
	}
}

template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecChargesQuadrupoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask, const bool& switched)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {
		// orientation vector of quadrupolar moment
		const double ejj_x = soa2._quadrupoles_e_x[j];
		const double ejj_y = soa2._quadrupoles_e_y[j];
		const double ejj_z = soa2._quadrupoles_e_z[j];

		// Cartesian distance
		const double c_dx = soa1._charges_r_x[i] - soa2._quadrupoles_r_x[j];
		const double c_dy = soa1._charges_r_y[i] - soa2._quadrupoles_r_y[j];
		const double c_dz = soa1._charges_r_z[i] - soa2._quadrupoles_r_z[j];
		const double c_dr2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;

		const double invdr2 = 1.0 / c_dr2;
		const double invdr = sqrt(invdr2);

		const double qQ05per4pie0 = 0.5 * soa1._charges_q[i]*soa2._quadrupoles_m[j];

		double costj = ejj_x * c_dx + ejj_y * c_dy + ejj_z * c_dz;
		costj *= invdr;

		const double qQinv4dr3 = qQ05per4pie0 * invdr * invdr2;
		const double upot = qQinv4dr3 * (3.0 * costj * costj - 1);

		const double partialRijInvdr = -3.0 * upot * invdr2;
		const double partialTjInvdr = 6.0 * costj * qQinv4dr3 * invdr;
		const double fac = costj * partialTjInvdr * invdr - partialRijInvdr;

		const double f_x = fac * c_dx - partialTjInvdr * ejj_x;
		const double f_y = fac * c_dy - partialTjInvdr * ejj_y;
		const double f_z = fac * c_dz - partialTjInvdr * ejj_z;

		const double m_dx = soa1._charges_m_r_x[i] - soa2._quadrupoles_m_r_x[j];
		const double m_dy = soa1._charges_m_r_y[i] - soa2._quadrupoles_m_r_y[j];
		const double m_dz = soa1._charges_m_r_z[i] - soa2._quadrupoles_m_r_z[j];
		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueConditionSwitched(m_dx, m_dy, m_dz, switched)) {
			_upotXpoles += upot;
			_virial += m_dx * f_x + m_dy * f_y + m_dz * f_z;
		}

		const double minuseXrij_x = ejj_z * c_dy - ejj_y * c_dz;
		const double minuseXrij_y = ejj_x * c_dz - ejj_z * c_dx;
		const double minuseXrij_z = ejj_y * c_dx - ejj_x * c_dy;

		const double M2_x = partialTjInvdr * minuseXrij_x;
		const double M2_y = partialTjInvdr * minuseXrij_y;
		const double M2_z = partialTjInvdr * minuseXrij_z;

		soa1._charges_f_x[i] += f_x;
		soa1._charges_f_y[i] += f_y;
		soa1._charges_f_z[i] += f_z;

		soa2._quadrupoles_f_x[j] -= f_x;
		soa2._quadrupoles_f_y[j] -= f_y;
		soa2._quadrupoles_f_z[j] -= f_z;

		soa2._quadrupoles_M_x[j] += M2_x;
		soa2._quadrupoles_M_y[j] += M2_y;
		soa2._quadrupoles_M_z[j] += M2_z;
	}
}

template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecDipolesQuadrupoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask, const bool& switched)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {

		// orientation vector of dipole moment
		const double eii_x = soa1._dipoles_e_x[i];
		const double eii_y = soa1._dipoles_e_y[i];
		const double eii_z = soa1._dipoles_e_z[i];

		// orientation vector of quadrupolar moment
		const double ejj_x = soa2._quadrupoles_e_x[j];
		const double ejj_y = soa2._quadrupoles_e_y[j];
		const double ejj_z = soa2._quadrupoles_e_z[j];

		// Cartesian distance
		const double c_dx = soa1._dipoles_r_x[i] - soa2._quadrupoles_r_x[j];
		const double c_dy = soa1._dipoles_r_y[i] - soa2._quadrupoles_r_y[j];
		const double c_dz = soa1._dipoles_r_z[i] - soa2._quadrupoles_r_z[j];
		const double c_dr2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;

		const double invdr2 = 1.0 / c_dr2;
		const double invdr = sqrt(invdr2);

		const double myqfac = 1.5 * soa1._dipoles_p[i] * soa2._quadrupoles_m[j] * invdr2 * invdr2;

		double costi = eii_x * c_dx + eii_y * c_dy + eii_z * c_dz;
		double costj = ejj_x * c_dx + ejj_y * c_dy + ejj_z * c_dz;
		costi *= invdr;
		costj *= invdr;
		const double cos2tj = costj * costj;
		const double cosgij = eii_x * ejj_x + eii_y * ejj_y + eii_z * ejj_z;

		//Potential
		// TODO: Check if upot has to be multiplied by -1 according to DISS_STOLL S.178.
		// This affects also the implementation in potforce.h
		const double upot = myqfac * (-costi * (5. * cos2tj - 1.) + 2. * cosgij * costj);

		const double partialRijInvdr1 = -4. * upot * invdr2;
		const double partialTiInvdr1 = myqfac * (-5. * cos2tj + 1.) * invdr;
		const double partialTjInvdr1 = myqfac * 2. * (-5. * costi * costj + cosgij) * invdr;
		const double partialGij = myqfac * 2. * costj;
		const double fac = -partialRijInvdr1 + (costi*partialTiInvdr1 + costj*partialTjInvdr1) * invdr;

		const double f_x = fac * c_dx - partialTiInvdr1 * eii_x - partialTjInvdr1 * ejj_x;
		const double f_y = fac * c_dy - partialTiInvdr1 * eii_y - partialTjInvdr1 * ejj_y;
		const double f_z = fac * c_dz - partialTiInvdr1 * eii_z - partialTjInvdr1 * ejj_z;

		const double m_dx = soa1._dipoles_m_r_x[i] - soa2._quadrupoles_m_r_x[j];
		const double m_dy = soa1._dipoles_m_r_y[i] - soa2._quadrupoles_m_r_y[j];
		const double m_dz = soa1._dipoles_m_r_z[i] - soa2._quadrupoles_m_r_z[j];

		if (MacroPolicy :: MacroscopicValueConditionSwitched(m_dx, m_dy, m_dz, switched)) {
			_upotXpoles += upot;
			_virial += m_dx * f_x + m_dy * f_y + m_dz * f_z;
		}

		double const partialGij_eiXej_x = partialGij * (eii_y * ejj_z - eii_z * ejj_y);
		double const partialGij_eiXej_y = partialGij * (eii_z * ejj_x - eii_x * ejj_z);
		double const partialGij_eiXej_z = partialGij * (eii_x * ejj_y - eii_y * ejj_x);

		double eXrij_x = eii_y * c_dz - eii_z * c_dy;
		double eXrij_y = eii_z * c_dx - eii_x * c_dz;
		double eXrij_z = eii_x * c_dy - eii_y * c_dx;

		double const M1_x = -partialTiInvdr1 * eXrij_x - partialGij_eiXej_x;
		double const M1_y = -partialTiInvdr1 * eXrij_y - partialGij_eiXej_y;
		double const M1_z = -partialTiInvdr1 * eXrij_z - partialGij_eiXej_z;

		eXrij_x = ejj_y * c_dz - ejj_z * c_dy;
		eXrij_y = ejj_z * c_dx - ejj_x * c_dz;
		eXrij_z = ejj_x * c_dy - ejj_y * c_dx;

		double const M2_x = -partialTjInvdr1 * eXrij_x + partialGij_eiXej_x;
		double const M2_y = -partialTjInvdr1 * eXrij_y + partialGij_eiXej_y;
		double const M2_z = -partialTjInvdr1 * eXrij_z + partialGij_eiXej_z;

		soa1._dipoles_f_x[i] += f_x;
		soa1._dipoles_f_y[i] += f_y;
		soa1._dipoles_f_z[i] += f_z;

		soa2._quadrupoles_f_x[j] -= f_x;
		soa2._quadrupoles_f_y[j] -= f_y;
		soa2._quadrupoles_f_z[j] -= f_z;

		soa1._dipoles_M_x[i] += M1_x;
		soa1._dipoles_M_y[i] += M1_y;
		soa1._dipoles_M_z[i] += M1_z;

		soa2._quadrupoles_M_x[j] += M2_x;
		soa2._quadrupoles_M_y[j] += M2_y;
		soa2._quadrupoles_M_z[j] += M2_z;
	}
}

template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecQuadrupoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {

		// orientation vector of quadrupolar moment Qi
		const double eii_x = soa1._quadrupoles_e_x[i];
		const double eii_y = soa1._quadrupoles_e_y[i];
		const double eii_z = soa1._quadrupoles_e_z[i];

		// orientation vector of quadrupolar moment Qj
		const double ejj_x = soa2._quadrupoles_e_x[j];
		const double ejj_y = soa2._quadrupoles_e_y[j];
		const double ejj_z = soa2._quadrupoles_e_z[j];

		// Cartesian distance
		const double c_dx = soa1._quadrupoles_r_x[i] - soa2._quadrupoles_r_x[j];
		const double c_dy = soa1._quadrupoles_r_y[i] - soa2._quadrupoles_r_y[j];
		const double c_dz = soa1._quadrupoles_r_z[i] - soa2._quadrupoles_r_z[j];
		const double c_dr2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;

		const double invdr2 = 1.0 / c_dr2;
		const double invdr = sqrt(invdr2);

		const double qfac = 0.75 * soa1._quadrupoles_m[i] * soa2._quadrupoles_m[j] * invdr2 * invdr2 * invdr;

		double costi = eii_x * c_dx + eii_y * c_dy + eii_z * c_dz;
		double costj = ejj_x * c_dx + ejj_y * c_dy + ejj_z * c_dz;
		costi *= invdr;
		costj *= invdr;
		const double cos2ti = costi * costi;
		const double cos2tj = costj * costj;
		const double cosgij = eii_x * ejj_x + eii_y * ejj_y + eii_z * ejj_z;

		const double term = cosgij - 5. * costi * costj;

		//Potential
		const double upot = qfac * (1. - 5. * (cos2ti + cos2tj) -
				15. * cos2ti * cos2tj + 2. * term * term);

		const double partialRijInvdr1 = -5. * upot * invdr2;
		const double partialTiInvdr1 = -qfac * 10. * (costi + 3. * costi * cos2tj + 2. * costj * term) * invdr;
		const double partialTjInvdr1 = -qfac * 10. * (costj + 3. * cos2ti * costj + 2. * costi * term) * invdr;
		const double partialGij = qfac * 4. * term;
		const double fac = -partialRijInvdr1 + (costi*partialTiInvdr1 + costj*partialTjInvdr1) * invdr;

		const double f_x = fac * c_dx - partialTiInvdr1 * eii_x - partialTjInvdr1 * ejj_x;
		const double f_y = fac * c_dy - partialTiInvdr1 * eii_y - partialTjInvdr1 * ejj_y;
		const double f_z = fac * c_dz - partialTiInvdr1 * eii_z - partialTjInvdr1 * ejj_z;

		const double m_dx = soa1._quadrupoles_m_r_x[i] - soa2._quadrupoles_m_r_x[j];
		const double m_dy = soa1._quadrupoles_m_r_y[i] - soa2._quadrupoles_m_r_y[j];
		const double m_dz = soa1._quadrupoles_m_r_z[i] - soa2._quadrupoles_m_r_z[j];

		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueCondition(m_dx, m_dy, m_dz)) {
			_upotXpoles += upot;
			_virial += m_dx * f_x + m_dy * f_y + m_dz * f_z;
		}

		double const partialGij_eiXej_x = partialGij * (eii_y * ejj_z - eii_z * ejj_y);
		double const partialGij_eiXej_y = partialGij * (eii_z * ejj_x - eii_x * ejj_z);
		double const partialGij_eiXej_z = partialGij * (eii_x * ejj_y - eii_y * ejj_x);

		double eXrij_x = eii_y * c_dz - eii_z * c_dy;
		double eXrij_y = eii_z * c_dx - eii_x * c_dz;
		double eXrij_z = eii_x * c_dy - eii_y * c_dx;

		double const M1_x = -partialTiInvdr1 * eXrij_x - partialGij_eiXej_x;
		double const M1_y = -partialTiInvdr1 * eXrij_y - partialGij_eiXej_y;
		double const M1_z = -partialTiInvdr1 * eXrij_z - partialGij_eiXej_z;

		eXrij_x = ejj_y * c_dz - ejj_z * c_dy;
		eXrij_y = ejj_z * c_dx - ejj_x * c_dz;
		eXrij_z = ejj_x * c_dy - ejj_y * c_dx;

		double const M2_x = -partialTjInvdr1 * eXrij_x + partialGij_eiXej_x;
		double const M2_y = -partialTjInvdr1 * eXrij_y + partialGij_eiXej_y;
		double const M2_z = -partialTjInvdr1 * eXrij_z + partialGij_eiXej_z;

		soa1._quadrupoles_f_x[i] += f_x;
		soa1._quadrupoles_f_y[i] += f_y;
		soa1._quadrupoles_f_z[i] += f_z;

		soa2._quadrupoles_f_x[j] -= f_x;
		soa2._quadrupoles_f_y[j] -= f_y;
		soa2._quadrupoles_f_z[j] -= f_z;

		soa1._quadrupoles_M_x[i] += M1_x;
		soa1._quadrupoles_M_y[i] += M1_y;
		soa1._quadrupoles_M_z[i] += M1_z;

		soa2._quadrupoles_M_x[j] += M2_x;
		soa2._quadrupoles_M_y[j] += M2_y;
		soa2._quadrupoles_M_z[j] += M2_z;
	}
}


#if VCP_VEC_TYPE==VCP_VEC_SSE3

const vcp_double_vec zero = vcp_simd_zerov();
const vcp_double_vec one = _mm_set1_pd(1.0);
const vcp_double_vec two = _mm_set1_pd(2.0);
const vcp_double_vec three = _mm_set1_pd(3.0);
const vcp_double_vec four = _mm_set1_pd(4.0);
const vcp_double_vec five = _mm_set1_pd(5.0);
const vcp_double_vec six = _mm_set1_pd(6.0);
const vcp_double_vec ten = _mm_set1_pd(10.0);
const vcp_double_vec _05 = _mm_set1_pd(0.5);
const vcp_double_vec _075 = _mm_set1_pd(0.75);
const vcp_double_vec _1pt5 = _mm_set1_pd(1.5);
const vcp_double_vec _15 = _mm_set1_pd(15.0);

template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyLJ(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& sum_upot6lj, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask,
		const vcp_double_vec& e1s1, const vcp_double_vec& e2s2,
		const size_t& id_j0, const size_t& id_j1, const size_t& id_i)
{
	const vcp_double_vec c_dx = _mm_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dxdx = _mm_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dydy = _mm_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dzdz = _mm_mul_pd(c_dz, c_dz);
	const vcp_double_vec c_dxdx_dydy = _mm_add_pd(c_dxdx, c_dydy);
	const vcp_double_vec c_r2 = _mm_add_pd(c_dxdx_dydy, c_dzdz);
	const vcp_double_vec r2_inv_unmasked = _mm_div_pd(one, c_r2);
	const vcp_double_vec r2_inv = vcp_simd_and(r2_inv_unmasked, forceMask);

	const vcp_double_vec eps_24 = _mm_unpacklo_pd(e1s1, e2s2);
	const vcp_double_vec sig2 = _mm_unpackhi_pd(e1s1, e2s2);
	const vcp_double_vec lj2 = _mm_mul_pd(sig2, r2_inv);
	const vcp_double_vec lj4 = _mm_mul_pd(lj2, lj2);
	const vcp_double_vec lj6 = _mm_mul_pd(lj4, lj2);
	const vcp_double_vec lj12 = _mm_mul_pd(lj6, lj6);
	const vcp_double_vec lj12m6 = _mm_sub_pd(lj12, lj6);
	const vcp_double_vec eps24r2inv = _mm_mul_pd(eps_24, r2_inv);
	const vcp_double_vec lj12lj12m6 = _mm_add_pd(lj12, lj12m6);
	const vcp_double_vec scale = _mm_mul_pd(eps24r2inv, lj12lj12m6);
	f_x = _mm_mul_pd(c_dx, scale);
	f_y = _mm_mul_pd(c_dy, scale);
	f_z = _mm_mul_pd(c_dz, scale);

	const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
	// Only go on if at least 1 macroscopic value has to be calculated.
	if (_mm_movemask_pd(macroMask) > 0) {
		const vcp_double_vec sh1 = _mm_load_sd(_shift6[id_i] + id_j0);
		const vcp_double_vec sh2 = _mm_load_sd(_shift6[id_i] + id_j1);
		const vcp_double_vec shift6 = _mm_unpacklo_pd(sh1, sh2);
		const vcp_double_vec upot = _mm_mul_pd(eps_24, lj12m6);
		const vcp_double_vec upot_sh = _mm_add_pd(shift6, upot);
		const vcp_double_vec upot_masked = vcp_simd_and(upot_sh, macroMask);
		sum_upot6lj = _mm_add_pd(sum_upot6lj, upot_masked);
		const vcp_double_vec vir_x = _mm_mul_pd(m_dx, f_x);
		const vcp_double_vec vir_y = _mm_mul_pd(m_dy, f_y);
		const vcp_double_vec vir_z = _mm_mul_pd(m_dz, f_z);
		const vcp_double_vec vir_xy = _mm_add_pd(vir_x, vir_y);
		const vcp_double_vec virial = _mm_add_pd(vir_xy, vir_z);
		const vcp_double_vec vir_masked = vcp_simd_and(virial, macroMask);
		sum_virial = _mm_add_pd(sum_virial, vir_masked);
	}
}

template<class MacroPolicy>
inline void _loopBodyCharge(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& qii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& qjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask)
{
	const vcp_double_vec c_dx = _mm_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm_add_pd(_mm_add_pd(c_dx2, c_dy2), c_dz2);
	const vcp_double_vec c_dr2_inv_unmasked = _mm_div_pd(one, c_dr2);
	const vcp_double_vec c_dr2_inv = vcp_simd_and(c_dr2_inv_unmasked, forceMask);
    const vcp_double_vec c_dr_inv = _mm_sqrt_pd(c_dr2_inv);

	const vcp_double_vec q1q2per4pie0 = _mm_mul_pd(qii, qjj);
	const vcp_double_vec upot = _mm_mul_pd(q1q2per4pie0, c_dr_inv);
	const vcp_double_vec fac = _mm_mul_pd(upot, c_dr2_inv);

	f_x = _mm_mul_pd(c_dx, fac);
	f_y = _mm_mul_pd(c_dy, fac);
	f_z = _mm_mul_pd(c_dz, fac);

	const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm_movemask_pd(macroMask) > 0) {
		const vcp_double_vec upot_masked = vcp_simd_and(upot, macroMask);
		sum_upotXpoles = _mm_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm_add_pd(_mm_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = vcp_simd_and(virial, macroMask);
		sum_virial = _mm_add_pd(sum_virial, virial_masked);

	}
}

template<class MacroPolicy>
inline void _loopBodyChargeDipole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& q,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& e_x, const vcp_double_vec& e_y, const vcp_double_vec& e_z,
		const vcp_double_vec& p,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched)
{
	const vcp_double_vec dx = _mm_sub_pd(r1_x, r2_x);
	const vcp_double_vec dy = _mm_sub_pd(r1_y, r2_y);
	const vcp_double_vec dz = _mm_sub_pd(r1_z, r2_z);

	const vcp_double_vec dx2 = _mm_mul_pd(dx, dx);
	const vcp_double_vec dy2 = _mm_mul_pd(dy, dy);
	const vcp_double_vec dz2 = _mm_mul_pd(dz, dz);

	const vcp_double_vec dr2 = _mm_add_pd(_mm_add_pd(dx2, dy2), dz2);

	const vcp_double_vec dr2_inv_unmasked = _mm_div_pd(one, dr2);
	const vcp_double_vec dr2_inv = vcp_simd_and(dr2_inv_unmasked, forceMask);
	const vcp_double_vec dr_inv = _mm_sqrt_pd(dr2_inv);
	const vcp_double_vec dr3_inv = _mm_mul_pd(dr2_inv, dr_inv);

	const vcp_double_vec re = _mm_add_pd(_mm_mul_pd(dx, e_x), _mm_add_pd(_mm_mul_pd(dy, e_y), _mm_mul_pd(dz, e_z)));

	const vcp_double_vec qpper4pie0 = _mm_mul_pd(q, p);
	const vcp_double_vec qpper4pie0dr3 = _mm_mul_pd(qpper4pie0, dr3_inv);

	const vcp_double_vec fac = _mm_mul_pd(dr2_inv, _mm_mul_pd(three, re));

	f_x = _mm_mul_pd(qpper4pie0dr3, _mm_sub_pd(e_x, _mm_mul_pd(dx, fac)));
	f_y = _mm_mul_pd(qpper4pie0dr3, _mm_sub_pd(e_y, _mm_mul_pd(dy, fac)));
	f_z = _mm_mul_pd(qpper4pie0dr3, _mm_sub_pd(e_z, _mm_mul_pd(dz, fac)));

	const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm_movemask_pd(macroMask) > 0)
	{
		const vcp_double_vec minusUpot_unmasked =  _mm_mul_pd(qpper4pie0dr3, re);
		const vcp_double_vec minusUpot = vcp_simd_and(minusUpot_unmasked, macroMask);
		sum_upotXpoles = _mm_sub_pd(sum_upotXpoles, minusUpot);

		const vcp_double_vec virial_x = _mm_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm_mul_pd(m_dz, f_z);

		const vcp_double_vec virial_unmasked = _mm_add_pd(_mm_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial = vcp_simd_and(virial_unmasked, macroMask);
		sum_virial = _mm_add_pd(sum_virial, virial);
	}

	const vcp_double_vec e_x_dy = _mm_mul_pd(e_x, dy);
	const vcp_double_vec e_x_dz = _mm_mul_pd(e_x, dz);
	const vcp_double_vec e_y_dx = _mm_mul_pd(e_y, dx);
	const vcp_double_vec e_y_dz = _mm_mul_pd(e_y, dz);
	const vcp_double_vec e_z_dx = _mm_mul_pd(e_z, dx);
	const vcp_double_vec e_z_dy = _mm_mul_pd(e_z, dy);

	const vcp_double_vec e_x_dy_minus_e_y_dx = _mm_sub_pd(e_x_dy, e_y_dx);
	const vcp_double_vec e_y_dz_minus_e_z_dy = _mm_sub_pd(e_y_dz, e_z_dy);
	const vcp_double_vec e_z_dx_minus_e_x_dz = _mm_sub_pd(e_z_dx, e_x_dz);

	M_x = _mm_mul_pd(qpper4pie0dr3, e_y_dz_minus_e_z_dy);
	M_y = _mm_mul_pd(qpper4pie0dr3, e_z_dx_minus_e_x_dz);
	M_z = _mm_mul_pd(qpper4pie0dr3, e_x_dy_minus_e_y_dx);
}

template<class MacroPolicy>
inline void _loopBodyDipole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& pii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& pjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial, vcp_double_vec& sum_myRF,
		const vcp_double_vec& forceMask,
		const vcp_double_vec& epsRFInvrc3)
{
	const vcp_double_vec dx = _mm_sub_pd(r1_x, r2_x);
	const vcp_double_vec dy = _mm_sub_pd(r1_y, r2_y);
	const vcp_double_vec dz = _mm_sub_pd(r1_z, r2_z);

	const vcp_double_vec dx2 = _mm_mul_pd(dx, dx);
	const vcp_double_vec dy2 = _mm_mul_pd(dy, dy);
	const vcp_double_vec dz2 = _mm_mul_pd(dz, dz);

	const vcp_double_vec dr2 = _mm_add_pd(_mm_add_pd(dx2, dy2), dz2);

	const vcp_double_vec dr2_inv_unmasked = _mm_div_pd(one, dr2);
	const vcp_double_vec dr2_inv = vcp_simd_and(dr2_inv_unmasked, forceMask);
	const vcp_double_vec dr_inv = _mm_sqrt_pd(dr2_inv);
	const vcp_double_vec dr2three_inv = _mm_mul_pd(three, dr2_inv);

	const vcp_double_vec p1p2 = vcp_simd_and(forceMask, _mm_mul_pd(pii, pjj));
	const vcp_double_vec p1p2per4pie0 = p1p2;
	const vcp_double_vec rffac = _mm_mul_pd(p1p2, epsRFInvrc3);

	const vcp_double_vec p1p2per4pie0r3 = _mm_mul_pd(p1p2per4pie0, _mm_mul_pd(dr_inv, dr2_inv));
	const vcp_double_vec p1p2threeper4pie0r5 = _mm_mul_pd(p1p2per4pie0r3, dr2three_inv);

	const vcp_double_vec e1e2 = _mm_add_pd(_mm_mul_pd(eii_x, ejj_x), _mm_add_pd(_mm_mul_pd(eii_y, ejj_y), _mm_mul_pd(eii_z, ejj_z)));
	const vcp_double_vec re1 = _mm_add_pd(_mm_mul_pd(dx, eii_x), _mm_add_pd(_mm_mul_pd(dy, eii_y), _mm_mul_pd(dz, eii_z)));
	const vcp_double_vec re2 = _mm_add_pd(_mm_mul_pd(dx, ejj_x), _mm_add_pd(_mm_mul_pd(dy, ejj_y), _mm_mul_pd(dz, ejj_z)));

	const vcp_double_vec re1threeperr2 = _mm_mul_pd(re1, dr2three_inv);
	const vcp_double_vec re2threeperr2 = _mm_mul_pd(re2, dr2three_inv);
	const vcp_double_vec re1re2perr2 = _mm_mul_pd(dr2_inv, _mm_mul_pd(re1, re2));

	const vcp_double_vec e1e2minus5re1re2perr2 = _mm_sub_pd(e1e2, _mm_mul_pd(five, re1re2perr2));


	f_x = _mm_mul_pd(p1p2threeper4pie0r5, _mm_add_pd(_mm_mul_pd(dx, e1e2minus5re1re2perr2), _mm_add_pd(_mm_mul_pd(eii_x, re2), _mm_mul_pd(ejj_x, re1))));
	f_y = _mm_mul_pd(p1p2threeper4pie0r5, _mm_add_pd(_mm_mul_pd(dy, e1e2minus5re1re2perr2), _mm_add_pd(_mm_mul_pd(eii_y, re2), _mm_mul_pd(ejj_y, re1))));
	f_z = _mm_mul_pd(p1p2threeper4pie0r5, _mm_add_pd(_mm_mul_pd(dz, e1e2minus5re1re2perr2), _mm_add_pd(_mm_mul_pd(eii_z, re2), _mm_mul_pd(ejj_z, re1))));

	const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm_movemask_pd(macroMask) > 0) {
		// can we precompute some of this?
		const vcp_double_vec upot = _mm_mul_pd(p1p2per4pie0r3, _mm_sub_pd(e1e2, _mm_mul_pd(three, re1re2perr2)));
		const vcp_double_vec upot_masked = vcp_simd_and(upot, macroMask);
		sum_upotXpoles = _mm_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm_add_pd(_mm_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = vcp_simd_and(virial, macroMask);
		sum_virial = _mm_add_pd(sum_virial, virial_masked);

		const vcp_double_vec myRF_masked =  vcp_simd_and(macroMask, _mm_mul_pd(rffac, e1e2));
		sum_myRF = _mm_add_pd(sum_myRF, myRF_masked);
	}

	const vcp_double_vec e1_x_e2_y = _mm_mul_pd(eii_x, ejj_y);
	const vcp_double_vec e1_x_e2_z = _mm_mul_pd(eii_x, ejj_z);
	const vcp_double_vec e1_y_e2_x = _mm_mul_pd(eii_y, ejj_x);
	const vcp_double_vec e1_y_e2_z = _mm_mul_pd(eii_y, ejj_z);
	const vcp_double_vec e1_z_e2_x = _mm_mul_pd(eii_z, ejj_x);
	const vcp_double_vec e1_z_e2_y = _mm_mul_pd(eii_z, ejj_y);

	const vcp_double_vec e1_x_e2_y_minus_e1_y_e2_x = _mm_sub_pd(e1_x_e2_y, e1_y_e2_x);
	const vcp_double_vec e1_y_e2_z_minus_e1_z_e2_y = _mm_sub_pd(e1_y_e2_z, e1_z_e2_y);
	const vcp_double_vec e1_z_e2_x_minus_e1_x_e2_z = _mm_sub_pd(e1_z_e2_x, e1_x_e2_z);

	M1_x = _mm_add_pd(_mm_mul_pd(p1p2per4pie0r3, _mm_sub_pd(_mm_mul_pd(re2threeperr2, _mm_sub_pd(_mm_mul_pd(eii_y, dz), _mm_mul_pd(eii_z, dy))), e1_y_e2_z_minus_e1_z_e2_y)), _mm_mul_pd(rffac, e1_y_e2_z_minus_e1_z_e2_y));
	M1_y = _mm_add_pd(_mm_mul_pd(p1p2per4pie0r3, _mm_sub_pd(_mm_mul_pd(re2threeperr2, _mm_sub_pd(_mm_mul_pd(eii_z, dx), _mm_mul_pd(eii_x, dz))), e1_z_e2_x_minus_e1_x_e2_z)), _mm_mul_pd(rffac, e1_z_e2_x_minus_e1_x_e2_z));
	M1_z = _mm_add_pd(_mm_mul_pd(p1p2per4pie0r3, _mm_sub_pd(_mm_mul_pd(re2threeperr2, _mm_sub_pd(_mm_mul_pd(eii_x, dy), _mm_mul_pd(eii_y, dx))), e1_x_e2_y_minus_e1_y_e2_x)), _mm_mul_pd(rffac, e1_x_e2_y_minus_e1_y_e2_x));

	M2_x = _mm_sub_pd(_mm_mul_pd(p1p2per4pie0r3, _mm_add_pd(_mm_mul_pd(re1threeperr2, _mm_sub_pd(_mm_mul_pd(ejj_y, dz), _mm_mul_pd(ejj_z, dy))), e1_y_e2_z_minus_e1_z_e2_y)), _mm_mul_pd(rffac, e1_y_e2_z_minus_e1_z_e2_y));
	M2_y = _mm_sub_pd(_mm_mul_pd(p1p2per4pie0r3, _mm_add_pd(_mm_mul_pd(re1threeperr2, _mm_sub_pd(_mm_mul_pd(ejj_z, dx), _mm_mul_pd(ejj_x, dz))), e1_z_e2_x_minus_e1_x_e2_z)), _mm_mul_pd(rffac, e1_z_e2_x_minus_e1_x_e2_z));
	M2_z = _mm_sub_pd(_mm_mul_pd(p1p2per4pie0r3, _mm_add_pd(_mm_mul_pd(re1threeperr2, _mm_sub_pd(_mm_mul_pd(ejj_x, dy), _mm_mul_pd(ejj_y, dx))), e1_x_e2_y_minus_e1_y_e2_x)), _mm_mul_pd(rffac, e1_x_e2_y_minus_e1_y_e2_x));
}

template<class MacroPolicy>
inline void _loopBodyChargeQuadrupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& q,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& m,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched) {

	const vcp_double_vec c_dx = _mm_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm_add_pd(_mm_add_pd(c_dx2, c_dy2), c_dz2);

	const vcp_double_vec invdr2_unmasked = _mm_div_pd(one, c_dr2);
	const vcp_double_vec invdr2 = vcp_simd_and(invdr2_unmasked, forceMask);
	const vcp_double_vec invdr = _mm_sqrt_pd(invdr2);

	const vcp_double_vec qQ05per4pie0 = _mm_mul_pd(_05, _mm_mul_pd(q, m));

	const vcp_double_vec ejj_xXdx = _mm_mul_pd(ejj_x, c_dx);
	const vcp_double_vec ejj_yXdy = _mm_mul_pd(ejj_y, c_dy);
	const vcp_double_vec ejj_zXdz = _mm_mul_pd(ejj_z, c_dz);
	vcp_double_vec costj = _mm_add_pd(ejj_xXdx, _mm_add_pd(ejj_yXdy, ejj_zXdz));
	costj = _mm_mul_pd(costj, invdr);

	const vcp_double_vec qQinv4dr3 = _mm_mul_pd(qQ05per4pie0, _mm_mul_pd(invdr, invdr2));
	vcp_double_vec part1 = _mm_mul_pd(three, _mm_mul_pd(costj, costj));
	const vcp_double_vec upot = _mm_mul_pd(qQinv4dr3, _mm_sub_pd(part1, one));

	/**********
	 * Force
	 **********/
	const vcp_double_vec minus_partialRijInvdr = _mm_mul_pd(three, _mm_mul_pd(upot, invdr2));
	const vcp_double_vec partialTjInvdr = _mm_mul_pd(_mm_mul_pd(six, costj), _mm_mul_pd(qQinv4dr3, invdr));

	part1 = _mm_mul_pd(costj, _mm_mul_pd(partialTjInvdr, invdr));
	const vcp_double_vec fac = _mm_add_pd(part1, minus_partialRijInvdr);

	f_x = _mm_sub_pd(_mm_mul_pd(fac, c_dx), _mm_mul_pd(partialTjInvdr, ejj_x));
	f_y = _mm_sub_pd(_mm_mul_pd(fac, c_dy), _mm_mul_pd(partialTjInvdr, ejj_y));
	f_z = _mm_sub_pd(_mm_mul_pd(fac, c_dz), _mm_mul_pd(partialTjInvdr, ejj_z));

	const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);

	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm_movemask_pd(macroMask) > 0) {
		// do we have to mask "upot"? It should already have been masked by "qQinv4dr3" which
		// itself is masked by "invdr2".
		const vcp_double_vec upot_masked = vcp_simd_and(upot, macroMask);
		sum_upotXpoles = _mm_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm_add_pd(_mm_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = vcp_simd_and(virial, macroMask);
		sum_virial = _mm_add_pd(sum_virial, virial_masked);
	}

	/**********
	 * Torque
	 **********/
	const vcp_double_vec minuseXrij_x = _mm_sub_pd(_mm_mul_pd(ejj_z, c_dy), _mm_mul_pd(ejj_y, c_dz));
	const vcp_double_vec minuseXrij_y = _mm_sub_pd(_mm_mul_pd(ejj_x, c_dz), _mm_mul_pd(ejj_z, c_dx));
	const vcp_double_vec minuseXrij_z = _mm_sub_pd(_mm_mul_pd(ejj_y, c_dx), _mm_mul_pd(ejj_x, c_dy));

	M_x = _mm_mul_pd(partialTjInvdr, minuseXrij_x);
	M_y = _mm_mul_pd(partialTjInvdr, minuseXrij_y);
	M_z = _mm_mul_pd(partialTjInvdr, minuseXrij_z);
}

template<class MacroPolicy>
inline void _loopBodyDipoleQuadrupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& p,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& m,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched) {


	const vcp_double_vec c_dx = _mm_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm_add_pd(_mm_add_pd(c_dx2, c_dy2), c_dz2);

	const vcp_double_vec invdr2_unmasked = _mm_div_pd(one, c_dr2);
	const vcp_double_vec invdr2 = vcp_simd_and(invdr2_unmasked, forceMask);
	const vcp_double_vec invdr = _mm_sqrt_pd(invdr2);

	const vcp_double_vec myqfac = _mm_mul_pd(_mm_mul_pd(_1pt5, _mm_mul_pd(p, m)), _mm_mul_pd(invdr2, invdr2));

	const vcp_double_vec eii_xXdx = _mm_mul_pd(eii_x, c_dx);
	const vcp_double_vec eii_yXdy = _mm_mul_pd(eii_y, c_dy);
	const vcp_double_vec eii_zXdz = _mm_mul_pd(eii_z, c_dz);
	vcp_double_vec costi = _mm_add_pd(eii_xXdx, _mm_add_pd(eii_yXdy, eii_zXdz));
	costi = _mm_mul_pd(costi, invdr);

	const vcp_double_vec ejj_xXdx = _mm_mul_pd(ejj_x, c_dx);
	const vcp_double_vec ejj_yXdy = _mm_mul_pd(ejj_y, c_dy);
	const vcp_double_vec ejj_zXdz = _mm_mul_pd(ejj_z, c_dz);
	vcp_double_vec costj = _mm_add_pd(ejj_xXdx, _mm_add_pd(ejj_yXdy, ejj_zXdz));
	costj = _mm_mul_pd(costj, invdr);

	const vcp_double_vec cos2tj = _mm_mul_pd(costj, costj);

	const vcp_double_vec eiiXejj_x = _mm_mul_pd(eii_x, ejj_x);
	const vcp_double_vec eiiXejj_y = _mm_mul_pd(eii_y, ejj_y);
	const vcp_double_vec eiiXejj_z = _mm_mul_pd(eii_z, ejj_z);
	const vcp_double_vec cosgij = _mm_add_pd(eiiXejj_x, _mm_add_pd(eiiXejj_y, eiiXejj_z));

	/************
	 * Potential
	 ************/
	// TODO: Check if upot has to be multiplied by -1 according to DISS_STOLL S.178.
	// This affects also the implementation in potforce.h
	const vcp_double_vec _5cos2tjminus1 = _mm_sub_pd(_mm_mul_pd(five, cos2tj), one);
	const vcp_double_vec _2costj = _mm_mul_pd(two, costj);

	vcp_double_vec part1 = _mm_mul_pd(costi, _5cos2tjminus1);
	vcp_double_vec part2 = _mm_mul_pd(_2costj, cosgij);

	vcp_double_vec const upot = _mm_mul_pd(myqfac, _mm_sub_pd(part2, part1));

	const vcp_double_vec myqfacXinvdr = _mm_mul_pd(myqfac, invdr);
	const vcp_double_vec minus_partialRijInvdr = _mm_mul_pd(four, _mm_mul_pd(upot, invdr2));
	const vcp_double_vec minus_partialTiInvdr = _mm_mul_pd(myqfacXinvdr, _5cos2tjminus1);

	part1 = _mm_mul_pd(five, _mm_mul_pd(costi, costj));
	part1 = _mm_sub_pd(part1, cosgij); // *-1!

	const vcp_double_vec minus_partialTjInvdr = _mm_mul_pd(myqfacXinvdr, _mm_mul_pd(two, part1));
	const vcp_double_vec partialGij = _mm_mul_pd(myqfac, _2costj);

	part1 = _mm_mul_pd(costi, minus_partialTiInvdr);
	part2 = _mm_mul_pd(costj, minus_partialTjInvdr);
	vcp_double_vec part3 = _mm_add_pd(part1, part2);
	const vcp_double_vec fac = _mm_sub_pd(minus_partialRijInvdr, _mm_mul_pd(part3, invdr));

	// Force components
	part1 = _mm_mul_pd(fac, c_dx);
	part2 = _mm_mul_pd(minus_partialTiInvdr, eii_x);
	part3 = _mm_mul_pd(minus_partialTjInvdr, ejj_x);
	f_x = _mm_add_pd(part1, _mm_add_pd(part2, part3));

	part1 = _mm_mul_pd(fac, c_dy);
	part2 = _mm_mul_pd(minus_partialTiInvdr, eii_y);
	part3 = _mm_mul_pd(minus_partialTjInvdr, ejj_y);
	f_y = _mm_add_pd(part1, _mm_add_pd(part2, part3));

	part1 = _mm_mul_pd(fac, c_dz);
	part2 = _mm_mul_pd(minus_partialTiInvdr, eii_z);
	part3 = _mm_mul_pd(minus_partialTjInvdr, ejj_z);
	f_z = _mm_add_pd(part1, _mm_add_pd(part2, part3));

	const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);

	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm_movemask_pd(macroMask) > 0) {
		// do we have to mask "upot"? It should already have been masked by "myqfac" which
		// itself is masked by "invdr2".
		const vcp_double_vec upot_masked = vcp_simd_and(upot, macroMask);
		sum_upotXpoles = _mm_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm_add_pd(_mm_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = vcp_simd_and(virial, macroMask);
		sum_virial = _mm_add_pd(sum_virial, virial_masked);
	}

	/**********
	 * Torque
	 **********/
	const vcp_double_vec eii_x_ejj_y = _mm_mul_pd(eii_x, ejj_y);
	const vcp_double_vec eii_x_ejj_z = _mm_mul_pd(eii_x, ejj_z);
	const vcp_double_vec eii_y_ejj_x = _mm_mul_pd(eii_y, ejj_x);
	const vcp_double_vec eii_y_ejj_z = _mm_mul_pd(eii_y, ejj_z);
	const vcp_double_vec eii_z_ejj_x = _mm_mul_pd(eii_z, ejj_x);
	const vcp_double_vec eii_z_ejj_y = _mm_mul_pd(eii_z, ejj_y);

	const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = _mm_sub_pd(eii_x_ejj_y, eii_y_ejj_x);
	const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = _mm_sub_pd(eii_y_ejj_z, eii_z_ejj_y);
	const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = _mm_sub_pd(eii_z_ejj_x, eii_x_ejj_z);

	const vcp_double_vec partialGij_eiXej_x = _mm_mul_pd(partialGij, eii_y_ejj_z_minus_eii_z_ejj_y);
	const vcp_double_vec partialGij_eiXej_y = _mm_mul_pd(partialGij, eii_z_ejj_x_minus_eii_x_ejj_z);
	const vcp_double_vec partialGij_eiXej_z = _mm_mul_pd(partialGij, eii_x_ejj_y_minus_eii_y_ejj_x);

	vcp_double_vec eXrij_x = _mm_sub_pd(_mm_mul_pd(eii_y, c_dz), _mm_mul_pd(eii_z, c_dy));
	vcp_double_vec eXrij_y = _mm_sub_pd(_mm_mul_pd(eii_z, c_dx), _mm_mul_pd(eii_x, c_dz));
	vcp_double_vec eXrij_z = _mm_sub_pd(_mm_mul_pd(eii_x, c_dy), _mm_mul_pd(eii_y, c_dx));

	M1_x = _mm_sub_pd(_mm_mul_pd(minus_partialTiInvdr, eXrij_x), partialGij_eiXej_x);
	M1_y = _mm_sub_pd(_mm_mul_pd(minus_partialTiInvdr, eXrij_y), partialGij_eiXej_y);
	M1_z = _mm_sub_pd(_mm_mul_pd(minus_partialTiInvdr, eXrij_z), partialGij_eiXej_z);

	eXrij_x = _mm_sub_pd(_mm_mul_pd(ejj_y, c_dz), _mm_mul_pd(ejj_z, c_dy));
	eXrij_y = _mm_sub_pd(_mm_mul_pd(ejj_z, c_dx), _mm_mul_pd(ejj_x, c_dz));
	eXrij_z = _mm_sub_pd(_mm_mul_pd(ejj_x, c_dy), _mm_mul_pd(ejj_y, c_dx));

	M2_x = _mm_add_pd(_mm_mul_pd(minus_partialTjInvdr, eXrij_x), partialGij_eiXej_x);
	M2_y = _mm_add_pd(_mm_mul_pd(minus_partialTjInvdr, eXrij_y), partialGij_eiXej_y);
	M2_z = _mm_add_pd(_mm_mul_pd(minus_partialTjInvdr, eXrij_z), partialGij_eiXej_z);
}

template<class MacroPolicy>
inline void _loopBodyQudarupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& mii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& mjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& Mii_x, vcp_double_vec& Mii_y, vcp_double_vec& Mii_z,
		vcp_double_vec& Mjj_x, vcp_double_vec& Mjj_y, vcp_double_vec& Mjj_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask)
{
	const vcp_double_vec c_dx = _mm_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm_add_pd(_mm_add_pd(c_dx2, c_dy2), c_dz2);

	const vcp_double_vec invdr2_unmasked = _mm_div_pd(one, c_dr2);
	const vcp_double_vec invdr2 = vcp_simd_and(invdr2_unmasked, forceMask);
	const vcp_double_vec invdr = _mm_sqrt_pd(invdr2);

	vcp_double_vec qfac = _mm_mul_pd(_075, invdr);
	qfac = _mm_mul_pd(qfac, _mm_mul_pd(mii, mjj));
	qfac = _mm_mul_pd(qfac, _mm_mul_pd(invdr2, invdr2));

	const vcp_double_vec eii_xXdx = _mm_mul_pd(eii_x, c_dx);
	const vcp_double_vec eii_yXdy = _mm_mul_pd(eii_y, c_dy);
	const vcp_double_vec eii_zXdz = _mm_mul_pd(eii_z, c_dz);
	vcp_double_vec costi = _mm_add_pd(eii_xXdx, _mm_add_pd(eii_yXdy, eii_zXdz));
	costi = _mm_mul_pd(costi, invdr);

	const vcp_double_vec ejj_xXdx = _mm_mul_pd(ejj_x, c_dx);
	const vcp_double_vec ejj_yXdy = _mm_mul_pd(ejj_y, c_dy);
	const vcp_double_vec ejj_zXdz = _mm_mul_pd(ejj_z, c_dz);
	vcp_double_vec costj = _mm_add_pd(ejj_xXdx, _mm_add_pd(ejj_yXdy, ejj_zXdz));
	costj = _mm_mul_pd(costj, invdr);

	const vcp_double_vec cos2ti = _mm_mul_pd(costi, costi);
	const vcp_double_vec cos2tj = _mm_mul_pd(costj, costj);

	const vcp_double_vec eiiXejj_x = _mm_mul_pd(eii_x, ejj_x);
	const vcp_double_vec eiiXejj_y = _mm_mul_pd(eii_y, ejj_y);
	const vcp_double_vec eiiXejj_z = _mm_mul_pd(eii_z, ejj_z);
	const vcp_double_vec cosgij = _mm_add_pd(eiiXejj_x, _mm_add_pd(eiiXejj_y, eiiXejj_z));

	vcp_double_vec term = _mm_mul_pd(five, _mm_mul_pd(costi, costj));
	term = _mm_sub_pd(cosgij, term);

	/************
	 * Potential
	 ************/
	vcp_double_vec part1 = _mm_mul_pd(five, _mm_add_pd(cos2ti, cos2tj));
	vcp_double_vec part2 = _mm_mul_pd(_15, _mm_mul_pd(cos2ti, cos2tj));
	vcp_double_vec part3 = _mm_mul_pd(two, _mm_mul_pd(term, term));
	vcp_double_vec upot = _mm_add_pd(part1, part2);
	upot = _mm_sub_pd(_mm_add_pd(one, part3), upot);
	upot = _mm_mul_pd(qfac, upot);

	/**********
	 * Force
	 **********/
	const vcp_double_vec minus_partialRijInvdr = _mm_mul_pd(five, _mm_mul_pd(upot, invdr2));

	// partialTiInvdr & partialTjInvdr
	part1 = _mm_mul_pd(qfac, _mm_mul_pd(ten, invdr));
	part2 = _mm_mul_pd(two, term);

	// partialTiInvdr only
	part3 = _mm_mul_pd(three, _mm_mul_pd(costi, cos2tj));
	vcp_double_vec part4 = _mm_add_pd(costi, _mm_add_pd(part3, _mm_mul_pd(part2, costj)));
	const vcp_double_vec minus_partialTiInvdr = _mm_mul_pd(part1, part4);

	// partialTjInvdr only
	part3 = _mm_mul_pd(three, _mm_mul_pd(costj, cos2ti));
	part4 = _mm_add_pd(costj, _mm_add_pd(part3, _mm_mul_pd(part2, costi)));
	const vcp_double_vec minus_partialTjInvdr = _mm_mul_pd(part1, part4);

	const vcp_double_vec partialGij = _mm_mul_pd(qfac, _mm_mul_pd(four, term));

	// fac
	part1 = _mm_mul_pd(minus_partialTiInvdr, costi);
	part2 = _mm_mul_pd(minus_partialTjInvdr, costj);
	part3 = _mm_mul_pd(_mm_add_pd(part1, part2), invdr);
	const vcp_double_vec fac = _mm_sub_pd(minus_partialRijInvdr, part3);

	// Force components
	part1 = _mm_mul_pd(fac, c_dx);
	part2 = _mm_mul_pd(minus_partialTiInvdr, eii_x);
	part3 = _mm_mul_pd(minus_partialTjInvdr, ejj_x);
	f_x = _mm_add_pd(part1, _mm_add_pd(part2, part3));

	part1 = _mm_mul_pd(fac, c_dy);
	part2 = _mm_mul_pd(minus_partialTiInvdr, eii_y);
	part3 = _mm_mul_pd(minus_partialTjInvdr, ejj_y);
	f_y = _mm_add_pd(part1, _mm_add_pd(part2, part3));

	part1 = _mm_mul_pd(fac, c_dz);
	part2 = _mm_mul_pd(minus_partialTiInvdr, eii_z);
	part3 = _mm_mul_pd(minus_partialTjInvdr, ejj_z);
	f_z = _mm_add_pd(part1, _mm_add_pd(part2, part3));

	const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm_movemask_pd(macroMask) > 0) {
		// do we have to mask "upot"? It should already have been masked by "qfac" ...
		const vcp_double_vec upot_masked = vcp_simd_and(upot, macroMask);
		sum_upotXpoles = _mm_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm_add_pd(_mm_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = vcp_simd_and(virial, macroMask);
		sum_virial = _mm_add_pd(sum_virial, virial_masked);
	}

	const vcp_double_vec eii_x_ejj_y = _mm_mul_pd(eii_x, ejj_y);
	const vcp_double_vec eii_x_ejj_z = _mm_mul_pd(eii_x, ejj_z);
	const vcp_double_vec eii_y_ejj_x = _mm_mul_pd(eii_y, ejj_x);
	const vcp_double_vec eii_y_ejj_z = _mm_mul_pd(eii_y, ejj_z);
	const vcp_double_vec eii_z_ejj_x = _mm_mul_pd(eii_z, ejj_x);
	const vcp_double_vec eii_z_ejj_y = _mm_mul_pd(eii_z, ejj_y);

	const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = _mm_sub_pd(eii_x_ejj_y, eii_y_ejj_x);
	const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = _mm_sub_pd(eii_y_ejj_z, eii_z_ejj_y);
	const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = _mm_sub_pd(eii_z_ejj_x, eii_x_ejj_z);

	const vcp_double_vec partialGij_eiXej_x = _mm_mul_pd(partialGij, eii_y_ejj_z_minus_eii_z_ejj_y);
	const vcp_double_vec partialGij_eiXej_y = _mm_mul_pd(partialGij, eii_z_ejj_x_minus_eii_x_ejj_z);
	const vcp_double_vec partialGij_eiXej_z = _mm_mul_pd(partialGij, eii_x_ejj_y_minus_eii_y_ejj_x);

	vcp_double_vec eXrij_x = _mm_sub_pd(_mm_mul_pd(eii_y, c_dz), _mm_mul_pd(eii_z, c_dy));
	vcp_double_vec eXrij_y = _mm_sub_pd(_mm_mul_pd(eii_z, c_dx), _mm_mul_pd(eii_x, c_dz));
	vcp_double_vec eXrij_z = _mm_sub_pd(_mm_mul_pd(eii_x, c_dy), _mm_mul_pd(eii_y, c_dx));

	Mii_x = _mm_sub_pd(_mm_mul_pd(minus_partialTiInvdr, eXrij_x), partialGij_eiXej_x);
	Mii_y = _mm_sub_pd(_mm_mul_pd(minus_partialTiInvdr, eXrij_y), partialGij_eiXej_y);
	Mii_z = _mm_sub_pd(_mm_mul_pd(minus_partialTiInvdr, eXrij_z), partialGij_eiXej_z);

	eXrij_x = _mm_sub_pd(_mm_mul_pd(ejj_y, c_dz), _mm_mul_pd(ejj_z, c_dy));
	eXrij_y = _mm_sub_pd(_mm_mul_pd(ejj_z, c_dx), _mm_mul_pd(ejj_x, c_dz));
	eXrij_z = _mm_sub_pd(_mm_mul_pd(ejj_x, c_dy), _mm_mul_pd(ejj_y, c_dx));

	Mjj_x = _mm_add_pd(_mm_mul_pd(minus_partialTjInvdr, eXrij_x), partialGij_eiXej_x);
	Mjj_y = _mm_add_pd(_mm_mul_pd(minus_partialTjInvdr, eXrij_y), partialGij_eiXej_y);
	Mjj_z = _mm_add_pd(_mm_mul_pd(minus_partialTjInvdr, eXrij_z), partialGij_eiXej_z);
}

#elif VCP_VEC_TYPE==VCP_VEC_AVX

//const vcp_double_vec minus_one = _mm256_set1_pd(-1.0);
const vcp_double_vec zero = vcp_simd_zerov();
const vcp_double_vec one = _mm256_set1_pd(1.0);
const vcp_double_vec two = _mm256_set1_pd(2.0);
const vcp_double_vec three = _mm256_set1_pd(3.0);
const vcp_double_vec four = _mm256_set1_pd(4.0);
const vcp_double_vec five = _mm256_set1_pd(5.0);
const vcp_double_vec six = _mm256_set1_pd(6.0);
const vcp_double_vec ten = _mm256_set1_pd(10.0);
const vcp_double_vec _05 = _mm256_set1_pd(0.5);
const vcp_double_vec _075 = _mm256_set1_pd(0.75);
const vcp_double_vec _1pt5 = _mm256_set1_pd(1.5);
const vcp_double_vec _15 = _mm256_set1_pd(15.0);

static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);

template<class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyLJ(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& sum_upot6lj, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& e0s0,
		const vcp_double_vec& e1s1, const vcp_double_vec& e2s2, const vcp_double_vec& e3s3,
		const size_t& id_j0, const size_t& id_j1,  const size_t& id_j2,  const size_t& id_j3,
		const size_t& id_i)
{
	const vcp_double_vec c_dx = _mm256_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm256_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm256_sub_pd(r1_z, r2_z);
	const vcp_double_vec c_dxdx = _mm256_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dydy = _mm256_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dzdz = _mm256_mul_pd(c_dz, c_dz);
	const vcp_double_vec c_dxdx_dydy = _mm256_add_pd(c_dxdx, c_dydy);
	const vcp_double_vec c_r2 = _mm256_add_pd(c_dxdx_dydy, c_dzdz);
	const vcp_double_vec r2_inv_unmasked = _mm256_div_pd(one, c_r2);
	const vcp_double_vec r2_inv = _mm256_and_pd(r2_inv_unmasked, forceMask);


	const vcp_double_vec e0e1 = _mm256_unpacklo_pd(e0s0, e1s1);
	const vcp_double_vec s0s1 = _mm256_unpackhi_pd(e0s0, e1s1);
	const vcp_double_vec e2e3 = _mm256_unpacklo_pd(e2s2, e3s3);
	const vcp_double_vec s2s3 = _mm256_unpackhi_pd(e2s2, e3s3);

	const vcp_double_vec eps_24 = _mm256_permute2f128_pd(e0e1, e2e3, 1<<5);
	const vcp_double_vec sig2 = _mm256_permute2f128_pd(s0s1, s2s3, 1<<5);

	const vcp_double_vec lj2 = _mm256_mul_pd(sig2, r2_inv);
	const vcp_double_vec lj4 = _mm256_mul_pd(lj2, lj2);
	const vcp_double_vec lj6 = _mm256_mul_pd(lj4, lj2);
	const vcp_double_vec lj12 = _mm256_mul_pd(lj6, lj6);
	const vcp_double_vec lj12m6 = _mm256_sub_pd(lj12, lj6);

	const vcp_double_vec eps24r2inv = _mm256_mul_pd(eps_24, r2_inv);
	const vcp_double_vec lj12lj12m6 = _mm256_add_pd(lj12, lj12m6);
	const vcp_double_vec scale = _mm256_mul_pd(eps24r2inv, lj12lj12m6);

	const vcp_double_vec fx = _mm256_mul_pd(c_dx, scale);
	const vcp_double_vec fy = _mm256_mul_pd(c_dy, scale);
	const vcp_double_vec fz = _mm256_mul_pd(c_dz, scale);


	const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);

	const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);

	const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);

	// Only go on if at least 1 macroscopic value has to be calculated.
	if (_mm256_movemask_pd(macroMask) > 0) {
		const vcp_double_vec sh0 = _mm256_maskload_pd(_shift6[id_i] + id_j0, memoryMask_first);
		const vcp_double_vec sh1 = _mm256_maskload_pd(_shift6[id_i] + id_j1, memoryMask_first);
		const vcp_double_vec sh2 = _mm256_maskload_pd(_shift6[id_i] + id_j2, memoryMask_first);
		const vcp_double_vec sh3 = _mm256_maskload_pd(_shift6[id_i] + id_j3, memoryMask_first);

		const vcp_double_vec sh0sh1 = _mm256_unpacklo_pd(sh0, sh1);
		const vcp_double_vec sh2sh3 = _mm256_unpacklo_pd(sh2, sh3);

		const vcp_double_vec shift6 = _mm256_permute2f128_pd(sh0sh1, sh2sh3, 1<<5);

		const vcp_double_vec upot = _mm256_mul_pd(eps_24, lj12m6);
		const vcp_double_vec upot_sh = _mm256_add_pd(shift6, upot);
		const vcp_double_vec upot_masked = _mm256_and_pd(upot_sh, macroMask);

		sum_upot6lj = _mm256_add_pd(sum_upot6lj, upot_masked);

		const vcp_double_vec vir_x = _mm256_mul_pd(m_dx, fx);
		const vcp_double_vec vir_y = _mm256_mul_pd(m_dy, fy);
		const vcp_double_vec vir_z = _mm256_mul_pd(m_dz, fz);

		const vcp_double_vec vir_xy = _mm256_add_pd(vir_x, vir_y);
		const vcp_double_vec virial = _mm256_add_pd(vir_xy, vir_z);
		const vcp_double_vec vir_masked = _mm256_and_pd(virial, macroMask);

		sum_virial = _mm256_add_pd(sum_virial, vir_masked);
	}

}



template<class MacroPolicy>
inline void _loopBodyCharge(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& qii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& qjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask)
{
	const vcp_double_vec c_dx = _mm256_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm256_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm256_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm256_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm256_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm256_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm256_add_pd(_mm256_add_pd(c_dx2, c_dy2), c_dz2);
	const vcp_double_vec c_dr2_inv_unmasked = _mm256_div_pd(one, c_dr2);
	const vcp_double_vec c_dr2_inv = _mm256_and_pd(c_dr2_inv_unmasked, forceMask);
	const vcp_double_vec c_dr_inv = _mm256_sqrt_pd(c_dr2_inv);

	const vcp_double_vec q1q2per4pie0 = _mm256_mul_pd(qii, qjj);
	const vcp_double_vec upot = _mm256_mul_pd(q1q2per4pie0, c_dr_inv);
	const vcp_double_vec fac = _mm256_mul_pd(upot, c_dr2_inv);

	f_x = _mm256_mul_pd(c_dx, fac);
	f_y = _mm256_mul_pd(c_dy, fac);
	f_z = _mm256_mul_pd(c_dz, fac);

	const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm256_movemask_pd(macroMask) > 0) {
		const vcp_double_vec upot_masked = _mm256_and_pd(upot, macroMask);
		sum_upotXpoles = _mm256_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm256_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm256_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm256_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm256_add_pd(_mm256_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = _mm256_and_pd(virial, macroMask);
		sum_virial = _mm256_add_pd(sum_virial, virial_masked);

	}
}

template<class MacroPolicy>
inline void _loopBodyChargeDipole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& q,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& e_x, const vcp_double_vec& e_y, const vcp_double_vec& e_z,
		const vcp_double_vec& p,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched)
{
	const vcp_double_vec dx = _mm256_sub_pd(r1_x, r2_x);
	const vcp_double_vec dy = _mm256_sub_pd(r1_y, r2_y);
	const vcp_double_vec dz = _mm256_sub_pd(r1_z, r2_z);

	const vcp_double_vec dx2 = _mm256_mul_pd(dx, dx);
	const vcp_double_vec dy2 = _mm256_mul_pd(dy, dy);
	const vcp_double_vec dz2 = _mm256_mul_pd(dz, dz);

	const vcp_double_vec dr2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);

	const vcp_double_vec dr2_inv_unmasked = _mm256_div_pd(one, dr2);
	const vcp_double_vec dr2_inv = _mm256_and_pd(dr2_inv_unmasked, forceMask);
	const vcp_double_vec dr_inv = _mm256_sqrt_pd(dr2_inv);
	const vcp_double_vec dr3_inv = _mm256_mul_pd(dr2_inv, dr_inv);

	const vcp_double_vec re = _mm256_add_pd(_mm256_mul_pd(dx, e_x), _mm256_add_pd(_mm256_mul_pd(dy, e_y), _mm256_mul_pd(dz, e_z)));

	const vcp_double_vec qpper4pie0 = _mm256_mul_pd(q, p);
	const vcp_double_vec qpper4pie0dr3 = _mm256_mul_pd(qpper4pie0, dr3_inv);

	// rename ?
	const vcp_double_vec fac = _mm256_mul_pd(dr2_inv, _mm256_mul_pd(three, re));

	f_x = _mm256_mul_pd(qpper4pie0dr3, _mm256_sub_pd(e_x, _mm256_mul_pd(dx, fac)));
	f_y = _mm256_mul_pd(qpper4pie0dr3, _mm256_sub_pd(e_y, _mm256_mul_pd(dy, fac)));
	f_z = _mm256_mul_pd(qpper4pie0dr3, _mm256_sub_pd(e_z, _mm256_mul_pd(dz, fac)));

	const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm256_movemask_pd(macroMask) > 0)
	{
		const vcp_double_vec minusUpot_unmasked =  _mm256_mul_pd(qpper4pie0dr3, re);
		const vcp_double_vec minusUpot = _mm256_and_pd(minusUpot_unmasked, macroMask);
		sum_upotXpoles = _mm256_sub_pd(sum_upotXpoles, minusUpot);

		const vcp_double_vec virial_x = _mm256_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm256_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm256_mul_pd(m_dz, f_z);

		const vcp_double_vec virial_unmasked = _mm256_add_pd(_mm256_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial = _mm256_and_pd(virial_unmasked, macroMask);
		sum_virial = _mm256_add_pd(sum_virial, virial);
	}

	const vcp_double_vec e_x_dy = _mm256_mul_pd(e_x, dy);
	const vcp_double_vec e_x_dz = _mm256_mul_pd(e_x, dz);
	const vcp_double_vec e_y_dx = _mm256_mul_pd(e_y, dx);
	const vcp_double_vec e_y_dz = _mm256_mul_pd(e_y, dz);
	const vcp_double_vec e_z_dx = _mm256_mul_pd(e_z, dx);
	const vcp_double_vec e_z_dy = _mm256_mul_pd(e_z, dy);

	const vcp_double_vec e_x_dy_minus_e_y_dx = _mm256_sub_pd(e_x_dy, e_y_dx);
	const vcp_double_vec e_y_dz_minus_e_z_dy = _mm256_sub_pd(e_y_dz, e_z_dy);
	const vcp_double_vec e_z_dx_minus_e_x_dz = _mm256_sub_pd(e_z_dx, e_x_dz);

	M_x = _mm256_mul_pd(qpper4pie0dr3, e_y_dz_minus_e_z_dy);
	M_y = _mm256_mul_pd(qpper4pie0dr3, e_z_dx_minus_e_x_dz);
	M_z = _mm256_mul_pd(qpper4pie0dr3, e_x_dy_minus_e_y_dx);
}

template<class MacroPolicy>
inline void _loopBodyDipole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& pii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& pjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial, vcp_double_vec& sum_myRF,
		const vcp_double_vec& forceMask,
		const vcp_double_vec& epsRFInvrc3)
{
	const vcp_double_vec dx = _mm256_sub_pd(r1_x, r2_x);
	const vcp_double_vec dy = _mm256_sub_pd(r1_y, r2_y);
	const vcp_double_vec dz = _mm256_sub_pd(r1_z, r2_z);

	const vcp_double_vec dx2 = _mm256_mul_pd(dx, dx);
	const vcp_double_vec dy2 = _mm256_mul_pd(dy, dy);
	const vcp_double_vec dz2 = _mm256_mul_pd(dz, dz);

	const vcp_double_vec dr2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);

	const vcp_double_vec dr2_inv_unmasked = _mm256_div_pd(one, dr2);
	const vcp_double_vec dr2_inv = _mm256_and_pd(dr2_inv_unmasked, forceMask);
	const vcp_double_vec dr_inv = _mm256_sqrt_pd(dr2_inv);
	const vcp_double_vec dr2three_inv = _mm256_mul_pd(three, dr2_inv);

	const vcp_double_vec p1p2 = _mm256_and_pd(forceMask, _mm256_mul_pd(pii, pjj));
	const vcp_double_vec p1p2per4pie0 = p1p2;
	const vcp_double_vec rffac = _mm256_mul_pd(p1p2, epsRFInvrc3);

	const vcp_double_vec p1p2per4pie0r3 = _mm256_mul_pd(p1p2per4pie0, _mm256_mul_pd(dr_inv, dr2_inv));
	const vcp_double_vec p1p2threeper4pie0r5 = _mm256_mul_pd(p1p2per4pie0r3, dr2three_inv);

	const vcp_double_vec e1e2 = _mm256_add_pd(_mm256_mul_pd(eii_x, ejj_x), _mm256_add_pd(_mm256_mul_pd(eii_y, ejj_y), _mm256_mul_pd(eii_z, ejj_z)));
	const vcp_double_vec re1 = _mm256_add_pd(_mm256_mul_pd(dx, eii_x), _mm256_add_pd(_mm256_mul_pd(dy, eii_y), _mm256_mul_pd(dz, eii_z)));
	const vcp_double_vec re2 = _mm256_add_pd(_mm256_mul_pd(dx, ejj_x), _mm256_add_pd(_mm256_mul_pd(dy, ejj_y), _mm256_mul_pd(dz, ejj_z)));

	const vcp_double_vec re1threeperr2 = _mm256_mul_pd(re1, dr2three_inv);
	const vcp_double_vec re2threeperr2 = _mm256_mul_pd(re2, dr2three_inv);
	const vcp_double_vec re1re2perr2 = _mm256_mul_pd(dr2_inv, _mm256_mul_pd(re1, re2));

	const vcp_double_vec e1e2minus5re1re2perr2 = _mm256_sub_pd(e1e2, _mm256_mul_pd(five, re1re2perr2));


	f_x = _mm256_mul_pd(p1p2threeper4pie0r5, _mm256_add_pd(_mm256_mul_pd(dx, e1e2minus5re1re2perr2), _mm256_add_pd(_mm256_mul_pd(eii_x, re2), _mm256_mul_pd(ejj_x, re1))));
	f_y = _mm256_mul_pd(p1p2threeper4pie0r5, _mm256_add_pd(_mm256_mul_pd(dy, e1e2minus5re1re2perr2), _mm256_add_pd(_mm256_mul_pd(eii_y, re2), _mm256_mul_pd(ejj_y, re1))));
	f_z = _mm256_mul_pd(p1p2threeper4pie0r5, _mm256_add_pd(_mm256_mul_pd(dz, e1e2minus5re1re2perr2), _mm256_add_pd(_mm256_mul_pd(eii_z, re2), _mm256_mul_pd(ejj_z, re1))));

	const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm256_movemask_pd(macroMask) > 0) {
		// can we precompute some of this?
		const vcp_double_vec upot = _mm256_mul_pd(p1p2per4pie0r3, _mm256_sub_pd(e1e2, _mm256_mul_pd(three, re1re2perr2)));
		const vcp_double_vec upot_masked = _mm256_and_pd(upot, macroMask);
		sum_upotXpoles = _mm256_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm256_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm256_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm256_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm256_add_pd(_mm256_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = _mm256_and_pd(virial, macroMask);
		sum_virial = _mm256_add_pd(sum_virial, virial_masked);

		const vcp_double_vec myRF_masked =  _mm256_and_pd(macroMask, _mm256_mul_pd(rffac, e1e2));
		sum_myRF = _mm256_add_pd(sum_myRF, myRF_masked);
	}

	const vcp_double_vec e1_x_e2_y = _mm256_mul_pd(eii_x, ejj_y);
	const vcp_double_vec e1_x_e2_z = _mm256_mul_pd(eii_x, ejj_z);
	const vcp_double_vec e1_y_e2_x = _mm256_mul_pd(eii_y, ejj_x);
	const vcp_double_vec e1_y_e2_z = _mm256_mul_pd(eii_y, ejj_z);
	const vcp_double_vec e1_z_e2_x = _mm256_mul_pd(eii_z, ejj_x);
	const vcp_double_vec e1_z_e2_y = _mm256_mul_pd(eii_z, ejj_y);

	const vcp_double_vec e1_x_e2_y_minus_e1_y_e2_x = _mm256_sub_pd(e1_x_e2_y, e1_y_e2_x);
	const vcp_double_vec e1_y_e2_z_minus_e1_z_e2_y = _mm256_sub_pd(e1_y_e2_z, e1_z_e2_y);
	const vcp_double_vec e1_z_e2_x_minus_e1_x_e2_z = _mm256_sub_pd(e1_z_e2_x, e1_x_e2_z);

	M1_x = _mm256_add_pd(_mm256_mul_pd(p1p2per4pie0r3, _mm256_sub_pd(_mm256_mul_pd(re2threeperr2, _mm256_sub_pd(_mm256_mul_pd(eii_y, dz), _mm256_mul_pd(eii_z, dy))), e1_y_e2_z_minus_e1_z_e2_y)), _mm256_mul_pd(rffac, e1_y_e2_z_minus_e1_z_e2_y));
	M1_y = _mm256_add_pd(_mm256_mul_pd(p1p2per4pie0r3, _mm256_sub_pd(_mm256_mul_pd(re2threeperr2, _mm256_sub_pd(_mm256_mul_pd(eii_z, dx), _mm256_mul_pd(eii_x, dz))), e1_z_e2_x_minus_e1_x_e2_z)), _mm256_mul_pd(rffac, e1_z_e2_x_minus_e1_x_e2_z));
	M1_z = _mm256_add_pd(_mm256_mul_pd(p1p2per4pie0r3, _mm256_sub_pd(_mm256_mul_pd(re2threeperr2, _mm256_sub_pd(_mm256_mul_pd(eii_x, dy), _mm256_mul_pd(eii_y, dx))), e1_x_e2_y_minus_e1_y_e2_x)), _mm256_mul_pd(rffac, e1_x_e2_y_minus_e1_y_e2_x));

	M2_x = _mm256_sub_pd(_mm256_mul_pd(p1p2per4pie0r3, _mm256_add_pd(_mm256_mul_pd(re1threeperr2, _mm256_sub_pd(_mm256_mul_pd(ejj_y, dz), _mm256_mul_pd(ejj_z, dy))), e1_y_e2_z_minus_e1_z_e2_y)), _mm256_mul_pd(rffac, e1_y_e2_z_minus_e1_z_e2_y));
	M2_y = _mm256_sub_pd(_mm256_mul_pd(p1p2per4pie0r3, _mm256_add_pd(_mm256_mul_pd(re1threeperr2, _mm256_sub_pd(_mm256_mul_pd(ejj_z, dx), _mm256_mul_pd(ejj_x, dz))), e1_z_e2_x_minus_e1_x_e2_z)), _mm256_mul_pd(rffac, e1_z_e2_x_minus_e1_x_e2_z));
	M2_z = _mm256_sub_pd(_mm256_mul_pd(p1p2per4pie0r3, _mm256_add_pd(_mm256_mul_pd(re1threeperr2, _mm256_sub_pd(_mm256_mul_pd(ejj_x, dy), _mm256_mul_pd(ejj_y, dx))), e1_x_e2_y_minus_e1_y_e2_x)), _mm256_mul_pd(rffac, e1_x_e2_y_minus_e1_y_e2_x));
}

template<class MacroPolicy>
inline void _loopBodyChargeQuadrupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& q,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& m,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched) {


	const vcp_double_vec c_dx = _mm256_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm256_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm256_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm256_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm256_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm256_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm256_add_pd(_mm256_add_pd(c_dx2, c_dy2), c_dz2);

	const vcp_double_vec invdr2_unmasked = _mm256_div_pd(one, c_dr2);
	const vcp_double_vec invdr2 = _mm256_and_pd(invdr2_unmasked, forceMask);
	const vcp_double_vec invdr = _mm256_sqrt_pd(invdr2);

	const vcp_double_vec qQ05per4pie0 = _mm256_mul_pd(_05, _mm256_mul_pd(q, m));

	const vcp_double_vec ejj_xXdx = _mm256_mul_pd(ejj_x, c_dx);
	const vcp_double_vec ejj_yXdy = _mm256_mul_pd(ejj_y, c_dy);
	const vcp_double_vec ejj_zXdz = _mm256_mul_pd(ejj_z, c_dz);
	vcp_double_vec costj = _mm256_add_pd(ejj_xXdx, _mm256_add_pd(ejj_yXdy, ejj_zXdz));
	costj = _mm256_mul_pd(costj, invdr);

	const vcp_double_vec qQinv4dr3 = _mm256_mul_pd(qQ05per4pie0, _mm256_mul_pd(invdr, invdr2));
	vcp_double_vec part1 = _mm256_mul_pd(three, _mm256_mul_pd(costj, costj));
	const vcp_double_vec upot = _mm256_mul_pd(qQinv4dr3, _mm256_sub_pd(part1, one));

	/**********
	 * Force
	 **********/
	const vcp_double_vec minus_partialRijInvdr = _mm256_mul_pd(three, _mm256_mul_pd(upot, invdr2));
	const vcp_double_vec partialTjInvdr = _mm256_mul_pd(_mm256_mul_pd(six, costj), _mm256_mul_pd(qQinv4dr3, invdr));

	part1 = _mm256_mul_pd(costj, _mm256_mul_pd(partialTjInvdr, invdr));
	const vcp_double_vec fac = _mm256_add_pd(part1, minus_partialRijInvdr);

	f_x = _mm256_sub_pd(_mm256_mul_pd(fac, c_dx), _mm256_mul_pd(partialTjInvdr, ejj_x));
	f_y = _mm256_sub_pd(_mm256_mul_pd(fac, c_dy), _mm256_mul_pd(partialTjInvdr, ejj_y));
	f_z = _mm256_sub_pd(_mm256_mul_pd(fac, c_dz), _mm256_mul_pd(partialTjInvdr, ejj_z));

	const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);

	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm256_movemask_pd(macroMask) > 0) {
		// do we have to mask "upot"? ...
		const vcp_double_vec upot_masked = _mm256_and_pd(upot, macroMask);
		sum_upotXpoles = _mm256_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm256_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm256_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm256_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm256_add_pd(_mm256_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = _mm256_and_pd(virial, macroMask);
		sum_virial = _mm256_add_pd(sum_virial, virial_masked);
	}

	/**********
	 * Torque
	 **********/
	const vcp_double_vec minuseXrij_x = _mm256_sub_pd(_mm256_mul_pd(ejj_z, c_dy), _mm256_mul_pd(ejj_y, c_dz));
	const vcp_double_vec minuseXrij_y = _mm256_sub_pd(_mm256_mul_pd(ejj_x, c_dz), _mm256_mul_pd(ejj_z, c_dx));
	const vcp_double_vec minuseXrij_z = _mm256_sub_pd(_mm256_mul_pd(ejj_y, c_dx), _mm256_mul_pd(ejj_x, c_dy));

	M_x = _mm256_mul_pd(partialTjInvdr, minuseXrij_x);
	M_y = _mm256_mul_pd(partialTjInvdr, minuseXrij_y);
	M_z = _mm256_mul_pd(partialTjInvdr, minuseXrij_z);
}

template<class MacroPolicy>
inline void _loopBodyDipoleQuadrupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& p,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& m,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched) {


	const vcp_double_vec c_dx = _mm256_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm256_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm256_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm256_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm256_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm256_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm256_add_pd(_mm256_add_pd(c_dx2, c_dy2), c_dz2);

	const vcp_double_vec invdr2_unmasked = _mm256_div_pd(one, c_dr2);
	const vcp_double_vec invdr2 = _mm256_and_pd(invdr2_unmasked, forceMask);
	const vcp_double_vec invdr = _mm256_sqrt_pd(invdr2);

	const vcp_double_vec myqfac = _mm256_mul_pd(_mm256_mul_pd(_1pt5, _mm256_mul_pd(p, m)), _mm256_mul_pd(invdr2, invdr2));

	const vcp_double_vec eii_xXdx = _mm256_mul_pd(eii_x, c_dx);
	const vcp_double_vec eii_yXdy = _mm256_mul_pd(eii_y, c_dy);
	const vcp_double_vec eii_zXdz = _mm256_mul_pd(eii_z, c_dz);
	vcp_double_vec costi = _mm256_add_pd(eii_xXdx, _mm256_add_pd(eii_yXdy, eii_zXdz));
	costi = _mm256_mul_pd(costi, invdr);

	const vcp_double_vec ejj_xXdx = _mm256_mul_pd(ejj_x, c_dx);
	const vcp_double_vec ejj_yXdy = _mm256_mul_pd(ejj_y, c_dy);
	const vcp_double_vec ejj_zXdz = _mm256_mul_pd(ejj_z, c_dz);
	vcp_double_vec costj = _mm256_add_pd(ejj_xXdx, _mm256_add_pd(ejj_yXdy, ejj_zXdz));
	costj = _mm256_mul_pd(costj, invdr);

	const vcp_double_vec cos2tj = _mm256_mul_pd(costj, costj);

	const vcp_double_vec eiiXejj_x = _mm256_mul_pd(eii_x, ejj_x);
	const vcp_double_vec eiiXejj_y = _mm256_mul_pd(eii_y, ejj_y);
	const vcp_double_vec eiiXejj_z = _mm256_mul_pd(eii_z, ejj_z);
	const vcp_double_vec cosgij = _mm256_add_pd(eiiXejj_x, _mm256_add_pd(eiiXejj_y, eiiXejj_z));

	/************
	 * Potential
	 ************/
	// TODO: Check if upot has to be multiplied by -1 according to DISS_STOLL S.178.
	// This affects also the implementation in potforce.h
	const vcp_double_vec _5cos2tjminus1 = _mm256_sub_pd(_mm256_mul_pd(five, cos2tj), one);
	const vcp_double_vec _2costj = _mm256_mul_pd(two, costj);

	vcp_double_vec part1 = _mm256_mul_pd(costi, _5cos2tjminus1);
	vcp_double_vec part2 = _mm256_mul_pd(_2costj, cosgij);

	vcp_double_vec const upot = _mm256_mul_pd(myqfac, _mm256_sub_pd(part2, part1));

	const vcp_double_vec myqfacXinvdr = _mm256_mul_pd(myqfac, invdr);
	const vcp_double_vec minus_partialRijInvdr = _mm256_mul_pd(four, _mm256_mul_pd(upot, invdr2));
	const vcp_double_vec minus_partialTiInvdr = _mm256_mul_pd(myqfacXinvdr, _5cos2tjminus1);

	part1 = _mm256_mul_pd(five, _mm256_mul_pd(costi, costj));
	part1 = _mm256_sub_pd(part1, cosgij); // *-1!

	const vcp_double_vec minus_partialTjInvdr = _mm256_mul_pd(myqfacXinvdr, _mm256_mul_pd(two, part1));
	const vcp_double_vec partialGij = _mm256_mul_pd(myqfac, _2costj);

	part1 = _mm256_mul_pd(costi, minus_partialTiInvdr);
	part2 = _mm256_mul_pd(costj, minus_partialTjInvdr);
	vcp_double_vec part3 = _mm256_add_pd(part1, part2);
	const vcp_double_vec fac = _mm256_sub_pd(minus_partialRijInvdr, _mm256_mul_pd(part3, invdr));

	// Force components
	part1 = _mm256_mul_pd(fac, c_dx);
	part2 = _mm256_mul_pd(minus_partialTiInvdr, eii_x);
	part3 = _mm256_mul_pd(minus_partialTjInvdr, ejj_x);
	f_x = _mm256_add_pd(part1, _mm256_add_pd(part2, part3));

	part1 = _mm256_mul_pd(fac, c_dy);
	part2 = _mm256_mul_pd(minus_partialTiInvdr, eii_y);
	part3 = _mm256_mul_pd(minus_partialTjInvdr, ejj_y);
	f_y = _mm256_add_pd(part1, _mm256_add_pd(part2, part3));

	part1 = _mm256_mul_pd(fac, c_dz);
	part2 = _mm256_mul_pd(minus_partialTiInvdr, eii_z);
	part3 = _mm256_mul_pd(minus_partialTjInvdr, ejj_z);
	f_z = _mm256_add_pd(part1, _mm256_add_pd(part2, part3));

	const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);

	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm256_movemask_pd(macroMask) > 0) {
		// do we have to mask "upot"? ...
		const vcp_double_vec upot_masked = _mm256_and_pd(upot, macroMask);
		sum_upotXpoles = _mm256_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm256_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm256_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm256_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm256_add_pd(_mm256_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = _mm256_and_pd(virial, macroMask);
		sum_virial = _mm256_add_pd(sum_virial, virial_masked);
	}

	/**********
	 * Torque
	 **********/
	const vcp_double_vec eii_x_ejj_y = _mm256_mul_pd(eii_x, ejj_y);
	const vcp_double_vec eii_x_ejj_z = _mm256_mul_pd(eii_x, ejj_z);
	const vcp_double_vec eii_y_ejj_x = _mm256_mul_pd(eii_y, ejj_x);
	const vcp_double_vec eii_y_ejj_z = _mm256_mul_pd(eii_y, ejj_z);
	const vcp_double_vec eii_z_ejj_x = _mm256_mul_pd(eii_z, ejj_x);
	const vcp_double_vec eii_z_ejj_y = _mm256_mul_pd(eii_z, ejj_y);

	const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = _mm256_sub_pd(eii_x_ejj_y, eii_y_ejj_x);
	const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = _mm256_sub_pd(eii_y_ejj_z, eii_z_ejj_y);
	const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = _mm256_sub_pd(eii_z_ejj_x, eii_x_ejj_z);

	const vcp_double_vec partialGij_eiXej_x = _mm256_mul_pd(partialGij, eii_y_ejj_z_minus_eii_z_ejj_y);
	const vcp_double_vec partialGij_eiXej_y = _mm256_mul_pd(partialGij, eii_z_ejj_x_minus_eii_x_ejj_z);
	const vcp_double_vec partialGij_eiXej_z = _mm256_mul_pd(partialGij, eii_x_ejj_y_minus_eii_y_ejj_x);

	vcp_double_vec eXrij_x = _mm256_sub_pd(_mm256_mul_pd(eii_y, c_dz), _mm256_mul_pd(eii_z, c_dy));
	vcp_double_vec eXrij_y = _mm256_sub_pd(_mm256_mul_pd(eii_z, c_dx), _mm256_mul_pd(eii_x, c_dz));
	vcp_double_vec eXrij_z = _mm256_sub_pd(_mm256_mul_pd(eii_x, c_dy), _mm256_mul_pd(eii_y, c_dx));

	M1_x = _mm256_sub_pd(_mm256_mul_pd(minus_partialTiInvdr, eXrij_x), partialGij_eiXej_x);
	M1_y = _mm256_sub_pd(_mm256_mul_pd(minus_partialTiInvdr, eXrij_y), partialGij_eiXej_y);
	M1_z = _mm256_sub_pd(_mm256_mul_pd(minus_partialTiInvdr, eXrij_z), partialGij_eiXej_z);

	eXrij_x = _mm256_sub_pd(_mm256_mul_pd(ejj_y, c_dz), _mm256_mul_pd(ejj_z, c_dy));
	eXrij_y = _mm256_sub_pd(_mm256_mul_pd(ejj_z, c_dx), _mm256_mul_pd(ejj_x, c_dz));
	eXrij_z = _mm256_sub_pd(_mm256_mul_pd(ejj_x, c_dy), _mm256_mul_pd(ejj_y, c_dx));

	M2_x = _mm256_add_pd(_mm256_mul_pd(minus_partialTjInvdr, eXrij_x), partialGij_eiXej_x);
	M2_y = _mm256_add_pd(_mm256_mul_pd(minus_partialTjInvdr, eXrij_y), partialGij_eiXej_y);
	M2_z = _mm256_add_pd(_mm256_mul_pd(minus_partialTjInvdr, eXrij_z), partialGij_eiXej_z);
}



template<class MacroPolicy>
inline void _loopBodyQudarupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& mii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& mjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& Mii_x, vcp_double_vec& Mii_y, vcp_double_vec& Mii_z,
		vcp_double_vec& Mjj_x, vcp_double_vec& Mjj_y, vcp_double_vec& Mjj_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask)
{
	const vcp_double_vec c_dx = _mm256_sub_pd(r1_x, r2_x);
	const vcp_double_vec c_dy = _mm256_sub_pd(r1_y, r2_y);
	const vcp_double_vec c_dz = _mm256_sub_pd(r1_z, r2_z);

	const vcp_double_vec c_dx2 = _mm256_mul_pd(c_dx, c_dx);
	const vcp_double_vec c_dy2 = _mm256_mul_pd(c_dy, c_dy);
	const vcp_double_vec c_dz2 = _mm256_mul_pd(c_dz, c_dz);

	const vcp_double_vec c_dr2 = _mm256_add_pd(_mm256_add_pd(c_dx2, c_dy2), c_dz2);

	const vcp_double_vec invdr2_unmasked = _mm256_div_pd(one, c_dr2);
	const vcp_double_vec invdr2 = _mm256_and_pd(invdr2_unmasked, forceMask);
	const vcp_double_vec invdr = _mm256_sqrt_pd(invdr2);

	vcp_double_vec qfac = _mm256_mul_pd(_075, invdr);
	qfac = _mm256_mul_pd(qfac, _mm256_mul_pd(mii, mjj));
	qfac = _mm256_mul_pd(qfac, _mm256_mul_pd(invdr2, invdr2));

	const vcp_double_vec eii_xXdx = _mm256_mul_pd(eii_x, c_dx);
	const vcp_double_vec eii_yXdy = _mm256_mul_pd(eii_y, c_dy);
	const vcp_double_vec eii_zXdz = _mm256_mul_pd(eii_z, c_dz);
	vcp_double_vec costi = _mm256_add_pd(eii_xXdx, _mm256_add_pd(eii_yXdy, eii_zXdz));
	costi = _mm256_mul_pd(costi, invdr);

	const vcp_double_vec ejj_xXdx = _mm256_mul_pd(ejj_x, c_dx);
	const vcp_double_vec ejj_yXdy = _mm256_mul_pd(ejj_y, c_dy);
	const vcp_double_vec ejj_zXdz = _mm256_mul_pd(ejj_z, c_dz);
	vcp_double_vec costj = _mm256_add_pd(ejj_xXdx, _mm256_add_pd(ejj_yXdy, ejj_zXdz));
	costj = _mm256_mul_pd(costj, invdr);

	const vcp_double_vec cos2ti = _mm256_mul_pd(costi, costi);
	const vcp_double_vec cos2tj = _mm256_mul_pd(costj, costj);

	const vcp_double_vec eiiXejj_x = _mm256_mul_pd(eii_x, ejj_x);
	const vcp_double_vec eiiXejj_y = _mm256_mul_pd(eii_y, ejj_y);
	const vcp_double_vec eiiXejj_z = _mm256_mul_pd(eii_z, ejj_z);
	const vcp_double_vec cosgij = _mm256_add_pd(eiiXejj_x, _mm256_add_pd(eiiXejj_y, eiiXejj_z));

	vcp_double_vec term = _mm256_mul_pd(five, _mm256_mul_pd(costi, costj));
	term = _mm256_sub_pd(cosgij, term);

	/************
	 * Potential
	 ************/
	vcp_double_vec part1 = _mm256_mul_pd(five, _mm256_add_pd(cos2ti, cos2tj));
	vcp_double_vec part2 = _mm256_mul_pd(_15, _mm256_mul_pd(cos2ti, cos2tj));
	vcp_double_vec part3 = _mm256_mul_pd(two, _mm256_mul_pd(term, term));
	vcp_double_vec upot = _mm256_add_pd(part1, part2);
	upot = _mm256_sub_pd(_mm256_add_pd(one, part3), upot);
	upot = _mm256_mul_pd(qfac, upot);

	/**********
	 * Force
	 **********/
	const vcp_double_vec minus_partialRijInvdr = _mm256_mul_pd(five, _mm256_mul_pd(upot, invdr2));

	// partialTiInvdr & partialTjInvdr
	part1 = _mm256_mul_pd(qfac, _mm256_mul_pd(ten, invdr));
	part2 = _mm256_mul_pd(two, term);

	// partialTiInvdr only
	part3 = _mm256_mul_pd(three, _mm256_mul_pd(costi, cos2tj));
	vcp_double_vec part4 = _mm256_add_pd(costi, _mm256_add_pd(part3, _mm256_mul_pd(part2, costj)));
	const vcp_double_vec minus_partialTiInvdr = _mm256_mul_pd(part1, part4);

	// partialTjInvdr only
	part3 = _mm256_mul_pd(three, _mm256_mul_pd(costj, cos2ti));
	part4 = _mm256_add_pd(costj, _mm256_add_pd(part3, _mm256_mul_pd(part2, costi)));
	const vcp_double_vec minus_partialTjInvdr = _mm256_mul_pd(part1, part4);

	const vcp_double_vec partialGij = _mm256_mul_pd(qfac, _mm256_mul_pd(four, term));

	// fac
	part1 = _mm256_mul_pd(minus_partialTiInvdr, costi);
	part2 = _mm256_mul_pd(minus_partialTjInvdr, costj);
	part3 = _mm256_mul_pd(_mm256_add_pd(part1, part2), invdr);
	const vcp_double_vec fac = _mm256_sub_pd(minus_partialRijInvdr, part3);

	// Force components
	part1 = _mm256_mul_pd(fac, c_dx);
	part2 = _mm256_mul_pd(minus_partialTiInvdr, eii_x);
	part3 = _mm256_mul_pd(minus_partialTjInvdr, ejj_x);
	f_x = _mm256_add_pd(part1, _mm256_add_pd(part2, part3));

	part1 = _mm256_mul_pd(fac, c_dy);
	part2 = _mm256_mul_pd(minus_partialTiInvdr, eii_y);
	part3 = _mm256_mul_pd(minus_partialTjInvdr, ejj_y);
	f_y = _mm256_add_pd(part1, _mm256_add_pd(part2, part3));

	part1 = _mm256_mul_pd(fac, c_dz);
	part2 = _mm256_mul_pd(minus_partialTiInvdr, eii_z);
	part3 = _mm256_mul_pd(minus_partialTjInvdr, ejj_z);
	f_z = _mm256_add_pd(part1, _mm256_add_pd(part2, part3));

	const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
	// Check if we have to add the macroscopic values up for at least one of this pairs
	if (_mm256_movemask_pd(macroMask) > 0) {
		// do we have to mask "upot"...
		const vcp_double_vec upot_masked = _mm256_and_pd(upot, macroMask);
		sum_upotXpoles = _mm256_add_pd(sum_upotXpoles, upot_masked);

		const vcp_double_vec virial_x = _mm256_mul_pd(m_dx, f_x);
		const vcp_double_vec virial_y = _mm256_mul_pd(m_dy, f_y);
		const vcp_double_vec virial_z = _mm256_mul_pd(m_dz, f_z);

		const vcp_double_vec virial = _mm256_add_pd(_mm256_add_pd(virial_x, virial_y), virial_z);
		const vcp_double_vec virial_masked = _mm256_and_pd(virial, macroMask);
		sum_virial = _mm256_add_pd(sum_virial, virial_masked);
	}

	const vcp_double_vec eii_x_ejj_y = _mm256_mul_pd(eii_x, ejj_y);
	const vcp_double_vec eii_x_ejj_z = _mm256_mul_pd(eii_x, ejj_z);
	const vcp_double_vec eii_y_ejj_x = _mm256_mul_pd(eii_y, ejj_x);
	const vcp_double_vec eii_y_ejj_z = _mm256_mul_pd(eii_y, ejj_z);
	const vcp_double_vec eii_z_ejj_x = _mm256_mul_pd(eii_z, ejj_x);
	const vcp_double_vec eii_z_ejj_y = _mm256_mul_pd(eii_z, ejj_y);

	const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = _mm256_sub_pd(eii_x_ejj_y, eii_y_ejj_x);
	const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = _mm256_sub_pd(eii_y_ejj_z, eii_z_ejj_y);
	const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = _mm256_sub_pd(eii_z_ejj_x, eii_x_ejj_z);

	const vcp_double_vec partialGij_eiXej_x = _mm256_mul_pd(partialGij, eii_y_ejj_z_minus_eii_z_ejj_y);
	const vcp_double_vec partialGij_eiXej_y = _mm256_mul_pd(partialGij, eii_z_ejj_x_minus_eii_x_ejj_z);
	const vcp_double_vec partialGij_eiXej_z = _mm256_mul_pd(partialGij, eii_x_ejj_y_minus_eii_y_ejj_x);

	vcp_double_vec eXrij_x = _mm256_sub_pd(_mm256_mul_pd(eii_y, c_dz), _mm256_mul_pd(eii_z, c_dy));
	vcp_double_vec eXrij_y = _mm256_sub_pd(_mm256_mul_pd(eii_z, c_dx), _mm256_mul_pd(eii_x, c_dz));
	vcp_double_vec eXrij_z = _mm256_sub_pd(_mm256_mul_pd(eii_x, c_dy), _mm256_mul_pd(eii_y, c_dx));

	Mii_x = _mm256_sub_pd(_mm256_mul_pd(minus_partialTiInvdr, eXrij_x), partialGij_eiXej_x);
	Mii_y = _mm256_sub_pd(_mm256_mul_pd(minus_partialTiInvdr, eXrij_y), partialGij_eiXej_y);
	Mii_z = _mm256_sub_pd(_mm256_mul_pd(minus_partialTiInvdr, eXrij_z), partialGij_eiXej_z);

	eXrij_x = _mm256_sub_pd(_mm256_mul_pd(ejj_y, c_dz), _mm256_mul_pd(ejj_z, c_dy));
	eXrij_y = _mm256_sub_pd(_mm256_mul_pd(ejj_z, c_dx), _mm256_mul_pd(ejj_x, c_dz));
	eXrij_z = _mm256_sub_pd(_mm256_mul_pd(ejj_x, c_dy), _mm256_mul_pd(ejj_y, c_dx));

	Mjj_x = _mm256_add_pd(_mm256_mul_pd(minus_partialTjInvdr, eXrij_x), partialGij_eiXej_x);
	Mjj_y = _mm256_add_pd(_mm256_mul_pd(minus_partialTjInvdr, eXrij_y), partialGij_eiXej_y);
	Mjj_z = _mm256_add_pd(_mm256_mul_pd(minus_partialTjInvdr, eXrij_z), partialGij_eiXej_z);
}

#endif

#if VCP_VEC_TYPE==VCP_VEC_SSE3

inline void hSum_Add_Store( double * const mem_addr, const vcp_double_vec & a ) {
	_mm_store_sd(
			mem_addr,
			_mm_add_sd(_mm_hadd_pd(a, a),
			_mm_load_sd(mem_addr)));
}

#elif VCP_VEC_TYPE==VCP_VEC_AVX


inline void hSum_Add_Store( double * const mem_addr, const vcp_double_vec & a ) {
	const vcp_double_vec a_t1 = _mm256_permute2f128_pd(a, a, 0x1);
	const vcp_double_vec a_t2 = _mm256_hadd_pd(a, a_t1);
	const vcp_double_vec a_t3 = _mm256_hadd_pd(a_t2, a_t2);
	_mm256_maskstore_pd(
			mem_addr,
			memoryMask_first,
			_mm256_add_pd(
				a_t3,
				_mm256_maskload_pd(mem_addr, memoryMask_first)
			)

	);
}

#endif



template<class ForcePolicy>
#if VCP_VEC_TYPE==VCP_NOVEC
	unsigned long
#else
	vcp_double_vec
#endif
inline VectorizedCellProcessor::calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		double* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z
#if VCP_VEC_TYPE!=VCP_NOVEC
		, const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z
#endif
		) {

#if VCP_VEC_TYPE==VCP_NOVEC

	unsigned long compute_molecule = 0;

	for (size_t j = ForcePolicy :: InitJ(i_center_idx); j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];
		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		const signed long forceMask = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? (~0l) : 0l;
		compute_molecule |= forceMask;
		*(soa2_center_dist_lookup + j) = forceMask;
	}

	return compute_molecule;

#elif VCP_VEC_TYPE==VCP_VEC_SSE3

	vcp_double_vec compute_molecule = vcp_simd_zerov();

	// Iterate over centers of second cell
	size_t j = ForcePolicy :: InitJ(i_center_idx);
	for (; j < end_j; j+=2) {
		const vcp_double_vec m2_r_x = _mm_load_pd(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = _mm_load_pd(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = _mm_load_pd(soa2_m_r_z + j);

		const vcp_double_vec m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

		const vcp_double_vec m_dx2 = _mm_mul_pd(m_dx, m_dx);
		const vcp_double_vec m_dy2 = _mm_mul_pd(m_dy, m_dy);
		const vcp_double_vec m_dz2 = _mm_mul_pd(m_dz, m_dz);

		const vcp_double_vec m_r2 = _mm_add_pd(_mm_add_pd(m_dx2, m_dy2), m_dz2);

		const vcp_double_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD);
		_mm_store_pd(soa2_center_dist_lookup + j, forceMask);
		compute_molecule = _mm_or_pd(compute_molecule, forceMask);
	}

	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];

		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		// can we do this nicer?
		const signed long forceMask_l = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;
		// this casting via void* is required for gcc
		const void* forceMask_tmp = reinterpret_cast<const void*>(&forceMask_l);
		double forceMask = *reinterpret_cast<double const* /*const*/>(forceMask_tmp);

		*(soa2_center_dist_lookup + j) = forceMask;
		const vcp_double_vec forceMask_128 = _mm_set1_pd(forceMask);
		compute_molecule = _mm_or_pd(compute_molecule, forceMask_128);
	}

	return compute_molecule;

#elif VCP_VEC_TYPE==VCP_VEC_AVX

	vcp_double_vec compute_molecule = vcp_simd_zerov();

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	vcp_double_vec initJ_mask = ForcePolicy::InitJ_Mask(i_center_idx);
	for (; j < end_j; j+=4) {
		const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_m_r_z + j);

		const vcp_double_vec m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

		const vcp_double_vec m_dx2 = _mm256_mul_pd(m_dx, m_dx);
		const vcp_double_vec m_dy2 = _mm256_mul_pd(m_dy, m_dy);
		const vcp_double_vec m_dz2 = _mm256_mul_pd(m_dz, m_dz);

		const vcp_double_vec m_r2 = _mm256_add_pd(_mm256_add_pd(m_dx2, m_dy2), m_dz2);

		const vcp_double_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);
		_mm256_store_pd(soa2_center_dist_lookup + j, forceMask);
		compute_molecule = _mm256_or_pd(compute_molecule, forceMask);
	}

	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];

		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		// can we do this nicer?
		signed long forceMask_l;
		// DetectSingleCell() = false for SingleCellDistinctPolicy and CellPairPolicy, true for SingleCellPolicy
		if (ForcePolicy::DetectSingleCell()) {
			forceMask_l = (ForcePolicy::Condition(m_r2, cutoffRadiusSquare) && j > i_center_idx) ? ~0l : 0l;
		} else {
			forceMask_l = ForcePolicy::Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;
		}

//			this casting via void* is required for gcc
		void* forceMask_tmp = reinterpret_cast<void*>(&forceMask_l);
		double forceMask = *reinterpret_cast<double const* /*const*/>(forceMask_tmp);

		*(soa2_center_dist_lookup + j) = forceMask;
		const vcp_double_vec forceMask_256 = _mm256_set1_pd(forceMask);
		compute_molecule = _mm256_or_pd(compute_molecule, forceMask_256);
	}

	return compute_molecule;

#endif
}

template<class ForcePolicy, class MacroPolicy>
void VectorizedCellProcessor :: _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
#if VCP_VEC_TYPE==VCP_NOVEC
	// For the unvectorized version, we only have to iterate over all pairs of
	// LJ centers and apply the unvectorized loop body.
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

	for (size_t i = 0; i < soa1._mol_num; ++i) {
		// Computation of LJ interaction
		unsigned long compute_molecule_lj = calcDistLookup<ForcePolicy>(soa1, i, i_ljc_idx, soa2._ljc_num, _LJcutoffRadiusSquare,
				soa2._ljc_dist_lookup, soa2._ljc_m_r_x, soa2._ljc_m_r_y, soa2._ljc_m_r_z);
		unsigned long compute_molecule_charge  = calcDistLookup<ForcePolicy>(soa1, i, i_charge_idx , soa2._charges_num, _cutoffRadiusSquare,
				soa2._charges_dist_lookup, soa2._charges_m_r_x, soa2._charges_m_r_y, soa2._charges_m_r_z);
		unsigned long compute_molecule_dipole  = calcDistLookup<ForcePolicy>(soa1, i, i_dipole_idx , soa2._dipoles_num, _cutoffRadiusSquare,
				soa2._dipoles_dist_lookup, soa2._dipoles_m_r_x, soa2._dipoles_m_r_y, soa2._dipoles_m_r_z);
		unsigned long compute_molecule_quadrupole  = calcDistLookup<ForcePolicy>(soa1, i, i_quadrupole_idx, soa2._quadrupoles_num, _cutoffRadiusSquare,
				soa2._quadrupoles_dist_lookup, soa2._quadrupoles_m_r_x, soa2._quadrupoles_m_r_y, soa2._quadrupoles_m_r_z);

		// Computation of LJ-site interactions

		if (!compute_molecule_lj) {
			i_ljc_idx += soa1._mol_ljc_num[i];
		}
		else {
			for (int local_i = 0; local_i < soa1._mol_ljc_num[i]; local_i++) {
				for (size_t j = ForcePolicy :: InitJ(i_ljc_idx); j < soa2._ljc_num; ++j) {
					_loopBodyNovecLJ<MacroPolicy>(soa1, i_ljc_idx, soa2, j, soa2._ljc_dist_lookup + j);
				}
				i_ljc_idx++;
			}
		}

		// Computation of site interactions with charges

		if (!compute_molecule_charge) {
			i_charge_idx += soa1._mol_charges_num[i];
			i_dipole_charge_idx += soa1._mol_dipoles_num[i];
			i_quadrupole_charge_idx += soa1._mol_quadrupoles_num[i];
		}
		else {
			for (int local_i = 0; local_i < soa1._mol_charges_num[i]; local_i++) {
				for (size_t j = ForcePolicy :: InitJ(i_charge_idx + local_i); j < soa2._charges_num; ++j) {
					_loopBodyNovecCharges<MacroPolicy>(soa1, i_charge_idx + local_i, soa2, j, soa2._charges_dist_lookup + j);
				}
			}

			for (int local_i = 0; local_i < soa1._mol_dipoles_num[i]; local_i++)
			{
				for (size_t j = ForcePolicy :: InitJ(i_charge_idx); j < soa2._charges_num; ++j) {
					_loopBodyNovecChargesDipoles<MacroPolicy>(soa2, j, soa1, i_dipole_charge_idx, soa2._charges_dist_lookup + j, true); // true invertiert macroscopicValueCondition, false nicht.
				}
				i_dipole_charge_idx++;
			}

			for (int local_i = 0; local_i < soa1._mol_quadrupoles_num[i]; local_i++)
			{
				for (size_t j = ForcePolicy :: InitJ(i_charge_idx); j < soa2._charges_num; ++j) {
					_loopBodyNovecChargesQuadrupoles<MacroPolicy>(soa2, j, soa1, i_quadrupole_charge_idx, soa2._charges_dist_lookup + j, true); // true invertiert macroscopicValueCondition, false nicht.
				}
				i_quadrupole_charge_idx++;
			}

			i_charge_idx += soa1._mol_charges_num[i];
		}

		// Computation of site interactions with dipoles

		if (!compute_molecule_dipole) {
			i_dipole_idx += soa1._mol_dipoles_num[i];
			i_charge_dipole_idx += soa1._mol_charges_num[i];
			i_quadrupole_dipole_idx += soa1._mol_quadrupoles_num[i];
		}
		else {
			for (int local_i = 0; local_i < soa1._mol_dipoles_num[i]; local_i++) {
				for (size_t j = ForcePolicy :: InitJ(i_dipole_idx + local_i); j < soa2._dipoles_num; ++j) {
					_loopBodyNovecDipoles<MacroPolicy>(soa1, i_dipole_idx + local_i, soa2, j, soa2._dipoles_dist_lookup + j);
				}
			}

			for (int local_i = 0; local_i < soa1._mol_charges_num[i]; local_i++) {
				for (size_t j = ForcePolicy::InitJ(i_dipole_idx); j < soa2._dipoles_num; ++j) {
					_loopBodyNovecChargesDipoles<MacroPolicy>(soa1, i_charge_dipole_idx, soa2, j, soa2._dipoles_dist_lookup + j, false);
				}
				i_charge_dipole_idx++;
			}

			for (int local_i = 0; local_i < soa1._mol_quadrupoles_num[i]; local_i++) {
				for (size_t j = ForcePolicy::InitJ(i_dipole_idx); j < soa2._dipoles_num; ++j) {
					_loopBodyNovecDipolesQuadrupoles<MacroPolicy>(soa2, j, soa1, i_quadrupole_dipole_idx, soa2._dipoles_dist_lookup + j, true);
				}
				i_quadrupole_dipole_idx++;

			}

			i_dipole_idx += soa1._mol_dipoles_num[i];
		}


		// Computation of site interactions with quadrupoles

		if (!compute_molecule_quadrupole) {
			i_quadrupole_idx += soa1._mol_quadrupoles_num[i];
			i_charge_quadrupole_idx += soa1._mol_charges_num[i];
			i_dipole_quadrupole_idx += soa1._mol_dipoles_num[i];
		}
		else {
			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1._mol_quadrupoles_num[i]; local_i++) {
				for (size_t j = ForcePolicy :: InitJ(i_quadrupole_idx + local_i); j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecQuadrupoles<MacroPolicy>(soa1, i_quadrupole_idx + local_i , soa2, j, soa2._quadrupoles_dist_lookup + j);
				}
			}

			for (int local_i = 0; local_i < soa1._mol_charges_num[i]; local_i++) {
				for (size_t j = ForcePolicy :: InitJ(i_quadrupole_idx); j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecChargesQuadrupoles<MacroPolicy>(soa1, i_charge_quadrupole_idx, soa2, j, soa2._quadrupoles_dist_lookup + j, false);
				}
				i_charge_quadrupole_idx++;
			}

			for (int local_i = 0; local_i < soa1._mol_dipoles_num[i]; local_i++) {
				for (size_t j = ForcePolicy :: InitJ(i_quadrupole_idx); j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecDipolesQuadrupoles<MacroPolicy>(soa1, i_dipole_quadrupole_idx, soa2, j, soa2._quadrupoles_dist_lookup + j, false);
				}

				i_dipole_quadrupole_idx++;
			}

			i_quadrupole_idx += soa1._mol_quadrupoles_num[i];

		}
	}


#else
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
	const size_t * const soa2_ljc_id = soa2._ljc_id;

	double* const soa2_ljc_dist_lookup = soa2._ljc_dist_lookup;

	// Pointer for charges
	const double * const soa1_charges_r_x = soa1._charges_r_x;
	const double * const soa1_charges_r_y = soa1._charges_r_y;
	const double * const soa1_charges_r_z = soa1._charges_r_z;
	double * const soa1_charges_f_x = soa1._charges_f_x;
	double * const soa1_charges_f_y = soa1._charges_f_y;
	double * const soa1_charges_f_z = soa1._charges_f_z;
	const double * const soa1_charges_q = soa1._charges_q;
	const int * const soa1_mol_charges_num = soa1._mol_charges_num;

	const double * const soa2_charges_m_r_x = soa2._charges_m_r_x;
	const double * const soa2_charges_m_r_y = soa2._charges_m_r_y;
	const double * const soa2_charges_m_r_z = soa2._charges_m_r_z;
	const double * const soa2_charges_r_x = soa2._charges_r_x;
	const double * const soa2_charges_r_y = soa2._charges_r_y;
	const double * const soa2_charges_r_z = soa2._charges_r_z;
	double * const soa2_charges_f_x = soa2._charges_f_x;
	double * const soa2_charges_f_y = soa2._charges_f_y;
	double * const soa2_charges_f_z = soa2._charges_f_z;
	const double * const soa2_charges_q = soa2._charges_q;

	double* const soa2_charges_dist_lookup = soa2._charges_dist_lookup;

	// Pointer for dipoles
	const double * const soa1_dipoles_r_x = soa1._dipoles_r_x;
	const double * const soa1_dipoles_r_y = soa1._dipoles_r_y;
	const double * const soa1_dipoles_r_z = soa1._dipoles_r_z;
	double * const soa1_dipoles_f_x = soa1._dipoles_f_x;
	double * const soa1_dipoles_f_y = soa1._dipoles_f_y;
	double * const soa1_dipoles_f_z = soa1._dipoles_f_z;
	const double * const soa1_dipoles_p = soa1._dipoles_p;
	const double * const soa1_dipoles_e_x = soa1._dipoles_e_x;
	const double * const soa1_dipoles_e_y = soa1._dipoles_e_y;
	const double * const soa1_dipoles_e_z = soa1._dipoles_e_z;
	double * const soa1_dipoles_M_x = soa1._dipoles_M_x;
	double * const soa1_dipoles_M_y = soa1._dipoles_M_y;
	double * const soa1_dipoles_M_z = soa1._dipoles_M_z;
	const int * const soa1_mol_dipoles_num = soa1._mol_dipoles_num;

	const double * const soa2_dipoles_m_r_x = soa2._dipoles_m_r_x;
	const double * const soa2_dipoles_m_r_y = soa2._dipoles_m_r_y;
	const double * const soa2_dipoles_m_r_z = soa2._dipoles_m_r_z;
	const double * const soa2_dipoles_r_x = soa2._dipoles_r_x;
	const double * const soa2_dipoles_r_y = soa2._dipoles_r_y;
	const double * const soa2_dipoles_r_z = soa2._dipoles_r_z;
	double * const soa2_dipoles_f_x = soa2._dipoles_f_x;
	double * const soa2_dipoles_f_y = soa2._dipoles_f_y;
	double * const soa2_dipoles_f_z = soa2._dipoles_f_z;
	const double * const soa2_dipoles_p = soa2._dipoles_p;
	const double * const soa2_dipoles_e_x = soa2._dipoles_e_x;
	const double * const soa2_dipoles_e_y = soa2._dipoles_e_y;
	const double * const soa2_dipoles_e_z = soa2._dipoles_e_z;
	double * const soa2_dipoles_M_x = soa2._dipoles_M_x;
	double * const soa2_dipoles_M_y = soa2._dipoles_M_y;
	double * const soa2_dipoles_M_z = soa2._dipoles_M_z;

	double* const soa2_dipoles_dist_lookup = soa2._dipoles_dist_lookup;

	// Pointer for quadrupoles
	const double * const soa1_quadrupoles_r_x = soa1._quadrupoles_r_x;
	const double * const soa1_quadrupoles_r_y = soa1._quadrupoles_r_y;
	const double * const soa1_quadrupoles_r_z = soa1._quadrupoles_r_z;
	double * const soa1_quadrupoles_f_x = soa1._quadrupoles_f_x;
	double * const soa1_quadrupoles_f_y = soa1._quadrupoles_f_y;
	double * const soa1_quadrupoles_f_z = soa1._quadrupoles_f_z;
	const double * const soa1_quadrupoles_m = soa1._quadrupoles_m;
	const double * const soa1_quadrupoles_e_x = soa1._quadrupoles_e_x;
	const double * const soa1_quadrupoles_e_y = soa1._quadrupoles_e_y;
	const double * const soa1_quadrupoles_e_z = soa1._quadrupoles_e_z;
	double * const soa1_quadrupoles_M_x = soa1._quadrupoles_M_x;
	double * const soa1_quadrupoles_M_y = soa1._quadrupoles_M_y;
	double * const soa1_quadrupoles_M_z = soa1._quadrupoles_M_z;
	const int * const soa1_mol_quadrupoles_num = soa1._mol_quadrupoles_num;

	const double * const soa2_quadrupoles_m_r_x = soa2._quadrupoles_m_r_x;
	const double * const soa2_quadrupoles_m_r_y = soa2._quadrupoles_m_r_y;
	const double * const soa2_quadrupoles_m_r_z = soa2._quadrupoles_m_r_z;
	const double * const soa2_quadrupoles_r_x = soa2._quadrupoles_r_x;
	const double * const soa2_quadrupoles_r_y = soa2._quadrupoles_r_y;
	const double * const soa2_quadrupoles_r_z = soa2._quadrupoles_r_z;
	double * const soa2_quadrupoles_f_x = soa2._quadrupoles_f_x;
	double * const soa2_quadrupoles_f_y = soa2._quadrupoles_f_y;
	double * const soa2_quadrupoles_f_z = soa2._quadrupoles_f_z;
	const double * const soa2_quadrupoles_m = soa2._quadrupoles_m;
	const double * const soa2_quadrupoles_e_x = soa2._quadrupoles_e_x;
	const double * const soa2_quadrupoles_e_y = soa2._quadrupoles_e_y;
	const double * const soa2_quadrupoles_e_z = soa2._quadrupoles_e_z;
	double * const soa2_quadrupoles_M_x = soa2._quadrupoles_M_x;
	double * const soa2_quadrupoles_M_y = soa2._quadrupoles_M_y;
	double * const soa2_quadrupoles_M_z = soa2._quadrupoles_M_z;

	double* const soa2_quadrupoles_dist_lookup = soa2._quadrupoles_dist_lookup;


#if VCP_VEC_TYPE==VCP_VEC_SSE3

	const vcp_double_vec ones = _mm_castsi128_pd(_mm_set_epi32(~0, ~0, ~0, ~0));

	vcp_double_vec sum_upot6lj = vcp_simd_zerov();
	vcp_double_vec sum_upotXpoles = vcp_simd_zerov();
	vcp_double_vec sum_virial = vcp_simd_zerov();
	vcp_double_vec sum_myRF = vcp_simd_zerov();

	const vcp_double_vec rc2 = _mm_set1_pd(_LJcutoffRadiusSquare);
	const vcp_double_vec cutoffRadiusSquare = _mm_set1_pd(_cutoffRadiusSquare);
	const vcp_double_vec epsRFInvrc3 = _mm_loaddup_pd(&_epsRFInvrc3);

	/*
	 *  "var" & (~1):
	 *  Making sure that result is a multiple of 2. If "var" can not be divided by 2
	 *  the result is: "var" decremented by one.
	 */
	const size_t end_ljc_j = soa2._ljc_num & (~1);
	const size_t end_charges_j = soa2._charges_num & (~1);
	const size_t end_dipoles_j = soa2._dipoles_num & (~1);
	const size_t end_quadrupoles_j = soa2._quadrupoles_num & (~1);

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

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._mol_num; ++i) {
		const vcp_double_vec m1_r_x = _mm_loaddup_pd(soa1_mol_pos_x + i);
		const vcp_double_vec m1_r_y = _mm_loaddup_pd(soa1_mol_pos_y + i);
		const vcp_double_vec m1_r_z = _mm_loaddup_pd(soa1_mol_pos_z + i);


		const vcp_double_vec compute_molecule_ljc = calcDistLookup<ForcePolicy>(soa1, i, i_ljc_idx, soa2._ljc_num, _LJcutoffRadiusSquare,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_double_vec compute_molecule_charges = calcDistLookup<ForcePolicy>(soa1, i, i_charge_idx, soa2._charges_num, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_double_vec compute_molecule_dipoles = calcDistLookup<ForcePolicy>(soa1, i, i_dipole_idx, soa2._dipoles_num, _cutoffRadiusSquare,
				soa2_dipoles_dist_lookup, soa2_dipoles_m_r_x, soa2_dipoles_m_r_y, soa2_dipoles_m_r_z,
				cutoffRadiusSquare,	end_dipoles_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_double_vec compute_molecule_quadrupoles = calcDistLookup<ForcePolicy>(soa1, i, i_quadrupole_idx, soa2._quadrupoles_num, _cutoffRadiusSquare,
				soa2_quadrupoles_dist_lookup, soa2_quadrupoles_m_r_x, soa2_quadrupoles_m_r_y, soa2_quadrupoles_m_r_z,
				cutoffRadiusSquare, end_quadrupoles_j, m1_r_x, m1_r_y, m1_r_z);

		if (!_mm_movemask_pd(compute_molecule_ljc)) {
			i_ljc_idx += soa1_mol_ljc_num[i];
		}
		else {

			// LJ force computation
			for (int local_i = 0; local_i < soa1_mol_ljc_num[i]; local_i++) {
				vcp_double_vec sum_fx1 = vcp_simd_zerov();
				vcp_double_vec sum_fy1 = vcp_simd_zerov();
				vcp_double_vec sum_fz1 = vcp_simd_zerov();
				const vcp_double_vec c_r_x1 = _mm_loaddup_pd(soa1_ljc_r_x + i_ljc_idx);
				const vcp_double_vec c_r_y1 = _mm_loaddup_pd(soa1_ljc_r_y + i_ljc_idx);
				const vcp_double_vec c_r_z1 = _mm_loaddup_pd(soa1_ljc_r_z + i_ljc_idx);
				// Iterate over each pair of centers in the second cell.
				size_t j = ForcePolicy::InitJ(i_ljc_idx);
				for (; j < end_ljc_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_ljc_dist_lookup + j);
					// Only go on if at least 1 of the forces has to be calculated.
					if (_mm_movemask_pd(forceMask) > 0) {

						const vcp_double_vec c_r_x2 = _mm_load_pd(soa2_ljc_r_x + j);
						const vcp_double_vec c_r_y2 = _mm_load_pd(soa2_ljc_r_y + j);
						const vcp_double_vec c_r_z2 = _mm_load_pd(soa2_ljc_r_z + j);

						const size_t id_j0 = soa2_ljc_id[j];
						const size_t id_j1 = soa2_ljc_id[j + 1];
						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						const vcp_double_vec e1s1 = _mm_load_pd(_eps_sig[id_i] + 2 * id_j0);
						const vcp_double_vec e2s2 = _mm_load_pd(_eps_sig[id_i] + 2 * id_j1);

						const vcp_double_vec m_r_x2 = _mm_load_pd(soa2_ljc_m_r_x + j);
						const vcp_double_vec m_r_y2 = _mm_load_pd(soa2_ljc_m_r_y + j);
						const vcp_double_vec m_r_z2 = _mm_load_pd(soa2_ljc_m_r_z + j);

						vcp_double_vec fx, fy, fz;

						_loopBodyLJ<MacroPolicy>(
							m1_r_x, m1_r_y, m1_r_z, c_r_x1, c_r_y1, c_r_z1,
							m_r_x2, m_r_y2, m_r_z2, c_r_x2, c_r_y2, c_r_z2,
							fx, fy, fz,
							sum_upot6lj, sum_virial,
							forceMask,
							e1s1, e2s2,
							id_j0, id_j1, id_i);

						const vcp_double_vec old_fx2 = _mm_load_pd(soa2_ljc_f_x + j);
						const vcp_double_vec new_fx2 = _mm_sub_pd(old_fx2, fx);
						_mm_store_pd(soa2_ljc_f_x + j, new_fx2);
						const vcp_double_vec old_fy2 = _mm_load_pd(soa2_ljc_f_y + j);
						const vcp_double_vec new_fy2 = _mm_sub_pd(old_fy2, fy);
						_mm_store_pd(soa2_ljc_f_y + j, new_fy2);
						const vcp_double_vec old_fz2 = _mm_load_pd(soa2_ljc_f_z + j);
						const vcp_double_vec new_fz2 = _mm_sub_pd(old_fz2, fz);
						_mm_store_pd(soa2_ljc_f_z + j, new_fz2);
						sum_fx1 = _mm_add_pd(sum_fx1, fx);
						sum_fy1 = _mm_add_pd(sum_fy1, fy);
						sum_fz1 = _mm_add_pd(sum_fz1, fz);
					}
				}

				hSum_Add_Store(soa1_ljc_f_x + i_ljc_idx, sum_fx1);
				hSum_Add_Store(soa1_ljc_f_y + i_ljc_idx, sum_fy1);
				hSum_Add_Store(soa1_ljc_f_z + i_ljc_idx, sum_fz1);

				// Unvectorized calculation for leftover pairs.
				switch (soa2._ljc_num & 1) {
					case 1: {
						_loopBodyNovecLJ<MacroPolicy>(soa1, i_ljc_idx, soa2, end_ljc_j, soa2_ljc_dist_lookup + j);
					}
					break;
				}

				i_ljc_idx++;
			}
		}

		// Computation of site interactions with charges

		if (!_mm_movemask_pd(compute_molecule_charges)) {
			i_charge_idx += soa1_mol_charges_num[i];
			i_dipole_charge_idx += soa1_mol_dipoles_num[i];
			i_quadrupole_charge_idx += soa1_mol_quadrupoles_num[i];
		}
		else {
			// Computation of charge-charge interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++) {

				const vcp_double_vec q1 = _mm_loaddup_pd(soa1_charges_q + i_charge_idx + local_i);
				const vcp_double_vec r1_x = _mm_loaddup_pd(soa1_charges_r_x + i_charge_idx + local_i);
				const vcp_double_vec r1_y = _mm_loaddup_pd(soa1_charges_r_y + i_charge_idx + local_i);
				const vcp_double_vec r1_z = _mm_loaddup_pd(soa1_charges_r_z + i_charge_idx + local_i);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_charge_idx + local_i);
				for (; j < end_charges_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {
						const vcp_double_vec q2 = _mm_load_pd(soa2_charges_q + j);

						const vcp_double_vec r2_x = _mm_load_pd(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = _mm_load_pd(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = _mm_load_pd(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z;
						_loopBodyCharge<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								sum_upotXpoles, sum_virial,
								forceMask);



						sum_f1_x = _mm_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_add_pd(sum_f1_z, f_z);

						vcp_double_vec f2_x = _mm_load_pd(soa2_charges_f_x + j);
						vcp_double_vec f2_y = _mm_load_pd(soa2_charges_f_y + j);
						vcp_double_vec f2_z = _mm_load_pd(soa2_charges_f_z + j);

						f2_x = _mm_sub_pd(f2_x, f_x);
						f2_y = _mm_sub_pd(f2_y, f_y);
						f2_z = _mm_sub_pd(f2_z, f_z);

						_mm_store_pd(soa2_charges_f_x + j, f2_x);
						_mm_store_pd(soa2_charges_f_y + j, f2_y);
						_mm_store_pd(soa2_charges_f_z + j, f2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx + local_i, sum_f1_z);


				// End iteration over centers with possible left over center
				for (; j < soa2._charges_num; ++j) {
					_loopBodyNovecCharges<MacroPolicy>(soa1, i_charge_idx + local_i, soa2, j, soa2_charges_dist_lookup + j);
				}

			}

			// Computation of dipole-charge interactions

			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++)
			{
				const vcp_double_vec p = _mm_loaddup_pd(soa1_dipoles_p + i_dipole_charge_idx);
				const vcp_double_vec e_x = _mm_loaddup_pd(soa1_dipoles_e_x + i_dipole_charge_idx);
				const vcp_double_vec e_y = _mm_loaddup_pd(soa1_dipoles_e_y + i_dipole_charge_idx);
				const vcp_double_vec e_z = _mm_loaddup_pd(soa1_dipoles_e_z + i_dipole_charge_idx);
				const vcp_double_vec r1_x = _mm_loaddup_pd(soa1_dipoles_r_x + i_dipole_charge_idx);
				const vcp_double_vec r1_y = _mm_loaddup_pd(soa1_dipoles_r_y + i_dipole_charge_idx);
				const vcp_double_vec r1_z = _mm_loaddup_pd(soa1_dipoles_r_z + i_dipole_charge_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M_x = vcp_simd_zerov();
				vcp_double_vec sum_M_y = vcp_simd_zerov();
				vcp_double_vec sum_M_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {

						const vcp_double_vec q = _mm_load_pd(soa2_charges_q + j);

						const vcp_double_vec r2_x = _mm_load_pd(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = _mm_load_pd(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = _mm_load_pd(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeDipole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = _mm_sub_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_sub_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_sub_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_charges_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_charges_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_charges_f_z + j);

						sum_f2_x = _mm_add_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_add_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_add_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_charges_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_charges_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_charges_f_z + j, sum_f2_z);

						// Store torque

						sum_M_x = _mm_add_pd(sum_M_x, M_x);
						sum_M_y = _mm_add_pd(sum_M_y, M_y);
						sum_M_z = _mm_add_pd(sum_M_z, M_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_charge_idx, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_charge_idx, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_charge_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_charge_idx, sum_M_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_charge_idx, sum_M_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_charge_idx, sum_M_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._charges_num; ++j) {
					_loopBodyNovecChargesDipoles<MacroPolicy>(soa2, j, soa1, i_dipole_charge_idx, soa2_charges_dist_lookup + j, true);
				}

				i_dipole_charge_idx++;
			}

			// Computation of quadrupole-charge interactions

			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++)
			{
				const vcp_double_vec m = _mm_loaddup_pd(soa1_quadrupoles_m + i_quadrupole_charge_idx);
				const vcp_double_vec e_x = _mm_loaddup_pd(soa1_quadrupoles_e_x + i_quadrupole_charge_idx);
				const vcp_double_vec e_y = _mm_loaddup_pd(soa1_quadrupoles_e_y + i_quadrupole_charge_idx);
				const vcp_double_vec e_z = _mm_loaddup_pd(soa1_quadrupoles_e_z + i_quadrupole_charge_idx);
				const vcp_double_vec r1_x = _mm_loaddup_pd(soa1_quadrupoles_r_x + i_quadrupole_charge_idx);
				const vcp_double_vec r1_y = _mm_loaddup_pd(soa1_quadrupoles_r_y + i_quadrupole_charge_idx);
				const vcp_double_vec r1_z = _mm_loaddup_pd(soa1_quadrupoles_r_z + i_quadrupole_charge_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {

						const vcp_double_vec q = _mm_load_pd(soa2_charges_q + j);

						const vcp_double_vec r2_x = _mm_load_pd(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = _mm_load_pd(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = _mm_load_pd(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeQuadrupole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = _mm_sub_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_sub_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_sub_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_charges_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_charges_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_charges_f_z + j);

						sum_f2_x = _mm_add_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_add_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_add_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_charges_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_charges_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_charges_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm_add_pd(sum_M1_x, M_x);
						sum_M1_y = _mm_add_pd(sum_M1_y, M_y);
						sum_M1_z = _mm_add_pd(sum_M1_z, M_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_charge_idx, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_charge_idx, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_charge_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_charge_idx, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_charge_idx, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_charge_idx, sum_M1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._charges_num; ++j) {
					_loopBodyNovecChargesQuadrupoles<MacroPolicy>(soa2, j, soa1, i_quadrupole_charge_idx, soa2_charges_dist_lookup + j, true);
				}
				i_quadrupole_charge_idx++;
			}

			i_charge_idx += soa1_mol_charges_num[i];
		}

		// Computation of site interactions with dipoles

		// Continue with next molecule if no force has to be calculated
		if (!_mm_movemask_pd(compute_molecule_dipoles)) {
			i_dipole_idx += soa1_mol_dipoles_num[i];
			i_charge_dipole_idx += soa1_mol_charges_num[i];
			i_quadrupole_dipole_idx += soa1_mol_quadrupoles_num[i];
		}
		else {
			// Computation of dipole-dipole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++) {

				const vcp_double_vec p1 = _mm_loaddup_pd(soa1_dipoles_p + i_dipole_idx + local_i);
				const vcp_double_vec e1_x = _mm_loaddup_pd(soa1_dipoles_e_x + i_dipole_idx + local_i);
				const vcp_double_vec e1_y = _mm_loaddup_pd(soa1_dipoles_e_y + i_dipole_idx + local_i);
				const vcp_double_vec e1_z = _mm_loaddup_pd(soa1_dipoles_e_z + i_dipole_idx + local_i);
				const vcp_double_vec r1_x = _mm_loaddup_pd(soa1_dipoles_r_x + i_dipole_idx + local_i);
				const vcp_double_vec r1_y = _mm_loaddup_pd(soa1_dipoles_r_y + i_dipole_idx + local_i);
				const vcp_double_vec r1_z = _mm_loaddup_pd(soa1_dipoles_r_z + i_dipole_idx + local_i);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_dipole_idx + local_i);
				for (; j < end_dipoles_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {
						const vcp_double_vec p2 = _mm_load_pd(soa2_dipoles_p + j);
						const vcp_double_vec e2_x = _mm_load_pd(soa2_dipoles_e_x + j);
						const vcp_double_vec e2_y = _mm_load_pd(soa2_dipoles_e_y + j);
						const vcp_double_vec e2_z = _mm_load_pd(soa2_dipoles_e_z + j);
						const vcp_double_vec r2_x = _mm_load_pd(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = _mm_load_pd(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = _mm_load_pd(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipole<MacroPolicy>(
							m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, p1,
							m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p2,
							f_x, f_y, f_z,
							M1_x, M1_y, M1_z,
							M2_x, M2_y, M2_z,
							sum_upotXpoles, sum_virial, sum_myRF,
							forceMask,
							epsRFInvrc3);

						// Store forces

						sum_f1_x = _mm_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_dipoles_f_z + j);

						sum_f2_x = _mm_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_sub_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_dipoles_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_dipoles_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm_load_pd(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm_load_pd(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm_load_pd(soa2_dipoles_M_z + j);

						sum_M2_x = _mm_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm_add_pd(sum_M2_z, M2_z);

						_mm_store_pd(soa2_dipoles_M_x + j, sum_M2_x);
						_mm_store_pd(soa2_dipoles_M_y + j, sum_M2_y);
						_mm_store_pd(soa2_dipoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_idx + local_i, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_idx + local_i, sum_M1_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_idx + local_i, sum_M1_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_idx + local_i, sum_M1_z);


				// End iteration over centers with possible left over center
				for (; j < soa2._dipoles_num; ++j) {
					_loopBodyNovecDipoles<MacroPolicy>(soa1, i_dipole_idx + local_i, soa2, j, soa2_dipoles_dist_lookup + j);
				}

			}

			// Computation of charge-dipole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{

				const vcp_double_vec q = _mm_loaddup_pd(soa1_charges_q + i_charge_dipole_idx);
				const vcp_double_vec r1_x = _mm_loaddup_pd(soa1_charges_r_x + i_charge_dipole_idx);
				const vcp_double_vec r1_y = _mm_loaddup_pd(soa1_charges_r_y + i_charge_dipole_idx);
				const vcp_double_vec r1_z = _mm_loaddup_pd(soa1_charges_r_z + i_charge_dipole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_dipole_idx);
				for (; j < end_dipoles_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {

						const vcp_double_vec p = _mm_load_pd(soa2_dipoles_p + j);

						const vcp_double_vec e_x = _mm_load_pd(soa2_dipoles_e_x + j);
						const vcp_double_vec e_y = _mm_load_pd(soa2_dipoles_e_y + j);
						const vcp_double_vec e_z = _mm_load_pd(soa2_dipoles_e_z + j);

						const vcp_double_vec r2_x = _mm_load_pd(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = _mm_load_pd(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = _mm_load_pd(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeDipole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = _mm_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_dipoles_f_z + j);

						sum_f2_x = _mm_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_sub_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_dipoles_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_dipoles_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						vcp_double_vec sum_M2_x = _mm_load_pd(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm_load_pd(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm_load_pd(soa2_dipoles_M_z + j);

						sum_M2_x = _mm_add_pd(sum_M2_x, M_x);
						sum_M2_y = _mm_add_pd(sum_M2_y, M_y);
						sum_M2_z = _mm_add_pd(sum_M2_z, M_z);

						_mm_store_pd(soa2_dipoles_M_x + j, sum_M2_x);
						_mm_store_pd(soa2_dipoles_M_y + j, sum_M2_y);
						_mm_store_pd(soa2_dipoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_dipole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_dipole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_dipole_idx, sum_f1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._dipoles_num; ++j) {
					_loopBodyNovecChargesDipoles<MacroPolicy>(soa1, i_charge_dipole_idx, soa2, j, soa2_dipoles_dist_lookup + j, false);
				}

				i_charge_dipole_idx++;
			}

			// Computation of quadrupole-dipole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++) {

				const vcp_double_vec m = _mm_loaddup_pd(soa1_quadrupoles_m + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_x = _mm_loaddup_pd(soa1_quadrupoles_e_x + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_y = _mm_loaddup_pd(soa1_quadrupoles_e_y + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_z = _mm_loaddup_pd(soa1_quadrupoles_e_z + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_x = _mm_loaddup_pd(soa1_quadrupoles_r_x + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_y = _mm_loaddup_pd(soa1_quadrupoles_r_y + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_z = _mm_loaddup_pd(soa1_quadrupoles_r_z + i_quadrupole_dipole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_dipole_idx);
				for (; j < end_dipoles_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {
						const vcp_double_vec p = _mm_load_pd(soa2_dipoles_p + j);
						const vcp_double_vec e2_x = _mm_load_pd(soa2_dipoles_e_x + j);
						const vcp_double_vec e2_y = _mm_load_pd(soa2_dipoles_e_y + j);
						const vcp_double_vec e2_z = _mm_load_pd(soa2_dipoles_e_z + j);
						const vcp_double_vec r2_x = _mm_load_pd(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = _mm_load_pd(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = _mm_load_pd(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipoleQuadrupole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, m,
								f_x, f_y, f_z, M2_x, M2_y, M2_z, M1_x, M1_y, M1_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = _mm_sub_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_sub_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_sub_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_dipoles_f_z + j);

						sum_f2_x = _mm_add_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_add_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_add_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_dipoles_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_dipoles_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm_load_pd(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm_load_pd(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm_load_pd(soa2_dipoles_M_z + j);

						sum_M2_x = _mm_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm_add_pd(sum_M2_z, M2_z);

						_mm_store_pd(soa2_dipoles_M_x + j, sum_M2_x);
						_mm_store_pd(soa2_dipoles_M_y + j, sum_M2_y);
						_mm_store_pd(soa2_dipoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_dipole_idx, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_dipole_idx, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_dipole_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_dipole_idx, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_dipole_idx, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_dipole_idx, sum_M1_z);


				// End iteration over centers with possible left over center
				for (; j < soa2._dipoles_num; ++j) {
					_loopBodyNovecDipolesQuadrupoles<MacroPolicy>(soa2, j, soa1, i_quadrupole_dipole_idx, soa2_dipoles_dist_lookup + j, true);
				}

				i_quadrupole_dipole_idx++;

			}

			i_dipole_idx += soa1_mol_dipoles_num[i];
		}

		// Computation of site interactions with quadrupoles

		if (!_mm_movemask_pd(compute_molecule_quadrupoles)) {
			i_quadrupole_idx += soa1_mol_quadrupoles_num[i];
			i_charge_quadrupole_idx += soa1_mol_charges_num[i];
			i_dipole_quadrupole_idx += soa1_mol_dipoles_num[i];
		}
		else {
			// Computation of quadrupole-quadrupole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++)
			{
				const vcp_double_vec mii = _mm_loaddup_pd(soa1_quadrupoles_m + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_x = _mm_loaddup_pd(soa1_quadrupoles_e_x + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_y = _mm_loaddup_pd(soa1_quadrupoles_e_y + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_z = _mm_loaddup_pd(soa1_quadrupoles_e_z + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_x = _mm_loaddup_pd(soa1_quadrupoles_r_x + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_y = _mm_loaddup_pd(soa1_quadrupoles_r_y + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_z = _mm_loaddup_pd(soa1_quadrupoles_r_z + i_quadrupole_idx + local_i);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_quadrupole_idx + local_i);
				for (; j < end_quadrupoles_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {

						const vcp_double_vec mjj = _mm_load_pd(soa2_quadrupoles_m + j);
						const vcp_double_vec ejj_x = _mm_load_pd(soa2_quadrupoles_e_x + j);
						const vcp_double_vec ejj_y = _mm_load_pd(soa2_quadrupoles_e_y + j);
						const vcp_double_vec ejj_z = _mm_load_pd(soa2_quadrupoles_e_z + j);
						const vcp_double_vec rjj_x = _mm_load_pd(soa2_quadrupoles_r_x + j);
						const vcp_double_vec rjj_y = _mm_load_pd(soa2_quadrupoles_r_y + j);
						const vcp_double_vec rjj_z = _mm_load_pd(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_quadrupoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyQudarupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, mii,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, mjj,
								f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								forceMask);

						// Store forces

						sum_f1_x = _mm_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_quadrupoles_f_z + j);

						sum_f2_x = _mm_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_sub_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_quadrupoles_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_quadrupoles_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm_load_pd(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm_load_pd(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm_load_pd(soa2_quadrupoles_M_z + j);

						sum_M2_x = _mm_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm_add_pd(sum_M2_z, M2_z);

						_mm_store_pd(soa2_quadrupoles_M_x + j, sum_M2_x);
						_mm_store_pd(soa2_quadrupoles_M_y + j, sum_M2_y);
						_mm_store_pd(soa2_quadrupoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_idx + local_i, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_idx + local_i, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_idx + local_i, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_idx + local_i, sum_M1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecQuadrupoles<MacroPolicy>(soa1, i_quadrupole_idx + local_i, soa2, j, soa2_quadrupoles_dist_lookup + j);
				}

			}

			// Computation of charge-quadrupole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{
				const vcp_double_vec q = _mm_loaddup_pd(soa1_charges_q + i_charge_quadrupole_idx);
				const vcp_double_vec r1_x = _mm_loaddup_pd(soa1_charges_r_x + i_charge_quadrupole_idx);
				const vcp_double_vec r1_y = _mm_loaddup_pd(soa1_charges_r_y + i_charge_quadrupole_idx);
				const vcp_double_vec r1_z = _mm_loaddup_pd(soa1_charges_r_z + i_charge_quadrupole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_quadrupole_idx);
				for (; j < end_quadrupoles_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {

						const vcp_double_vec m = _mm_load_pd(soa2_quadrupoles_m + j);
						const vcp_double_vec e_x = _mm_load_pd(soa2_quadrupoles_e_x + j);
						const vcp_double_vec e_y = _mm_load_pd(soa2_quadrupoles_e_y + j);
						const vcp_double_vec e_z = _mm_load_pd(soa2_quadrupoles_e_z + j);
						const vcp_double_vec r2_x = _mm_load_pd(soa2_quadrupoles_r_x + j);
						const vcp_double_vec r2_y = _mm_load_pd(soa2_quadrupoles_r_y + j);
						const vcp_double_vec r2_z = _mm_load_pd(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_quadrupoles_m_r_z + j);


						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeQuadrupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = _mm_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_quadrupoles_f_z + j);

						sum_f2_x = _mm_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_sub_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_quadrupoles_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_quadrupoles_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque
						vcp_double_vec sum_M_x = _mm_load_pd(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M_y = _mm_load_pd(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M_z = _mm_load_pd(soa2_quadrupoles_M_z + j);

						sum_M_x = _mm_add_pd(sum_M_x, M_x);
						sum_M_y = _mm_add_pd(sum_M_y, M_y);
						sum_M_z = _mm_add_pd(sum_M_z, M_z);

						_mm_store_pd(soa2_quadrupoles_M_x + j, sum_M_x);
						_mm_store_pd(soa2_quadrupoles_M_y + j, sum_M_y);
						_mm_store_pd(soa2_quadrupoles_M_z + j, sum_M_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_quadrupole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_quadrupole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_quadrupole_idx, sum_f1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecChargesQuadrupoles<MacroPolicy>(soa1, i_charge_quadrupole_idx, soa2, j, soa2_quadrupoles_dist_lookup + j, false);
				}

				i_charge_quadrupole_idx++;
			}

			// Computation of dipole-quadrupole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++)
			{
				const vcp_double_vec p = _mm_loaddup_pd(soa1_dipoles_p + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_x = _mm_loaddup_pd(soa1_dipoles_e_x + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_y = _mm_loaddup_pd(soa1_dipoles_e_y + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_z = _mm_loaddup_pd(soa1_dipoles_e_z + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_x = _mm_loaddup_pd(soa1_dipoles_r_x + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_y = _mm_loaddup_pd(soa1_dipoles_r_y + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_z = _mm_loaddup_pd(soa1_dipoles_r_z + i_dipole_quadrupole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_quadrupole_idx);
				for (; j < end_quadrupoles_j; j += 2) {
					const vcp_double_vec forceMask = _mm_load_pd(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {

						const vcp_double_vec m = _mm_load_pd(soa2_quadrupoles_m + j);
						const vcp_double_vec ejj_x = _mm_load_pd(soa2_quadrupoles_e_x + j);
						const vcp_double_vec ejj_y = _mm_load_pd(soa2_quadrupoles_e_y + j);
						const vcp_double_vec ejj_z = _mm_load_pd(soa2_quadrupoles_e_z + j);
						const vcp_double_vec rjj_x = _mm_load_pd(soa2_quadrupoles_r_x + j);
						const vcp_double_vec rjj_y = _mm_load_pd(soa2_quadrupoles_r_y + j);
						const vcp_double_vec rjj_z = _mm_load_pd(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm_load_pd(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm_load_pd(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm_load_pd(soa2_quadrupoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipoleQuadrupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, p,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, m,
								f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = _mm_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm_load_pd(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm_load_pd(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm_load_pd(soa2_quadrupoles_f_z + j);

						sum_f2_x = _mm_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm_sub_pd(sum_f2_z, f_z);

						_mm_store_pd(soa2_quadrupoles_f_x + j, sum_f2_x);
						_mm_store_pd(soa2_quadrupoles_f_y + j, sum_f2_y);
						_mm_store_pd(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm_load_pd(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm_load_pd(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm_load_pd(soa2_quadrupoles_M_z + j);

						sum_M2_x = _mm_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm_add_pd(sum_M2_z, M2_z);

						_mm_store_pd(soa2_quadrupoles_M_x + j, sum_M2_x);
						_mm_store_pd(soa2_quadrupoles_M_y + j, sum_M2_y);
						_mm_store_pd(soa2_quadrupoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_quadrupole_idx, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_quadrupole_idx, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_quadrupole_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_quadrupole_idx, sum_M1_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_quadrupole_idx, sum_M1_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_quadrupole_idx, sum_M1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecDipolesQuadrupoles<MacroPolicy>(soa1, i_dipole_quadrupole_idx, soa2, j, soa2_quadrupoles_dist_lookup + j, false);
				}
				i_dipole_quadrupole_idx++;

			}

			i_quadrupole_idx += soa1_mol_quadrupoles_num[i];
		}
	}

	hSum_Add_Store(&_upot6lj, sum_upot6lj);
	hSum_Add_Store(&_upotXpoles, sum_upotXpoles);
	hSum_Add_Store(&_virial, sum_virial);
	hSum_Add_Store(&_myRF, _mm_sub_pd(zero, sum_myRF));


#elif VCP_VEC_TYPE==VCP_VEC_AVX

	static const vcp_double_vec ones = _mm256_castsi256_pd(_mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0));
	static const __m256i memoryMask_first_second = _mm256_set_epi32(0, 0, 0, 0, 1<<31, 0, 1<<31, 0);

	vcp_double_vec sum_upot6lj = vcp_simd_zerov();
	vcp_double_vec sum_upotXpoles = vcp_simd_zerov();
	vcp_double_vec sum_virial = vcp_simd_zerov();
	vcp_double_vec sum_myRF = vcp_simd_zerov();

	const vcp_double_vec cutoffRadiusSquare = _mm256_set1_pd(_cutoffRadiusSquare);
	const vcp_double_vec epsRFInvrc3 = _mm256_broadcast_sd(&_epsRFInvrc3);
	const vcp_double_vec rc2 = _mm256_set1_pd(_LJcutoffRadiusSquare);

	const size_t end_ljc_j = soa2._ljc_num & (~3);
	const size_t end_charges_j = soa2._charges_num & (~3);
	const size_t end_dipoles_j = soa2._dipoles_num & (~3);
	const size_t end_quadrupoles_j = soa2._quadrupoles_num & (~3);

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

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._mol_num; ++i) {
		const vcp_double_vec m1_r_x = _mm256_broadcast_sd(soa1_mol_pos_x + i);
		const vcp_double_vec m1_r_y = _mm256_broadcast_sd(soa1_mol_pos_y + i);
		const vcp_double_vec m1_r_z = _mm256_broadcast_sd(soa1_mol_pos_z + i);

		// Iterate over centers of second cell
		const vcp_double_vec compute_molecule_ljc = calcDistLookup<ForcePolicy>(soa1, i, i_ljc_idx, soa2._ljc_num, _LJcutoffRadiusSquare,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_double_vec compute_molecule_charges = calcDistLookup<ForcePolicy>(soa1, i, i_charge_idx, soa2._charges_num, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_double_vec compute_molecule_dipoles = calcDistLookup<ForcePolicy>(soa1, i, i_dipole_idx, soa2._dipoles_num, _cutoffRadiusSquare,
				soa2_dipoles_dist_lookup, soa2_dipoles_m_r_x, soa2_dipoles_m_r_y, soa2_dipoles_m_r_z,
				cutoffRadiusSquare,	end_dipoles_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_double_vec compute_molecule_quadrupoles = calcDistLookup<ForcePolicy>(soa1, i, i_quadrupole_idx, soa2._quadrupoles_num, _cutoffRadiusSquare,
				soa2_quadrupoles_dist_lookup, soa2_quadrupoles_m_r_x, soa2_quadrupoles_m_r_y, soa2_quadrupoles_m_r_z,
				cutoffRadiusSquare, end_quadrupoles_j, m1_r_x, m1_r_y, m1_r_z);

		if (!_mm256_movemask_pd(compute_molecule_ljc)) {
			i_ljc_idx += soa1_mol_ljc_num[i];
		}
		else {
			// LJ force computation
			for (int local_i = 0; local_i < soa1_mol_ljc_num[i]; local_i++) {
				vcp_double_vec sum_fx1 = vcp_simd_zerov();
				vcp_double_vec sum_fy1 = vcp_simd_zerov();
				vcp_double_vec sum_fz1 = vcp_simd_zerov();
				const vcp_double_vec c_r_x1 = _mm256_broadcast_sd(soa1_ljc_r_x + i_ljc_idx);
				const vcp_double_vec c_r_y1 = _mm256_broadcast_sd(soa1_ljc_r_y + i_ljc_idx);
				const vcp_double_vec c_r_z1 = _mm256_broadcast_sd(soa1_ljc_r_z + i_ljc_idx);
				// Iterate over each pair of centers in the second cell.
				size_t j = ForcePolicy::InitJ(i_ljc_idx);
				for (; j < end_ljc_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_ljc_dist_lookup + j);
					// Only go on if at least 1 of the forces has to be calculated.
					if (_mm256_movemask_pd(forceMask) > 0) {
						const vcp_double_vec c_r_x2 = _mm256_load_pd(soa2_ljc_r_x + j);
						const vcp_double_vec c_r_y2 = _mm256_load_pd(soa2_ljc_r_y + j);
						const vcp_double_vec c_r_z2 = _mm256_load_pd(soa2_ljc_r_z + j);


						const size_t id_i = soa1_ljc_id[i_ljc_idx];

						const size_t id_j0 = soa2_ljc_id[j];
						const size_t id_j1 = soa2_ljc_id[j + 1];
						const size_t id_j2 = soa2_ljc_id[j + 2];
						const size_t id_j3 = soa2_ljc_id[j + 3];
						const vcp_double_vec e0s0 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j0, memoryMask_first_second);
						const vcp_double_vec e1s1 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j1, memoryMask_first_second);
						const vcp_double_vec e2s2 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j2, memoryMask_first_second);
						const vcp_double_vec e3s3 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j3, memoryMask_first_second);

						const vcp_double_vec m_r_x2 = _mm256_load_pd(soa2_ljc_m_r_x + j);
						const vcp_double_vec m_r_y2 = _mm256_load_pd(soa2_ljc_m_r_y + j);
						const vcp_double_vec m_r_z2 = _mm256_load_pd(soa2_ljc_m_r_z + j);
						vcp_double_vec fx, fy, fz;
						//begin loop_body
						_loopBodyLJ<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z,
								c_r_x1,  c_r_y1, c_r_z1,
								m_r_x2, m_r_y2, m_r_z2,
								c_r_x2, c_r_y2, c_r_z2,
								fx, fy, fz,
								sum_upot6lj, sum_virial,
								forceMask, e0s0,
								e1s1, e2s2, e3s3,
								id_j0, id_j1, id_j2, id_j3,
								id_i);

						//end loop_body_lj

						const vcp_double_vec old_fx2 = _mm256_load_pd(soa2_ljc_f_x + j);
						const vcp_double_vec new_fx2 = _mm256_sub_pd(old_fx2, fx);
						_mm256_store_pd(soa2_ljc_f_x + j, new_fx2);
						const vcp_double_vec old_fy2 = _mm256_load_pd(soa2_ljc_f_y + j);
						const vcp_double_vec new_fy2 = _mm256_sub_pd(old_fy2, fy);
						_mm256_store_pd(soa2_ljc_f_y + j, new_fy2);
						const vcp_double_vec old_fz2 = _mm256_load_pd(soa2_ljc_f_z + j);
						const vcp_double_vec new_fz2 = _mm256_sub_pd(old_fz2, fz);
						_mm256_store_pd(soa2_ljc_f_z + j, new_fz2);
						sum_fx1 = _mm256_add_pd(sum_fx1, fx);
						sum_fy1 = _mm256_add_pd(sum_fy1, fy);
						sum_fz1 = _mm256_add_pd(sum_fz1, fz);
					}
				}

				hSum_Add_Store(soa1_ljc_f_x + i_ljc_idx, sum_fx1);
				hSum_Add_Store(soa1_ljc_f_y + i_ljc_idx, sum_fy1);
				hSum_Add_Store(soa1_ljc_f_z + i_ljc_idx, sum_fz1);

				// Unvectorized calculation for leftover pairs.
				for (; j < soa2._ljc_num; ++j) {
					_loopBodyNovecLJ<MacroPolicy>(soa1, i_ljc_idx, soa2, j, soa2_ljc_dist_lookup + j);
				}

				i_ljc_idx++;
			}
		}

		// Computation of site interactions with charges

		if (!_mm256_movemask_pd(compute_molecule_charges)) {
			i_charge_idx += soa1_mol_charges_num[i];
			i_dipole_charge_idx += soa1_mol_dipoles_num[i];
			i_quadrupole_charge_idx += soa1_mol_quadrupoles_num[i];
		}
		else {
			// Computation of charge-charge interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++) {

				const vcp_double_vec q1 = _mm256_broadcast_sd(soa1_charges_q + i_charge_idx + local_i);
				const vcp_double_vec r1_x = _mm256_broadcast_sd(soa1_charges_r_x + i_charge_idx + local_i);
				const vcp_double_vec r1_y = _mm256_broadcast_sd(soa1_charges_r_y + i_charge_idx + local_i);
				const vcp_double_vec r1_z = _mm256_broadcast_sd(soa1_charges_r_z + i_charge_idx + local_i);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_charge_idx + local_i);
				for (; j < end_charges_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {
						const vcp_double_vec q2 = _mm256_load_pd(soa2_charges_q + j);

						const vcp_double_vec r2_x = _mm256_load_pd(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = _mm256_load_pd(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = _mm256_load_pd(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z;
						_loopBodyCharge<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								sum_upotXpoles, sum_virial,
								forceMask);



						sum_f1_x = _mm256_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_add_pd(sum_f1_z, f_z);

						vcp_double_vec f2_x = _mm256_load_pd(soa2_charges_f_x + j);
						vcp_double_vec f2_y = _mm256_load_pd(soa2_charges_f_y + j);
						vcp_double_vec f2_z = _mm256_load_pd(soa2_charges_f_z + j);

						f2_x = _mm256_sub_pd(f2_x, f_x);
						f2_y = _mm256_sub_pd(f2_y, f_y);
						f2_z = _mm256_sub_pd(f2_z, f_z);

						_mm256_store_pd(soa2_charges_f_x + j, f2_x);
						_mm256_store_pd(soa2_charges_f_y + j, f2_y);
						_mm256_store_pd(soa2_charges_f_z + j, f2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx + local_i, sum_f1_z);


				// End iteration over centers with possible left over center
				for (; j < soa2._charges_num; ++j) {
					_loopBodyNovecCharges<MacroPolicy>(soa1, i_charge_idx + local_i, soa2, j, soa2_charges_dist_lookup + j);
				}

			}

			// Computation of dipole-charge interactions

			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++)
			{
				const vcp_double_vec p = _mm256_broadcast_sd(soa1_dipoles_p + i_dipole_charge_idx);
				const vcp_double_vec e_x = _mm256_broadcast_sd(soa1_dipoles_e_x + i_dipole_charge_idx);
				const vcp_double_vec e_y = _mm256_broadcast_sd(soa1_dipoles_e_y + i_dipole_charge_idx);
				const vcp_double_vec e_z = _mm256_broadcast_sd(soa1_dipoles_e_z + i_dipole_charge_idx);
				const vcp_double_vec r1_x = _mm256_broadcast_sd(soa1_dipoles_r_x + i_dipole_charge_idx);
				const vcp_double_vec r1_y = _mm256_broadcast_sd(soa1_dipoles_r_y + i_dipole_charge_idx);
				const vcp_double_vec r1_z = _mm256_broadcast_sd(soa1_dipoles_r_z + i_dipole_charge_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M_x = vcp_simd_zerov();
				vcp_double_vec sum_M_y = vcp_simd_zerov();
				vcp_double_vec sum_M_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {

						const vcp_double_vec q = _mm256_load_pd(soa2_charges_q + j);

						const vcp_double_vec r2_x = _mm256_load_pd(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = _mm256_load_pd(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = _mm256_load_pd(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeDipole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = _mm256_sub_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_sub_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_sub_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_charges_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_charges_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_charges_f_z + j);

						sum_f2_x = _mm256_add_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_add_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_add_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_charges_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_charges_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_charges_f_z + j, sum_f2_z);

						// Store torque

						sum_M_x = _mm256_add_pd(sum_M_x, M_x);
						sum_M_y = _mm256_add_pd(sum_M_y, M_y);
						sum_M_z = _mm256_add_pd(sum_M_z, M_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_charge_idx, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_charge_idx, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_charge_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_charge_idx, sum_M_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_charge_idx, sum_M_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_charge_idx, sum_M_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._charges_num; ++j) {
					_loopBodyNovecChargesDipoles<MacroPolicy>(soa2, j, soa1, i_dipole_charge_idx, soa2_charges_dist_lookup + j, true);
				}

				i_dipole_charge_idx++;
			}

			// Computation of quadrupole-charge interactions

			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++)
			{
				const vcp_double_vec m = _mm256_broadcast_sd(soa1_quadrupoles_m + i_quadrupole_charge_idx);
				const vcp_double_vec e_x = _mm256_broadcast_sd(soa1_quadrupoles_e_x + i_quadrupole_charge_idx);
				const vcp_double_vec e_y = _mm256_broadcast_sd(soa1_quadrupoles_e_y + i_quadrupole_charge_idx);
				const vcp_double_vec e_z = _mm256_broadcast_sd(soa1_quadrupoles_e_z + i_quadrupole_charge_idx);
				const vcp_double_vec r1_x = _mm256_broadcast_sd(soa1_quadrupoles_r_x + i_quadrupole_charge_idx);
				const vcp_double_vec r1_y = _mm256_broadcast_sd(soa1_quadrupoles_r_y + i_quadrupole_charge_idx);
				const vcp_double_vec r1_z = _mm256_broadcast_sd(soa1_quadrupoles_r_z + i_quadrupole_charge_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {

						const vcp_double_vec q = _mm256_load_pd(soa2_charges_q + j);

						const vcp_double_vec r2_x = _mm256_load_pd(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = _mm256_load_pd(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = _mm256_load_pd(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeQuadrupole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = _mm256_sub_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_sub_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_sub_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_charges_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_charges_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_charges_f_z + j);

						sum_f2_x = _mm256_add_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_add_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_add_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_charges_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_charges_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_charges_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm256_add_pd(sum_M1_x, M_x);
						sum_M1_y = _mm256_add_pd(sum_M1_y, M_y);
						sum_M1_z = _mm256_add_pd(sum_M1_z, M_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_charge_idx, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_charge_idx, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_charge_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_charge_idx, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_charge_idx, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_charge_idx, sum_M1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._charges_num; ++j) {
					_loopBodyNovecChargesQuadrupoles<MacroPolicy>(soa2, j, soa1, i_quadrupole_charge_idx, soa2_charges_dist_lookup + j, true);
				}
				i_quadrupole_charge_idx++;
			}

			i_charge_idx += soa1_mol_charges_num[i];
		}

		// Computation of site interactions with dipoles

		// Continue with next molecule if no force has to be calculated
		if (!_mm256_movemask_pd(compute_molecule_dipoles)) {
			i_dipole_idx += soa1_mol_dipoles_num[i];
			i_charge_dipole_idx += soa1_mol_charges_num[i];
			i_quadrupole_dipole_idx += soa1_mol_quadrupoles_num[i];
		}
		else {
			// Computation of dipole-dipole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++) {

				const vcp_double_vec p1 = _mm256_broadcast_sd(soa1_dipoles_p + i_dipole_idx + local_i);
				const vcp_double_vec e1_x = _mm256_broadcast_sd(soa1_dipoles_e_x + i_dipole_idx + local_i);
				const vcp_double_vec e1_y = _mm256_broadcast_sd(soa1_dipoles_e_y + i_dipole_idx + local_i);
				const vcp_double_vec e1_z = _mm256_broadcast_sd(soa1_dipoles_e_z + i_dipole_idx + local_i);
				const vcp_double_vec r1_x = _mm256_broadcast_sd(soa1_dipoles_r_x + i_dipole_idx + local_i);
				const vcp_double_vec r1_y = _mm256_broadcast_sd(soa1_dipoles_r_y + i_dipole_idx + local_i);
				const vcp_double_vec r1_z = _mm256_broadcast_sd(soa1_dipoles_r_z + i_dipole_idx + local_i);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_dipole_idx + local_i);
				for (; j < end_dipoles_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {
						const vcp_double_vec p2 = _mm256_load_pd(soa2_dipoles_p + j);
						const vcp_double_vec e2_x = _mm256_load_pd(soa2_dipoles_e_x + j);
						const vcp_double_vec e2_y = _mm256_load_pd(soa2_dipoles_e_y + j);
						const vcp_double_vec e2_z = _mm256_load_pd(soa2_dipoles_e_z + j);
						const vcp_double_vec r2_x = _mm256_load_pd(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = _mm256_load_pd(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = _mm256_load_pd(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipole<MacroPolicy>(
							m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, p1,
							m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p2,
							f_x, f_y, f_z,
							M1_x, M1_y, M1_z,
							M2_x, M2_y, M2_z,
							sum_upotXpoles, sum_virial, sum_myRF,
							forceMask,
							epsRFInvrc3);

						// Store forces

						sum_f1_x = _mm256_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_dipoles_f_z + j);

						sum_f2_x = _mm256_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_sub_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_dipoles_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_dipoles_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm256_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm256_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm256_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm256_load_pd(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm256_load_pd(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm256_load_pd(soa2_dipoles_M_z + j);

						sum_M2_x = _mm256_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm256_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm256_add_pd(sum_M2_z, M2_z);

						_mm256_store_pd(soa2_dipoles_M_x + j, sum_M2_x);
						_mm256_store_pd(soa2_dipoles_M_y + j, sum_M2_y);
						_mm256_store_pd(soa2_dipoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_idx + local_i, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_idx + local_i, sum_M1_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_idx + local_i, sum_M1_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_idx + local_i, sum_M1_z);


				// End iteration over centers with possible left over center
				for (; j < soa2._dipoles_num; ++j) {
					_loopBodyNovecDipoles<MacroPolicy>(soa1, i_dipole_idx + local_i, soa2, j, soa2_dipoles_dist_lookup + j);
				}

			}

			// Computation of charge-dipole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{

				const vcp_double_vec q = _mm256_broadcast_sd(soa1_charges_q + i_charge_dipole_idx);
				const vcp_double_vec r1_x = _mm256_broadcast_sd(soa1_charges_r_x + i_charge_dipole_idx);
				const vcp_double_vec r1_y = _mm256_broadcast_sd(soa1_charges_r_y + i_charge_dipole_idx);
				const vcp_double_vec r1_z = _mm256_broadcast_sd(soa1_charges_r_z + i_charge_dipole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_dipole_idx);
				for (; j < end_dipoles_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {

						const vcp_double_vec p = _mm256_load_pd(soa2_dipoles_p + j);

						const vcp_double_vec e_x = _mm256_load_pd(soa2_dipoles_e_x + j);
						const vcp_double_vec e_y = _mm256_load_pd(soa2_dipoles_e_y + j);
						const vcp_double_vec e_z = _mm256_load_pd(soa2_dipoles_e_z + j);

						const vcp_double_vec r2_x = _mm256_load_pd(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = _mm256_load_pd(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = _mm256_load_pd(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeDipole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = _mm256_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_dipoles_f_z + j);

						sum_f2_x = _mm256_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_sub_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_dipoles_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_dipoles_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						vcp_double_vec sum_M2_x = _mm256_load_pd(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm256_load_pd(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm256_load_pd(soa2_dipoles_M_z + j);

						sum_M2_x = _mm256_add_pd(sum_M2_x, M_x);
						sum_M2_y = _mm256_add_pd(sum_M2_y, M_y);
						sum_M2_z = _mm256_add_pd(sum_M2_z, M_z);

						_mm256_store_pd(soa2_dipoles_M_x + j, sum_M2_x);
						_mm256_store_pd(soa2_dipoles_M_y + j, sum_M2_y);
						_mm256_store_pd(soa2_dipoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_dipole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_dipole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_dipole_idx, sum_f1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._dipoles_num; ++j) {
					_loopBodyNovecChargesDipoles<MacroPolicy>(soa1, i_charge_dipole_idx, soa2, j, soa2_dipoles_dist_lookup + j, false);
				}

				i_charge_dipole_idx++;
			}

			// Computation of quadrupole-dipole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++) {

				const vcp_double_vec m = _mm256_broadcast_sd(soa1_quadrupoles_m + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_x = _mm256_broadcast_sd(soa1_quadrupoles_e_x + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_y = _mm256_broadcast_sd(soa1_quadrupoles_e_y + i_quadrupole_dipole_idx);
				const vcp_double_vec e1_z = _mm256_broadcast_sd(soa1_quadrupoles_e_z + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_x = _mm256_broadcast_sd(soa1_quadrupoles_r_x + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_y = _mm256_broadcast_sd(soa1_quadrupoles_r_y + i_quadrupole_dipole_idx);
				const vcp_double_vec r1_z = _mm256_broadcast_sd(soa1_quadrupoles_r_z + i_quadrupole_dipole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_dipole_idx);
				for (; j < end_dipoles_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {
						const vcp_double_vec p = _mm256_load_pd(soa2_dipoles_p + j);
						const vcp_double_vec e2_x = _mm256_load_pd(soa2_dipoles_e_x + j);
						const vcp_double_vec e2_y = _mm256_load_pd(soa2_dipoles_e_y + j);
						const vcp_double_vec e2_z = _mm256_load_pd(soa2_dipoles_e_z + j);
						const vcp_double_vec r2_x = _mm256_load_pd(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = _mm256_load_pd(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = _mm256_load_pd(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipoleQuadrupole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, m,
								f_x, f_y, f_z, M2_x, M2_y, M2_z, M1_x, M1_y, M1_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = _mm256_sub_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_sub_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_sub_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_dipoles_f_z + j);

						sum_f2_x = _mm256_add_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_add_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_add_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_dipoles_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_dipoles_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm256_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm256_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm256_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm256_load_pd(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm256_load_pd(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm256_load_pd(soa2_dipoles_M_z + j);

						sum_M2_x = _mm256_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm256_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm256_add_pd(sum_M2_z, M2_z);

						_mm256_store_pd(soa2_dipoles_M_x + j, sum_M2_x);
						_mm256_store_pd(soa2_dipoles_M_y + j, sum_M2_y);
						_mm256_store_pd(soa2_dipoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_dipole_idx, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_dipole_idx, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_dipole_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_dipole_idx, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_dipole_idx, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_dipole_idx, sum_M1_z);


				// End iteration over centers with possible left over center
				for (; j < soa2._dipoles_num; ++j) {
					_loopBodyNovecDipolesQuadrupoles<MacroPolicy>(soa2, j, soa1, i_quadrupole_dipole_idx, soa2_dipoles_dist_lookup + j, true);
				}

				i_quadrupole_dipole_idx++;

			}

			i_dipole_idx += soa1_mol_dipoles_num[i];
		}

		// Computation of site interactions with quadrupoles

		if (!_mm256_movemask_pd(compute_molecule_quadrupoles)) {
			i_quadrupole_idx += soa1_mol_quadrupoles_num[i];
			i_charge_quadrupole_idx += soa1_mol_charges_num[i];
			i_dipole_quadrupole_idx += soa1_mol_dipoles_num[i];
		}
		else {
			// Computation of quadrupole-quadrupole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_quadrupoles_num[i]; local_i++)
			{
				const vcp_double_vec mii = _mm256_broadcast_sd(soa1_quadrupoles_m + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_x = _mm256_broadcast_sd(soa1_quadrupoles_e_x + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_y = _mm256_broadcast_sd(soa1_quadrupoles_e_y + i_quadrupole_idx + local_i);
				const vcp_double_vec eii_z = _mm256_broadcast_sd(soa1_quadrupoles_e_z + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_x = _mm256_broadcast_sd(soa1_quadrupoles_r_x + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_y = _mm256_broadcast_sd(soa1_quadrupoles_r_y + i_quadrupole_idx + local_i);
				const vcp_double_vec rii_z = _mm256_broadcast_sd(soa1_quadrupoles_r_z + i_quadrupole_idx + local_i);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_quadrupole_idx + local_i);
				for (; j < end_quadrupoles_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {

						const vcp_double_vec mjj = _mm256_load_pd(soa2_quadrupoles_m + j);
						const vcp_double_vec ejj_x = _mm256_load_pd(soa2_quadrupoles_e_x + j);
						const vcp_double_vec ejj_y = _mm256_load_pd(soa2_quadrupoles_e_y + j);
						const vcp_double_vec ejj_z = _mm256_load_pd(soa2_quadrupoles_e_z + j);
						const vcp_double_vec rjj_x = _mm256_load_pd(soa2_quadrupoles_r_x + j);
						const vcp_double_vec rjj_y = _mm256_load_pd(soa2_quadrupoles_r_y + j);
						const vcp_double_vec rjj_z = _mm256_load_pd(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_quadrupoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyQudarupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, mii,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, mjj,
								f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								forceMask);

						// Store forces

						sum_f1_x = _mm256_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_quadrupoles_f_z + j);

						sum_f2_x = _mm256_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_sub_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_quadrupoles_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_quadrupoles_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm256_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm256_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm256_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm256_load_pd(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm256_load_pd(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm256_load_pd(soa2_quadrupoles_M_z + j);

						sum_M2_x = _mm256_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm256_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm256_add_pd(sum_M2_z, M2_z);

						_mm256_store_pd(soa2_quadrupoles_M_x + j, sum_M2_x);
						_mm256_store_pd(soa2_quadrupoles_M_y + j, sum_M2_y);
						_mm256_store_pd(soa2_quadrupoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_quadrupoles_f_x + i_quadrupole_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_quadrupoles_f_y + i_quadrupole_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_quadrupoles_f_z + i_quadrupole_idx + local_i, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_quadrupoles_M_x + i_quadrupole_idx + local_i, sum_M1_x);
				hSum_Add_Store(soa1_quadrupoles_M_y + i_quadrupole_idx + local_i, sum_M1_y);
				hSum_Add_Store(soa1_quadrupoles_M_z + i_quadrupole_idx + local_i, sum_M1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecQuadrupoles<MacroPolicy>(soa1, i_quadrupole_idx + local_i, soa2, j, soa2_quadrupoles_dist_lookup + j);
				}

			}

			// Computation of charge-quadrupole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{
				const vcp_double_vec q = _mm256_broadcast_sd(soa1_charges_q + i_charge_quadrupole_idx);
				const vcp_double_vec r1_x = _mm256_broadcast_sd(soa1_charges_r_x + i_charge_quadrupole_idx);
				const vcp_double_vec r1_y = _mm256_broadcast_sd(soa1_charges_r_y + i_charge_quadrupole_idx);
				const vcp_double_vec r1_z = _mm256_broadcast_sd(soa1_charges_r_z + i_charge_quadrupole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_quadrupole_idx);
				for (; j < end_quadrupoles_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {

						const vcp_double_vec m = _mm256_load_pd(soa2_quadrupoles_m + j);
						const vcp_double_vec e_x = _mm256_load_pd(soa2_quadrupoles_e_x + j);
						const vcp_double_vec e_y = _mm256_load_pd(soa2_quadrupoles_e_y + j);
						const vcp_double_vec e_z = _mm256_load_pd(soa2_quadrupoles_e_z + j);
						const vcp_double_vec r2_x = _mm256_load_pd(soa2_quadrupoles_r_x + j);
						const vcp_double_vec r2_y = _mm256_load_pd(soa2_quadrupoles_r_y + j);
						const vcp_double_vec r2_z = _mm256_load_pd(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_quadrupoles_m_r_z + j);


						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeQuadrupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = _mm256_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_quadrupoles_f_z + j);

						sum_f2_x = _mm256_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_sub_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_quadrupoles_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_quadrupoles_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque
						vcp_double_vec sum_M_x = _mm256_load_pd(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M_y = _mm256_load_pd(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M_z = _mm256_load_pd(soa2_quadrupoles_M_z + j);

						sum_M_x = _mm256_add_pd(sum_M_x, M_x);
						sum_M_y = _mm256_add_pd(sum_M_y, M_y);
						sum_M_z = _mm256_add_pd(sum_M_z, M_z);

						_mm256_store_pd(soa2_quadrupoles_M_x + j, sum_M_x);
						_mm256_store_pd(soa2_quadrupoles_M_y + j, sum_M_y);
						_mm256_store_pd(soa2_quadrupoles_M_z + j, sum_M_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_quadrupole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_quadrupole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_quadrupole_idx, sum_f1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecChargesQuadrupoles<MacroPolicy>(soa1, i_charge_quadrupole_idx, soa2, j, soa2_quadrupoles_dist_lookup + j, false);
				}

				i_charge_quadrupole_idx++;
			}

			// Computation of dipole-quadrupole interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_dipoles_num[i]; local_i++)
			{
				const vcp_double_vec p = _mm256_broadcast_sd(soa1_dipoles_p + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_x = _mm256_broadcast_sd(soa1_dipoles_e_x + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_y = _mm256_broadcast_sd(soa1_dipoles_e_y + i_dipole_quadrupole_idx);
				const vcp_double_vec eii_z = _mm256_broadcast_sd(soa1_dipoles_e_z + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_x = _mm256_broadcast_sd(soa1_dipoles_r_x + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_y = _mm256_broadcast_sd(soa1_dipoles_r_y + i_dipole_quadrupole_idx);
				const vcp_double_vec rii_z = _mm256_broadcast_sd(soa1_dipoles_r_z + i_dipole_quadrupole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_quadrupole_idx);
				for (; j < end_quadrupoles_j; j += 4) {
					const vcp_double_vec forceMask = _mm256_load_pd(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {

						const vcp_double_vec m = _mm256_load_pd(soa2_quadrupoles_m + j);
						const vcp_double_vec ejj_x = _mm256_load_pd(soa2_quadrupoles_e_x + j);
						const vcp_double_vec ejj_y = _mm256_load_pd(soa2_quadrupoles_e_y + j);
						const vcp_double_vec ejj_z = _mm256_load_pd(soa2_quadrupoles_e_z + j);
						const vcp_double_vec rjj_x = _mm256_load_pd(soa2_quadrupoles_r_x + j);
						const vcp_double_vec rjj_y = _mm256_load_pd(soa2_quadrupoles_r_y + j);
						const vcp_double_vec rjj_z = _mm256_load_pd(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = _mm256_load_pd(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = _mm256_load_pd(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = _mm256_load_pd(soa2_quadrupoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipoleQuadrupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, p,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, m,
								f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = _mm256_add_pd(sum_f1_x, f_x);
						sum_f1_y = _mm256_add_pd(sum_f1_y, f_y);
						sum_f1_z = _mm256_add_pd(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = _mm256_load_pd(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = _mm256_load_pd(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = _mm256_load_pd(soa2_quadrupoles_f_z + j);

						sum_f2_x = _mm256_sub_pd(sum_f2_x, f_x);
						sum_f2_y = _mm256_sub_pd(sum_f2_y, f_y);
						sum_f2_z = _mm256_sub_pd(sum_f2_z, f_z);

						_mm256_store_pd(soa2_quadrupoles_f_x + j, sum_f2_x);
						_mm256_store_pd(soa2_quadrupoles_f_y + j, sum_f2_y);
						_mm256_store_pd(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = _mm256_add_pd(sum_M1_x, M1_x);
						sum_M1_y = _mm256_add_pd(sum_M1_y, M1_y);
						sum_M1_z = _mm256_add_pd(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = _mm256_load_pd(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M2_y = _mm256_load_pd(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M2_z = _mm256_load_pd(soa2_quadrupoles_M_z + j);

						sum_M2_x = _mm256_add_pd(sum_M2_x, M2_x);
						sum_M2_y = _mm256_add_pd(sum_M2_y, M2_y);
						sum_M2_z = _mm256_add_pd(sum_M2_z, M2_z);

						_mm256_store_pd(soa2_quadrupoles_M_x + j, sum_M2_x);
						_mm256_store_pd(soa2_quadrupoles_M_y + j, sum_M2_y);
						_mm256_store_pd(soa2_quadrupoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_dipoles_f_x + i_dipole_quadrupole_idx, sum_f1_x);
				hSum_Add_Store(soa1_dipoles_f_y + i_dipole_quadrupole_idx, sum_f1_y);
				hSum_Add_Store(soa1_dipoles_f_z + i_dipole_quadrupole_idx, sum_f1_z);

				// Add old torques and summed calculated torques for center 1
				hSum_Add_Store(soa1_dipoles_M_x + i_dipole_quadrupole_idx, sum_M1_x);
				hSum_Add_Store(soa1_dipoles_M_y + i_dipole_quadrupole_idx, sum_M1_y);
				hSum_Add_Store(soa1_dipoles_M_z + i_dipole_quadrupole_idx, sum_M1_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._quadrupoles_num; ++j) {
					_loopBodyNovecDipolesQuadrupoles<MacroPolicy>(soa1, i_dipole_quadrupole_idx, soa2, j, soa2_quadrupoles_dist_lookup + j, false);
				}

				i_dipole_quadrupole_idx++;

			}

			i_quadrupole_idx += soa1_mol_quadrupoles_num[i];
		}
	}

	hSum_Add_Store(&_upot6lj, sum_upot6lj);
	hSum_Add_Store(&_upotXpoles, sum_upotXpoles);
	hSum_Add_Store(&_virial, sum_virial);
	hSum_Add_Store(&_myRF, _mm256_sub_pd(zero, sum_myRF));

#endif
#endif
} // void LennardJonesCellHandler::CalculatePairs_(LJSoA & soa1, LJSoA & soa2)

void VectorizedCellProcessor::processCell(ParticleCell & c) {
	assert(c.getCellDataSoA());
	if (c.isHaloCell() || (c.getCellDataSoA()->_mol_num < 2)) {
		return;
	}

	_calculatePairs<SingleCellPolicy_, AllMacroPolicy_>(*(c.getCellDataSoA()), *(c.getCellDataSoA()));
}

void VectorizedCellProcessor::processCellPair(ParticleCell & c1, ParticleCell & c2) {
	assert(&c1 != &c2);
	assert(c1.getCellDataSoA());
	assert(c2.getCellDataSoA());

	if ((c1.getCellDataSoA()->_mol_num == 0) || (c2.getCellDataSoA()->_mol_num == 0)) {
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
