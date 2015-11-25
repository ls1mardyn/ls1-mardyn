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
		CellProcessor(cutoffRadius, LJcutoffRadius), _domain(domain),
		// maybe move the following to somewhere else:
		_epsRFInvrc3(2. * (domain.getepsilonRF() - 1.) / ((cutoffRadius * cutoffRadius * cutoffRadius) * (2. * domain.getepsilonRF() + 1.))), 
		_compIDs(), _eps_sig(), _shift6(), _upot6lj(0.0), _upotXpoles(0.0), _virial(0.0), _myRF(0.0) {

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

#if VCP_VEC_TYPE==VCP_VEC_SSE3 or VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_NOVEC
	//const vcp_double_vec minus_one = vcp_simd_set1(-1.0); //currently not used, would produce warning
	const vcp_double_vec zero = vcp_simd_zerov();
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
		const vcp_double_vec c_dx = vcp_simd_sub(r1_x, r2_x);
		const vcp_double_vec c_dy = vcp_simd_sub(r1_y, r2_y);
		const vcp_double_vec c_dz = vcp_simd_sub(r1_z, r2_z);

		const vcp_double_vec c_dx2 = vcp_simd_mul(c_dx, c_dx);
		const vcp_double_vec c_dy2 = vcp_simd_mul(c_dy, c_dy);
		const vcp_double_vec c_dz2 = vcp_simd_mul(c_dz, c_dz);

		const vcp_double_vec c_dr2 = vcp_simd_add(vcp_simd_add(c_dx2, c_dy2), c_dz2);
		const vcp_double_vec c_dr2_inv_unmasked = vcp_simd_div(one, c_dr2);
		const vcp_double_vec c_dr2_inv = vcp_simd_applymask(c_dr2_inv_unmasked, forceMask);//masked
	    const vcp_double_vec c_dr_inv = vcp_simd_sqrt(c_dr2_inv);//masked

		const vcp_double_vec q1q2per4pie0 = vcp_simd_mul(qii, qjj);
		const vcp_double_vec upot = vcp_simd_mul(q1q2per4pie0, c_dr_inv);//masked
		const vcp_double_vec fac = vcp_simd_mul(upot, c_dr2_inv);//masked

		f_x = vcp_simd_mul(c_dx, fac);
		f_y = vcp_simd_mul(c_dy, fac);
		f_z = vcp_simd_mul(c_dz, fac);

		const vcp_double_vec m_dx = vcp_simd_sub(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = vcp_simd_sub(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = vcp_simd_sub(m1_r_z, m2_r_z);

		const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (vcp_simd_movemask(macroMask) > 0) {
			const vcp_double_vec upot_masked = vcp_simd_applymask(upot, macroMask);
			sum_upotXpoles = vcp_simd_add(sum_upotXpoles, upot_masked);

			const vcp_double_vec virial_x = vcp_simd_mul(m_dx, f_x);
			const vcp_double_vec virial_y = vcp_simd_mul(m_dy, f_y);
			const vcp_double_vec virial_z = vcp_simd_mul(m_dz, f_z);

			const vcp_double_vec virial = vcp_simd_add(vcp_simd_add(virial_x, virial_y), virial_z);
			const vcp_double_vec virial_masked = vcp_simd_applymask(virial, macroMask);
			sum_virial = vcp_simd_add(sum_virial, virial_masked);

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
		const vcp_double_vec dx = vcp_simd_sub(r1_x, r2_x);
		const vcp_double_vec dy = vcp_simd_sub(r1_y, r2_y);
		const vcp_double_vec dz = vcp_simd_sub(r1_z, r2_z);

		const vcp_double_vec dx2 = vcp_simd_mul(dx, dx);
		const vcp_double_vec dy2 = vcp_simd_mul(dy, dy);
		const vcp_double_vec dz2 = vcp_simd_mul(dz, dz);

		const vcp_double_vec dr2 = vcp_simd_add(vcp_simd_add(dx2, dy2), dz2);

		const vcp_double_vec dr2_inv_unmasked = vcp_simd_div(one, dr2);
		const vcp_double_vec dr2_inv = vcp_simd_applymask(dr2_inv_unmasked, forceMask);
		const vcp_double_vec dr_inv = vcp_simd_sqrt(dr2_inv);
		const vcp_double_vec dr3_inv = vcp_simd_mul(dr2_inv, dr_inv);

		const vcp_double_vec re = vcp_simd_add(vcp_simd_mul(dx, e_x), vcp_simd_add(vcp_simd_mul(dy, e_y), vcp_simd_mul(dz, e_z)));

		const vcp_double_vec qpper4pie0 = vcp_simd_mul(q, p);
		const vcp_double_vec qpper4pie0dr3 = vcp_simd_mul(qpper4pie0, dr3_inv);

		const vcp_double_vec fac = vcp_simd_mul(dr2_inv, vcp_simd_mul(three, re));

		f_x = vcp_simd_mul(qpper4pie0dr3, vcp_simd_sub(e_x, vcp_simd_mul(dx, fac)));
		f_y = vcp_simd_mul(qpper4pie0dr3, vcp_simd_sub(e_y, vcp_simd_mul(dy, fac)));
		f_z = vcp_simd_mul(qpper4pie0dr3, vcp_simd_sub(e_z, vcp_simd_mul(dz, fac)));

		const vcp_double_vec m_dx = vcp_simd_sub(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = vcp_simd_sub(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = vcp_simd_sub(m1_r_z, m2_r_z);

		const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);
		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (vcp_simd_movemask(macroMask) > 0)
		{
			const vcp_double_vec minusUpot_unmasked =  vcp_simd_mul(qpper4pie0dr3, re);
			const vcp_double_vec minusUpot = vcp_simd_applymask(minusUpot_unmasked, macroMask);
			sum_upotXpoles = vcp_simd_sub(sum_upotXpoles, minusUpot);

			const vcp_double_vec virial_x = vcp_simd_mul(m_dx, f_x);
			const vcp_double_vec virial_y = vcp_simd_mul(m_dy, f_y);
			const vcp_double_vec virial_z = vcp_simd_mul(m_dz, f_z);

			const vcp_double_vec virial_unmasked = vcp_simd_add(vcp_simd_add(virial_x, virial_y), virial_z);
			const vcp_double_vec virial = vcp_simd_applymask(virial_unmasked, macroMask);
			sum_virial = vcp_simd_add(sum_virial, virial);
		}

		const vcp_double_vec e_x_dy = vcp_simd_mul(e_x, dy);
		const vcp_double_vec e_x_dz = vcp_simd_mul(e_x, dz);
		const vcp_double_vec e_y_dx = vcp_simd_mul(e_y, dx);
		const vcp_double_vec e_y_dz = vcp_simd_mul(e_y, dz);
		const vcp_double_vec e_z_dx = vcp_simd_mul(e_z, dx);
		const vcp_double_vec e_z_dy = vcp_simd_mul(e_z, dy);

		const vcp_double_vec e_x_dy_minus_e_y_dx = vcp_simd_sub(e_x_dy, e_y_dx);
		const vcp_double_vec e_y_dz_minus_e_z_dy = vcp_simd_sub(e_y_dz, e_z_dy);
		const vcp_double_vec e_z_dx_minus_e_x_dz = vcp_simd_sub(e_z_dx, e_x_dz);

		M_x = vcp_simd_mul(qpper4pie0dr3, e_y_dz_minus_e_z_dy);
		M_y = vcp_simd_mul(qpper4pie0dr3, e_z_dx_minus_e_x_dz);
		M_z = vcp_simd_mul(qpper4pie0dr3, e_x_dy_minus_e_y_dx);
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
		const vcp_double_vec dx = vcp_simd_sub(r1_x, r2_x);
		const vcp_double_vec dy = vcp_simd_sub(r1_y, r2_y);
		const vcp_double_vec dz = vcp_simd_sub(r1_z, r2_z);

		const vcp_double_vec dx2 = vcp_simd_mul(dx, dx);
		const vcp_double_vec dy2 = vcp_simd_mul(dy, dy);
		const vcp_double_vec dz2 = vcp_simd_mul(dz, dz);

		const vcp_double_vec dr2 = vcp_simd_add(vcp_simd_add(dx2, dy2), dz2);

		const vcp_double_vec dr2_inv_unmasked = vcp_simd_div(one, dr2);
		const vcp_double_vec dr2_inv = vcp_simd_applymask(dr2_inv_unmasked, forceMask);
		const vcp_double_vec dr_inv = vcp_simd_sqrt(dr2_inv);
		const vcp_double_vec dr2three_inv = vcp_simd_mul(three, dr2_inv);

		const vcp_double_vec p1p2 = vcp_simd_applymask(vcp_simd_mul(pii, pjj), forceMask);
		const vcp_double_vec p1p2per4pie0 = p1p2;
		const vcp_double_vec rffac = vcp_simd_mul(p1p2, epsRFInvrc3);

		const vcp_double_vec p1p2per4pie0r3 = vcp_simd_mul(p1p2per4pie0, vcp_simd_mul(dr_inv, dr2_inv));
		const vcp_double_vec p1p2threeper4pie0r5 = vcp_simd_mul(p1p2per4pie0r3, dr2three_inv);

		const vcp_double_vec e1e2 = vcp_simd_add(vcp_simd_mul(eii_x, ejj_x), vcp_simd_add(vcp_simd_mul(eii_y, ejj_y), vcp_simd_mul(eii_z, ejj_z)));
		const vcp_double_vec re1 = vcp_simd_add(vcp_simd_mul(dx, eii_x), vcp_simd_add(vcp_simd_mul(dy, eii_y), vcp_simd_mul(dz, eii_z)));
		const vcp_double_vec re2 = vcp_simd_add(vcp_simd_mul(dx, ejj_x), vcp_simd_add(vcp_simd_mul(dy, ejj_y), vcp_simd_mul(dz, ejj_z)));

		const vcp_double_vec re1threeperr2 = vcp_simd_mul(re1, dr2three_inv);
		const vcp_double_vec re2threeperr2 = vcp_simd_mul(re2, dr2three_inv);
		const vcp_double_vec re1re2perr2 = vcp_simd_mul(dr2_inv, vcp_simd_mul(re1, re2));

		const vcp_double_vec e1e2minus5re1re2perr2 = vcp_simd_sub(e1e2, vcp_simd_mul(five, re1re2perr2));


		f_x = vcp_simd_mul(p1p2threeper4pie0r5, vcp_simd_add(vcp_simd_mul(dx, e1e2minus5re1re2perr2), vcp_simd_add(vcp_simd_mul(eii_x, re2), vcp_simd_mul(ejj_x, re1))));
		f_y = vcp_simd_mul(p1p2threeper4pie0r5, vcp_simd_add(vcp_simd_mul(dy, e1e2minus5re1re2perr2), vcp_simd_add(vcp_simd_mul(eii_y, re2), vcp_simd_mul(ejj_y, re1))));
		f_z = vcp_simd_mul(p1p2threeper4pie0r5, vcp_simd_add(vcp_simd_mul(dz, e1e2minus5re1re2perr2), vcp_simd_add(vcp_simd_mul(eii_z, re2), vcp_simd_mul(ejj_z, re1))));

		const vcp_double_vec m_dx = vcp_simd_sub(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = vcp_simd_sub(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = vcp_simd_sub(m1_r_z, m2_r_z);

		const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (vcp_simd_movemask(macroMask) > 0) {
			// can we precompute some of this?
			const vcp_double_vec upot = vcp_simd_mul(p1p2per4pie0r3, vcp_simd_sub(e1e2, vcp_simd_mul(three, re1re2perr2)));
			const vcp_double_vec upot_masked = vcp_simd_applymask(upot, macroMask);
			sum_upotXpoles = vcp_simd_add(sum_upotXpoles, upot_masked);

			const vcp_double_vec virial_x = vcp_simd_mul(m_dx, f_x);
			const vcp_double_vec virial_y = vcp_simd_mul(m_dy, f_y);
			const vcp_double_vec virial_z = vcp_simd_mul(m_dz, f_z);

			const vcp_double_vec virial = vcp_simd_add(vcp_simd_add(virial_x, virial_y), virial_z);
			const vcp_double_vec virial_masked = vcp_simd_applymask(virial, macroMask);
			sum_virial = vcp_simd_add(sum_virial, virial_masked);

			const vcp_double_vec myRF_masked =  vcp_simd_applymask(vcp_simd_mul(rffac, e1e2), macroMask);
			sum_myRF = vcp_simd_add(sum_myRF, myRF_masked);
		}

		const vcp_double_vec e1_x_e2_y = vcp_simd_mul(eii_x, ejj_y);
		const vcp_double_vec e1_x_e2_z = vcp_simd_mul(eii_x, ejj_z);
		const vcp_double_vec e1_y_e2_x = vcp_simd_mul(eii_y, ejj_x);
		const vcp_double_vec e1_y_e2_z = vcp_simd_mul(eii_y, ejj_z);
		const vcp_double_vec e1_z_e2_x = vcp_simd_mul(eii_z, ejj_x);
		const vcp_double_vec e1_z_e2_y = vcp_simd_mul(eii_z, ejj_y);

		const vcp_double_vec e1_x_e2_y_minus_e1_y_e2_x = vcp_simd_sub(e1_x_e2_y, e1_y_e2_x);
		const vcp_double_vec e1_y_e2_z_minus_e1_z_e2_y = vcp_simd_sub(e1_y_e2_z, e1_z_e2_y);
		const vcp_double_vec e1_z_e2_x_minus_e1_x_e2_z = vcp_simd_sub(e1_z_e2_x, e1_x_e2_z);

		M1_x = vcp_simd_add(vcp_simd_mul(p1p2per4pie0r3, vcp_simd_sub(vcp_simd_mul(re2threeperr2, vcp_simd_sub(vcp_simd_mul(eii_y, dz), vcp_simd_mul(eii_z, dy))), e1_y_e2_z_minus_e1_z_e2_y)), vcp_simd_mul(rffac, e1_y_e2_z_minus_e1_z_e2_y));
		M1_y = vcp_simd_add(vcp_simd_mul(p1p2per4pie0r3, vcp_simd_sub(vcp_simd_mul(re2threeperr2, vcp_simd_sub(vcp_simd_mul(eii_z, dx), vcp_simd_mul(eii_x, dz))), e1_z_e2_x_minus_e1_x_e2_z)), vcp_simd_mul(rffac, e1_z_e2_x_minus_e1_x_e2_z));
		M1_z = vcp_simd_add(vcp_simd_mul(p1p2per4pie0r3, vcp_simd_sub(vcp_simd_mul(re2threeperr2, vcp_simd_sub(vcp_simd_mul(eii_x, dy), vcp_simd_mul(eii_y, dx))), e1_x_e2_y_minus_e1_y_e2_x)), vcp_simd_mul(rffac, e1_x_e2_y_minus_e1_y_e2_x));

		M2_x = vcp_simd_sub(vcp_simd_mul(p1p2per4pie0r3, vcp_simd_add(vcp_simd_mul(re1threeperr2, vcp_simd_sub(vcp_simd_mul(ejj_y, dz), vcp_simd_mul(ejj_z, dy))), e1_y_e2_z_minus_e1_z_e2_y)), vcp_simd_mul(rffac, e1_y_e2_z_minus_e1_z_e2_y));
		M2_y = vcp_simd_sub(vcp_simd_mul(p1p2per4pie0r3, vcp_simd_add(vcp_simd_mul(re1threeperr2, vcp_simd_sub(vcp_simd_mul(ejj_z, dx), vcp_simd_mul(ejj_x, dz))), e1_z_e2_x_minus_e1_x_e2_z)), vcp_simd_mul(rffac, e1_z_e2_x_minus_e1_x_e2_z));
		M2_z = vcp_simd_sub(vcp_simd_mul(p1p2per4pie0r3, vcp_simd_add(vcp_simd_mul(re1threeperr2, vcp_simd_sub(vcp_simd_mul(ejj_x, dy), vcp_simd_mul(ejj_y, dx))), e1_x_e2_y_minus_e1_y_e2_x)), vcp_simd_mul(rffac, e1_x_e2_y_minus_e1_y_e2_x));
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

		const vcp_double_vec c_dx = vcp_simd_sub(r1_x, r2_x);
		const vcp_double_vec c_dy = vcp_simd_sub(r1_y, r2_y);
		const vcp_double_vec c_dz = vcp_simd_sub(r1_z, r2_z);

		const vcp_double_vec c_dx2 = vcp_simd_mul(c_dx, c_dx);
		const vcp_double_vec c_dy2 = vcp_simd_mul(c_dy, c_dy);
		const vcp_double_vec c_dz2 = vcp_simd_mul(c_dz, c_dz);

		const vcp_double_vec c_dr2 = vcp_simd_add(vcp_simd_add(c_dx2, c_dy2), c_dz2);

		const vcp_double_vec invdr2_unmasked = vcp_simd_div(one, c_dr2);
		const vcp_double_vec invdr2 = vcp_simd_applymask(invdr2_unmasked, forceMask);
		const vcp_double_vec invdr = vcp_simd_sqrt(invdr2);

		const vcp_double_vec qQ05per4pie0 = vcp_simd_mul(_05, vcp_simd_mul(q, m));

		const vcp_double_vec ejj_xXdx = vcp_simd_mul(ejj_x, c_dx);
		const vcp_double_vec ejj_yXdy = vcp_simd_mul(ejj_y, c_dy);
		const vcp_double_vec ejj_zXdz = vcp_simd_mul(ejj_z, c_dz);
		vcp_double_vec costj = vcp_simd_add(ejj_xXdx, vcp_simd_add(ejj_yXdy, ejj_zXdz));
		costj = vcp_simd_mul(costj, invdr);

		const vcp_double_vec qQinv4dr3 = vcp_simd_mul(qQ05per4pie0, vcp_simd_mul(invdr, invdr2));
		vcp_double_vec part1 = vcp_simd_mul(three, vcp_simd_mul(costj, costj));
		const vcp_double_vec upot = vcp_simd_mul(qQinv4dr3, vcp_simd_sub(part1, one));

		/**********
		 * Force
		 **********/
		const vcp_double_vec minus_partialRijInvdr = vcp_simd_mul(three, vcp_simd_mul(upot, invdr2));
		const vcp_double_vec partialTjInvdr = vcp_simd_mul(vcp_simd_mul(six, costj), vcp_simd_mul(qQinv4dr3, invdr));

		part1 = vcp_simd_mul(costj, vcp_simd_mul(partialTjInvdr, invdr));
		const vcp_double_vec fac = vcp_simd_add(part1, minus_partialRijInvdr);

		f_x = vcp_simd_sub(vcp_simd_mul(fac, c_dx), vcp_simd_mul(partialTjInvdr, ejj_x));
		f_y = vcp_simd_sub(vcp_simd_mul(fac, c_dy), vcp_simd_mul(partialTjInvdr, ejj_y));
		f_z = vcp_simd_sub(vcp_simd_mul(fac, c_dz), vcp_simd_mul(partialTjInvdr, ejj_z));

		const vcp_double_vec m_dx = vcp_simd_sub(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = vcp_simd_sub(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = vcp_simd_sub(m1_r_z, m2_r_z);

		const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);

		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (vcp_simd_movemask(macroMask) > 0) {
			// do we have to mask "upot"? It should already have been masked by "qQinv4dr3" which
			// itself is masked by "invdr2".
			const vcp_double_vec upot_masked = vcp_simd_applymask(upot, macroMask);
			sum_upotXpoles = vcp_simd_add(sum_upotXpoles, upot_masked);

			const vcp_double_vec virial_x = vcp_simd_mul(m_dx, f_x);
			const vcp_double_vec virial_y = vcp_simd_mul(m_dy, f_y);
			const vcp_double_vec virial_z = vcp_simd_mul(m_dz, f_z);

			const vcp_double_vec virial = vcp_simd_add(vcp_simd_add(virial_x, virial_y), virial_z);
			const vcp_double_vec virial_masked = vcp_simd_applymask(virial, macroMask);
			sum_virial = vcp_simd_add(sum_virial, virial_masked);
		}

		/**********
		 * Torque
		 **********/
		const vcp_double_vec minuseXrij_x = vcp_simd_sub(vcp_simd_mul(ejj_z, c_dy), vcp_simd_mul(ejj_y, c_dz));
		const vcp_double_vec minuseXrij_y = vcp_simd_sub(vcp_simd_mul(ejj_x, c_dz), vcp_simd_mul(ejj_z, c_dx));
		const vcp_double_vec minuseXrij_z = vcp_simd_sub(vcp_simd_mul(ejj_y, c_dx), vcp_simd_mul(ejj_x, c_dy));

		M_x = vcp_simd_mul(partialTjInvdr, minuseXrij_x);
		M_y = vcp_simd_mul(partialTjInvdr, minuseXrij_y);
		M_z = vcp_simd_mul(partialTjInvdr, minuseXrij_z);
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


		const vcp_double_vec c_dx = vcp_simd_sub(r1_x, r2_x);
		const vcp_double_vec c_dy = vcp_simd_sub(r1_y, r2_y);
		const vcp_double_vec c_dz = vcp_simd_sub(r1_z, r2_z);

		const vcp_double_vec c_dx2 = vcp_simd_mul(c_dx, c_dx);
		const vcp_double_vec c_dy2 = vcp_simd_mul(c_dy, c_dy);
		const vcp_double_vec c_dz2 = vcp_simd_mul(c_dz, c_dz);

		const vcp_double_vec c_dr2 = vcp_simd_add(vcp_simd_add(c_dx2, c_dy2), c_dz2);

		const vcp_double_vec invdr2_unmasked = vcp_simd_div(one, c_dr2);
		const vcp_double_vec invdr2 = vcp_simd_applymask(invdr2_unmasked, forceMask);
		const vcp_double_vec invdr = vcp_simd_sqrt(invdr2);

		const vcp_double_vec myqfac = vcp_simd_mul(vcp_simd_mul(_1pt5, vcp_simd_mul(p, m)), vcp_simd_mul(invdr2, invdr2));

		const vcp_double_vec eii_xXdx = vcp_simd_mul(eii_x, c_dx);
		const vcp_double_vec eii_yXdy = vcp_simd_mul(eii_y, c_dy);
		const vcp_double_vec eii_zXdz = vcp_simd_mul(eii_z, c_dz);
		vcp_double_vec costi = vcp_simd_add(eii_xXdx, vcp_simd_add(eii_yXdy, eii_zXdz));
		costi = vcp_simd_mul(costi, invdr);

		const vcp_double_vec ejj_xXdx = vcp_simd_mul(ejj_x, c_dx);
		const vcp_double_vec ejj_yXdy = vcp_simd_mul(ejj_y, c_dy);
		const vcp_double_vec ejj_zXdz = vcp_simd_mul(ejj_z, c_dz);
		vcp_double_vec costj = vcp_simd_add(ejj_xXdx, vcp_simd_add(ejj_yXdy, ejj_zXdz));
		costj = vcp_simd_mul(costj, invdr);

		const vcp_double_vec cos2tj = vcp_simd_mul(costj, costj);

		const vcp_double_vec eiiXejj_x = vcp_simd_mul(eii_x, ejj_x);
		const vcp_double_vec eiiXejj_y = vcp_simd_mul(eii_y, ejj_y);
		const vcp_double_vec eiiXejj_z = vcp_simd_mul(eii_z, ejj_z);
		const vcp_double_vec cosgij = vcp_simd_add(eiiXejj_x, vcp_simd_add(eiiXejj_y, eiiXejj_z));

		/************
		 * Potential
		 ************/
		// TODO: Check if upot has to be multiplied by -1 according to DISS_STOLL S.178.
		// This affects also the implementation in potforce.h
		const vcp_double_vec _5cos2tjminus1 = vcp_simd_sub(vcp_simd_mul(five, cos2tj), one);
		const vcp_double_vec _2costj = vcp_simd_mul(two, costj);

		vcp_double_vec part1 = vcp_simd_mul(costi, _5cos2tjminus1);
		vcp_double_vec part2 = vcp_simd_mul(_2costj, cosgij);

		vcp_double_vec const upot = vcp_simd_mul(myqfac, vcp_simd_sub(part2, part1));

		const vcp_double_vec myqfacXinvdr = vcp_simd_mul(myqfac, invdr);
		const vcp_double_vec minus_partialRijInvdr = vcp_simd_mul(four, vcp_simd_mul(upot, invdr2));
		const vcp_double_vec minus_partialTiInvdr = vcp_simd_mul(myqfacXinvdr, _5cos2tjminus1);

		part1 = vcp_simd_mul(five, vcp_simd_mul(costi, costj));
		part1 = vcp_simd_sub(part1, cosgij); // *-1!

		const vcp_double_vec minus_partialTjInvdr = vcp_simd_mul(myqfacXinvdr, vcp_simd_mul(two, part1));
		const vcp_double_vec partialGij = vcp_simd_mul(myqfac, _2costj);

		part1 = vcp_simd_mul(costi, minus_partialTiInvdr);
		part2 = vcp_simd_mul(costj, minus_partialTjInvdr);
		vcp_double_vec part3 = vcp_simd_add(part1, part2);
		const vcp_double_vec fac = vcp_simd_sub(minus_partialRijInvdr, vcp_simd_mul(part3, invdr));

		// Force components
		part1 = vcp_simd_mul(fac, c_dx);
		part2 = vcp_simd_mul(minus_partialTiInvdr, eii_x);
		part3 = vcp_simd_mul(minus_partialTjInvdr, ejj_x);
		f_x = vcp_simd_add(part1, vcp_simd_add(part2, part3));

		part1 = vcp_simd_mul(fac, c_dy);
		part2 = vcp_simd_mul(minus_partialTiInvdr, eii_y);
		part3 = vcp_simd_mul(minus_partialTjInvdr, ejj_y);
		f_y = vcp_simd_add(part1, vcp_simd_add(part2, part3));

		part1 = vcp_simd_mul(fac, c_dz);
		part2 = vcp_simd_mul(minus_partialTiInvdr, eii_z);
		part3 = vcp_simd_mul(minus_partialTjInvdr, ejj_z);
		f_z = vcp_simd_add(part1, vcp_simd_add(part2, part3));

		const vcp_double_vec m_dx = vcp_simd_sub(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = vcp_simd_sub(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = vcp_simd_sub(m1_r_z, m2_r_z);

		const vcp_double_vec macroMask = MacroPolicy::GetMacroMaskSwitched(forceMask, m_dx, m_dy, m_dz, switched);

		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (vcp_simd_movemask(macroMask) > 0) {
			// do we have to mask "upot"? It should already have been masked by "myqfac" which
			// itself is masked by "invdr2".
			const vcp_double_vec upot_masked = vcp_simd_applymask(upot, macroMask);
			sum_upotXpoles = vcp_simd_add(sum_upotXpoles, upot_masked);

			const vcp_double_vec virial_x = vcp_simd_mul(m_dx, f_x);
			const vcp_double_vec virial_y = vcp_simd_mul(m_dy, f_y);
			const vcp_double_vec virial_z = vcp_simd_mul(m_dz, f_z);

			const vcp_double_vec virial = vcp_simd_add(vcp_simd_add(virial_x, virial_y), virial_z);
			const vcp_double_vec virial_masked = vcp_simd_applymask(virial, macroMask);
			sum_virial = vcp_simd_add(sum_virial, virial_masked);
		}

		/**********
		 * Torque
		 **********/
		const vcp_double_vec eii_x_ejj_y = vcp_simd_mul(eii_x, ejj_y);
		const vcp_double_vec eii_x_ejj_z = vcp_simd_mul(eii_x, ejj_z);
		const vcp_double_vec eii_y_ejj_x = vcp_simd_mul(eii_y, ejj_x);
		const vcp_double_vec eii_y_ejj_z = vcp_simd_mul(eii_y, ejj_z);
		const vcp_double_vec eii_z_ejj_x = vcp_simd_mul(eii_z, ejj_x);
		const vcp_double_vec eii_z_ejj_y = vcp_simd_mul(eii_z, ejj_y);

		const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = vcp_simd_sub(eii_x_ejj_y, eii_y_ejj_x);
		const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = vcp_simd_sub(eii_y_ejj_z, eii_z_ejj_y);
		const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = vcp_simd_sub(eii_z_ejj_x, eii_x_ejj_z);

		const vcp_double_vec partialGij_eiXej_x = vcp_simd_mul(partialGij, eii_y_ejj_z_minus_eii_z_ejj_y);
		const vcp_double_vec partialGij_eiXej_y = vcp_simd_mul(partialGij, eii_z_ejj_x_minus_eii_x_ejj_z);
		const vcp_double_vec partialGij_eiXej_z = vcp_simd_mul(partialGij, eii_x_ejj_y_minus_eii_y_ejj_x);

		vcp_double_vec eXrij_x = vcp_simd_sub(vcp_simd_mul(eii_y, c_dz), vcp_simd_mul(eii_z, c_dy));
		vcp_double_vec eXrij_y = vcp_simd_sub(vcp_simd_mul(eii_z, c_dx), vcp_simd_mul(eii_x, c_dz));
		vcp_double_vec eXrij_z = vcp_simd_sub(vcp_simd_mul(eii_x, c_dy), vcp_simd_mul(eii_y, c_dx));

		M1_x = vcp_simd_sub(vcp_simd_mul(minus_partialTiInvdr, eXrij_x), partialGij_eiXej_x);
		M1_y = vcp_simd_sub(vcp_simd_mul(minus_partialTiInvdr, eXrij_y), partialGij_eiXej_y);
		M1_z = vcp_simd_sub(vcp_simd_mul(minus_partialTiInvdr, eXrij_z), partialGij_eiXej_z);

		eXrij_x = vcp_simd_sub(vcp_simd_mul(ejj_y, c_dz), vcp_simd_mul(ejj_z, c_dy));
		eXrij_y = vcp_simd_sub(vcp_simd_mul(ejj_z, c_dx), vcp_simd_mul(ejj_x, c_dz));
		eXrij_z = vcp_simd_sub(vcp_simd_mul(ejj_x, c_dy), vcp_simd_mul(ejj_y, c_dx));

		M2_x = vcp_simd_add(vcp_simd_mul(minus_partialTjInvdr, eXrij_x), partialGij_eiXej_x);
		M2_y = vcp_simd_add(vcp_simd_mul(minus_partialTjInvdr, eXrij_y), partialGij_eiXej_y);
		M2_z = vcp_simd_add(vcp_simd_mul(minus_partialTjInvdr, eXrij_z), partialGij_eiXej_z);
	}

	template<class MacroPolicy>
	inline void _loopBodyQuadrupole(
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
		const vcp_double_vec c_dx = vcp_simd_sub(r1_x, r2_x);
		const vcp_double_vec c_dy = vcp_simd_sub(r1_y, r2_y);
		const vcp_double_vec c_dz = vcp_simd_sub(r1_z, r2_z);

		const vcp_double_vec c_dx2 = vcp_simd_mul(c_dx, c_dx);
		const vcp_double_vec c_dy2 = vcp_simd_mul(c_dy, c_dy);
		const vcp_double_vec c_dz2 = vcp_simd_mul(c_dz, c_dz);

		const vcp_double_vec c_dr2 = vcp_simd_add(vcp_simd_add(c_dx2, c_dy2), c_dz2);

		const vcp_double_vec invdr2_unmasked = vcp_simd_div(one, c_dr2);
		const vcp_double_vec invdr2 = vcp_simd_applymask(invdr2_unmasked, forceMask);
		const vcp_double_vec invdr = vcp_simd_sqrt(invdr2);

		vcp_double_vec qfac = vcp_simd_mul(_075, invdr);
		qfac = vcp_simd_mul(qfac, vcp_simd_mul(mii, mjj));
		qfac = vcp_simd_mul(qfac, vcp_simd_mul(invdr2, invdr2));

		const vcp_double_vec eii_xXdx = vcp_simd_mul(eii_x, c_dx);
		const vcp_double_vec eii_yXdy = vcp_simd_mul(eii_y, c_dy);
		const vcp_double_vec eii_zXdz = vcp_simd_mul(eii_z, c_dz);
		vcp_double_vec costi = vcp_simd_add(eii_xXdx, vcp_simd_add(eii_yXdy, eii_zXdz));
		costi = vcp_simd_mul(costi, invdr);

		const vcp_double_vec ejj_xXdx = vcp_simd_mul(ejj_x, c_dx);
		const vcp_double_vec ejj_yXdy = vcp_simd_mul(ejj_y, c_dy);
		const vcp_double_vec ejj_zXdz = vcp_simd_mul(ejj_z, c_dz);
		vcp_double_vec costj = vcp_simd_add(ejj_xXdx, vcp_simd_add(ejj_yXdy, ejj_zXdz));
		costj = vcp_simd_mul(costj, invdr);

		const vcp_double_vec cos2ti = vcp_simd_mul(costi, costi);
		const vcp_double_vec cos2tj = vcp_simd_mul(costj, costj);

		const vcp_double_vec eiiXejj_x = vcp_simd_mul(eii_x, ejj_x);
		const vcp_double_vec eiiXejj_y = vcp_simd_mul(eii_y, ejj_y);
		const vcp_double_vec eiiXejj_z = vcp_simd_mul(eii_z, ejj_z);
		const vcp_double_vec cosgij = vcp_simd_add(eiiXejj_x, vcp_simd_add(eiiXejj_y, eiiXejj_z));

		vcp_double_vec term = vcp_simd_mul(five, vcp_simd_mul(costi, costj));
		term = vcp_simd_sub(cosgij, term);

		/************
		 * Potential
		 ************/
		vcp_double_vec part1 = vcp_simd_mul(five, vcp_simd_add(cos2ti, cos2tj));
		vcp_double_vec part2 = vcp_simd_mul(_15, vcp_simd_mul(cos2ti, cos2tj));
		vcp_double_vec part3 = vcp_simd_mul(two, vcp_simd_mul(term, term));
		vcp_double_vec upot = vcp_simd_add(part1, part2);
		upot = vcp_simd_sub(vcp_simd_add(one, part3), upot);
		upot = vcp_simd_mul(qfac, upot);

		/**********
		 * Force
		 **********/
		const vcp_double_vec minus_partialRijInvdr = vcp_simd_mul(five, vcp_simd_mul(upot, invdr2));

		// partialTiInvdr & partialTjInvdr
		part1 = vcp_simd_mul(qfac, vcp_simd_mul(ten, invdr));
		part2 = vcp_simd_mul(two, term);

		// partialTiInvdr only
		part3 = vcp_simd_mul(three, vcp_simd_mul(costi, cos2tj));
		vcp_double_vec part4 = vcp_simd_add(costi, vcp_simd_add(part3, vcp_simd_mul(part2, costj)));
		const vcp_double_vec minus_partialTiInvdr = vcp_simd_mul(part1, part4);

		// partialTjInvdr only
		part3 = vcp_simd_mul(three, vcp_simd_mul(costj, cos2ti));
		part4 = vcp_simd_add(costj, vcp_simd_add(part3, vcp_simd_mul(part2, costi)));
		const vcp_double_vec minus_partialTjInvdr = vcp_simd_mul(part1, part4);

		const vcp_double_vec partialGij = vcp_simd_mul(qfac, vcp_simd_mul(four, term));

		// fac
		part1 = vcp_simd_mul(minus_partialTiInvdr, costi);
		part2 = vcp_simd_mul(minus_partialTjInvdr, costj);
		part3 = vcp_simd_mul(vcp_simd_add(part1, part2), invdr);
		const vcp_double_vec fac = vcp_simd_sub(minus_partialRijInvdr, part3);

		// Force components
		part1 = vcp_simd_mul(fac, c_dx);
		part2 = vcp_simd_mul(minus_partialTiInvdr, eii_x);
		part3 = vcp_simd_mul(minus_partialTjInvdr, ejj_x);
		f_x = vcp_simd_add(part1, vcp_simd_add(part2, part3));

		part1 = vcp_simd_mul(fac, c_dy);
		part2 = vcp_simd_mul(minus_partialTiInvdr, eii_y);
		part3 = vcp_simd_mul(minus_partialTjInvdr, ejj_y);
		f_y = vcp_simd_add(part1, vcp_simd_add(part2, part3));

		part1 = vcp_simd_mul(fac, c_dz);
		part2 = vcp_simd_mul(minus_partialTiInvdr, eii_z);
		part3 = vcp_simd_mul(minus_partialTjInvdr, ejj_z);
		f_z = vcp_simd_add(part1, vcp_simd_add(part2, part3));

		const vcp_double_vec m_dx = vcp_simd_sub(m1_r_x, m2_r_x);
		const vcp_double_vec m_dy = vcp_simd_sub(m1_r_y, m2_r_y);
		const vcp_double_vec m_dz = vcp_simd_sub(m1_r_z, m2_r_z);

		const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
		// Check if we have to add the macroscopic values up for at least one of this pairs
		if (vcp_simd_movemask(macroMask) > 0) {
			// do we have to mask "upot"? It should already have been masked by "qfac" ...
			const vcp_double_vec upot_masked = vcp_simd_applymask(upot, macroMask);
			sum_upotXpoles = vcp_simd_add(sum_upotXpoles, upot_masked);

			const vcp_double_vec virial_x = vcp_simd_mul(m_dx, f_x);
			const vcp_double_vec virial_y = vcp_simd_mul(m_dy, f_y);
			const vcp_double_vec virial_z = vcp_simd_mul(m_dz, f_z);

			const vcp_double_vec virial = vcp_simd_add(vcp_simd_add(virial_x, virial_y), virial_z);
			const vcp_double_vec virial_masked = vcp_simd_applymask(virial, macroMask);
			sum_virial = vcp_simd_add(sum_virial, virial_masked);
		}

		const vcp_double_vec eii_x_ejj_y = vcp_simd_mul(eii_x, ejj_y);
		const vcp_double_vec eii_x_ejj_z = vcp_simd_mul(eii_x, ejj_z);
		const vcp_double_vec eii_y_ejj_x = vcp_simd_mul(eii_y, ejj_x);
		const vcp_double_vec eii_y_ejj_z = vcp_simd_mul(eii_y, ejj_z);
		const vcp_double_vec eii_z_ejj_x = vcp_simd_mul(eii_z, ejj_x);
		const vcp_double_vec eii_z_ejj_y = vcp_simd_mul(eii_z, ejj_y);

		const vcp_double_vec eii_x_ejj_y_minus_eii_y_ejj_x = vcp_simd_sub(eii_x_ejj_y, eii_y_ejj_x);
		const vcp_double_vec eii_y_ejj_z_minus_eii_z_ejj_y = vcp_simd_sub(eii_y_ejj_z, eii_z_ejj_y);
		const vcp_double_vec eii_z_ejj_x_minus_eii_x_ejj_z = vcp_simd_sub(eii_z_ejj_x, eii_x_ejj_z);

		const vcp_double_vec partialGij_eiXej_x = vcp_simd_mul(partialGij, eii_y_ejj_z_minus_eii_z_ejj_y);
		const vcp_double_vec partialGij_eiXej_y = vcp_simd_mul(partialGij, eii_z_ejj_x_minus_eii_x_ejj_z);
		const vcp_double_vec partialGij_eiXej_z = vcp_simd_mul(partialGij, eii_x_ejj_y_minus_eii_y_ejj_x);

		vcp_double_vec eXrij_x = vcp_simd_sub(vcp_simd_mul(eii_y, c_dz), vcp_simd_mul(eii_z, c_dy));
		vcp_double_vec eXrij_y = vcp_simd_sub(vcp_simd_mul(eii_z, c_dx), vcp_simd_mul(eii_x, c_dz));
		vcp_double_vec eXrij_z = vcp_simd_sub(vcp_simd_mul(eii_x, c_dy), vcp_simd_mul(eii_y, c_dx));

		Mii_x = vcp_simd_sub(vcp_simd_mul(minus_partialTiInvdr, eXrij_x), partialGij_eiXej_x);
		Mii_y = vcp_simd_sub(vcp_simd_mul(minus_partialTiInvdr, eXrij_y), partialGij_eiXej_y);
		Mii_z = vcp_simd_sub(vcp_simd_mul(minus_partialTiInvdr, eXrij_z), partialGij_eiXej_z);

		eXrij_x = vcp_simd_sub(vcp_simd_mul(ejj_y, c_dz), vcp_simd_mul(ejj_z, c_dy));
		eXrij_y = vcp_simd_sub(vcp_simd_mul(ejj_z, c_dx), vcp_simd_mul(ejj_x, c_dz));
		eXrij_z = vcp_simd_sub(vcp_simd_mul(ejj_x, c_dy), vcp_simd_mul(ejj_y, c_dx));

		Mjj_x = vcp_simd_add(vcp_simd_mul(minus_partialTjInvdr, eXrij_x), partialGij_eiXej_x);
		Mjj_y = vcp_simd_add(vcp_simd_mul(minus_partialTjInvdr, eXrij_y), partialGij_eiXej_y);
		Mjj_z = vcp_simd_add(vcp_simd_mul(minus_partialTjInvdr, eXrij_z), partialGij_eiXej_z);
	}
#endif

#if VCP_VEC_TYPE==VCP_VEC_SSE3 or VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_NOVEC

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
		const vcp_double_vec& eps_24, const vcp_double_vec& sig2,
		const vcp_double_vec& shift6)
{
	const vcp_double_vec c_dx = vcp_simd_sub(r1_x, r2_x);
	const vcp_double_vec c_dy = vcp_simd_sub(r1_y, r2_y);
	const vcp_double_vec c_dz = vcp_simd_sub(r1_z, r2_z);

	const vcp_double_vec c_dxdx = vcp_simd_mul(c_dx, c_dx);
	const vcp_double_vec c_dydy = vcp_simd_mul(c_dy, c_dy);
	const vcp_double_vec c_dzdz = vcp_simd_mul(c_dz, c_dz);
	const vcp_double_vec c_dxdx_dydy = vcp_simd_add(c_dxdx, c_dydy);
	const vcp_double_vec c_r2 = vcp_simd_add(c_dxdx_dydy, c_dzdz);
	const vcp_double_vec r2_inv_unmasked = vcp_simd_div(one, c_r2);
	const vcp_double_vec r2_inv = vcp_simd_applymask(r2_inv_unmasked, forceMask);


	const vcp_double_vec lj2 = vcp_simd_mul(sig2, r2_inv);
	const vcp_double_vec lj4 = vcp_simd_mul(lj2, lj2);
	const vcp_double_vec lj6 = vcp_simd_mul(lj4, lj2);
	const vcp_double_vec lj12 = vcp_simd_mul(lj6, lj6);
	const vcp_double_vec lj12m6 = vcp_simd_sub(lj12, lj6);

	const vcp_double_vec eps24r2inv = vcp_simd_mul(eps_24, r2_inv);
	const vcp_double_vec lj12lj12m6 = vcp_simd_add(lj12, lj12m6);
	const vcp_double_vec scale = vcp_simd_mul(eps24r2inv, lj12lj12m6);

	f_x = vcp_simd_mul(c_dx, scale);
	f_y = vcp_simd_mul(c_dy, scale);
	f_z = vcp_simd_mul(c_dz, scale);

	const vcp_double_vec m_dx = vcp_simd_sub(m1_r_x, m2_r_x);
	const vcp_double_vec m_dy = vcp_simd_sub(m1_r_y, m2_r_y);
	const vcp_double_vec m_dz = vcp_simd_sub(m1_r_z, m2_r_z);

	const vcp_double_vec macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);

	// Only go on if at least 1 macroscopic value has to be calculated.
	if (vcp_simd_movemask(macroMask) > 0) {

		const vcp_double_vec upot = vcp_simd_mul(eps_24, lj12m6);
		const vcp_double_vec upot_sh = vcp_simd_add(shift6, upot);
		const vcp_double_vec upot_masked = vcp_simd_applymask(upot_sh, macroMask);

		sum_upot6lj = vcp_simd_add(sum_upot6lj, upot_masked);

		const vcp_double_vec vir_x = vcp_simd_mul(m_dx, f_x);
		const vcp_double_vec vir_y = vcp_simd_mul(m_dy, f_y);
		const vcp_double_vec vir_z = vcp_simd_mul(m_dz, f_z);

		const vcp_double_vec vir_xy = vcp_simd_add(vir_x, vir_y);
		const vcp_double_vec virial = vcp_simd_add(vir_xy, vir_z);
		const vcp_double_vec vir_masked = vcp_simd_applymask(virial, macroMask);

		sum_virial = vcp_simd_add(sum_virial, vir_masked);
	}
}
#endif

#if VCP_VEC_TYPE==VCP_NOVEC
/**
 * sums up values in a and adds the result to *mem_addr
 */
inline void hSum_Add_Store( double * const mem_addr, const vcp_double_vec & a ) {
	(*mem_addr) += a; //there is just one value of a, so no second sum needed.
}

#elif VCP_VEC_TYPE==VCP_VEC_SSE3
/**
 * sums up values in a and adds the result to *mem_addr
 */
inline void hSum_Add_Store( double * const mem_addr, const vcp_double_vec & a ) {
	_mm_store_sd(
			mem_addr,
			_mm_add_sd(
				_mm_hadd_pd(a, a),
				_mm_load_sd(mem_addr)
			)
	);
}

#elif VCP_VEC_TYPE==VCP_VEC_AVX
static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
/**
 * sums up values in a and adds the result to *mem_addr
 */
inline void hSum_Add_Store( double * const mem_addr, const vcp_double_vec & a ) {
	const vcp_double_vec a_t1 = _mm256_permute2f128_pd(a, a, 0x1);
	const vcp_double_vec a_t2 = _mm256_hadd_pd(a, a_t1);
	const vcp_double_vec a_t3 = _mm256_hadd_pd(a_t2, a_t2);
	_mm256_maskstore_pd(
			mem_addr,
			memoryMask_first,
			vcp_simd_add(
				a_t3,
				vcp_simd_maskload(mem_addr, memoryMask_first)
			)
	);
}

#endif



template<class ForcePolicy>
vcp_doublesizedmask_vec
inline VectorizedCellProcessor::calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		double* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z
		, const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z
		) {

#if VCP_VEC_TYPE==VCP_NOVEC

	unsigned long compute_molecule = 0ul;

	for (size_t j = ForcePolicy :: InitJ(i_center_idx); j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];
		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		const unsigned long forceMask = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? (~0l) : 0l;
		*(soa2_center_dist_lookup + j) = forceMask;
		compute_molecule |= forceMask;

	}

	return compute_molecule;

#elif VCP_VEC_TYPE==VCP_VEC_SSE3

	vcp_double_vec compute_molecule = vcp_simd_zerov();

	// Iterate over centers of second cell
	size_t j = ForcePolicy :: InitJ(i_center_idx);
	for (; j < end_j; j+=VCP_VEC_SIZE) {//end_j is chosen accordingly when function is called. (teilbar durch VCP_VEC_SIZE)
		const vcp_double_vec m2_r_x = vcp_simd_load(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = vcp_simd_load(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = vcp_simd_load(soa2_m_r_z + j);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		const vcp_double_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD);
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
		const unsigned long forceMask_l = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;
		// this casting via void* is required for gcc
		const void* forceMask_tmp = reinterpret_cast<const void*>(&forceMask_l);
		double forceMask = *reinterpret_cast<double const* /*const*/>(forceMask_tmp);

		*(soa2_center_dist_lookup + j) = forceMask;
		const vcp_double_vec forceMask_vec = vcp_simd_set1(forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask_vec);
	}
	return compute_molecule;

#elif VCP_VEC_TYPE==VCP_VEC_AVX

	vcp_double_vec compute_molecule = vcp_simd_zerov();

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	vcp_double_vec initJ_mask = ForcePolicy::InitJ_Mask(i_center_idx);
	// Iterate over centers of second cell
	for (; j < end_j; j+=VCP_VEC_SIZE) {//end_j is chosen accordingly when function is called. (teilbar durch VCP_VEC_SIZE)
		const vcp_double_vec m2_r_x = vcp_simd_load(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = vcp_simd_load(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = vcp_simd_load(soa2_m_r_z + j);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		const vcp_double_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);
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
		unsigned long forceMask_l;
			// DetectSingleCell() = false for SingleCellDistinctPolicy and CellPairPolicy, true for SingleCellPolicy
			// we need this, since in contrast to sse3 we can no longer guarantee, that j>=i by default (j==i is handled by ForcePolicy::Condition).
			// however only one of the branches should be chosen by the compiler, since the class is known at compile time.
		if (ForcePolicy::DetectSingleCell()) {
			forceMask_l = (ForcePolicy::Condition(m_r2, cutoffRadiusSquare) && j > i_center_idx) ? ~0l : 0l;
		} else {
			forceMask_l = ForcePolicy::Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;
		}

		//this casting via void* is required for gcc
		void* forceMask_tmp = reinterpret_cast<void*>(&forceMask_l);
		double forceMask = *reinterpret_cast<double const* /*const*/>(forceMask_tmp);

		*(soa2_center_dist_lookup + j) = forceMask;
		const vcp_double_vec forceMask_vec = vcp_simd_set1(forceMask);
		compute_molecule = vcp_simd_or(compute_molecule, forceMask_vec);
	}
	return compute_molecule;

#endif
}

template<class ForcePolicy, class MacroPolicy>
void VectorizedCellProcessor :: _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
#if VCP_VEC_TYPE==VCP_NOVEC-1
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
		unsigned long compute_molecule_lj = calcDistLookup<ForcePolicy>(soa1, i, i_ljc_idx, soa2._ljc_num, _LJCutoffRadiusSquare,
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




	static const vcp_double_vec ones = vcp_simd_ones();

	vcp_double_vec sum_upot6lj = vcp_simd_zerov();
	vcp_double_vec sum_upotXpoles = vcp_simd_zerov();
	vcp_double_vec sum_virial = vcp_simd_zerov();
	vcp_double_vec sum_myRF = vcp_simd_zerov();

	const vcp_double_vec rc2 = vcp_simd_set1(_LJCutoffRadiusSquare);
	const vcp_double_vec cutoffRadiusSquare = vcp_simd_set1(_cutoffRadiusSquare);
	const vcp_double_vec epsRFInvrc3 = vcp_simd_broadcast(&_epsRFInvrc3);

	/*
	 *  "var" & (~VCP_VEC_SIZE_M1):
	 *  Making sure that result is a multiple of VCP_VEC_SIZE. If "var" can not be divided by VCP_VEC_SIZE
	 *  the result is: the next smaller multiple of VCP_VEC_SIZE.
	 */
	const size_t end_ljc_j = soa2._ljc_num & (~VCP_VEC_SIZE_M1);
	const size_t end_ljc_j_longloop = (soa2._ljc_num + VCP_VEC_SIZE_M1) & (~VCP_VEC_SIZE_M1);//this is ceil _ljc_num, VCP_VEC_SIZE
	const size_t end_charges_j = soa2._charges_num & (~VCP_VEC_SIZE_M1);
	const size_t end_charges_j_longloop = (soa2._charges_num + VCP_VEC_SIZE_M1) & (~VCP_VEC_SIZE_M1);//this is ceil _charges_num, VCP_VEC_SIZE
	const size_t end_dipoles_j = soa2._dipoles_num & (~VCP_VEC_SIZE_M1);
	const size_t end_dipoles_j_longloop = (soa2._dipoles_num + VCP_VEC_SIZE_M1) & (~VCP_VEC_SIZE_M1);//this is ceil _dipoles_num, VCP_VEC_SIZE
	const size_t end_quadrupoles_j = soa2._quadrupoles_num & (~VCP_VEC_SIZE_M1);
	const size_t end_quadrupoles_j_longloop = (soa2._quadrupoles_num + VCP_VEC_SIZE_M1) & (~VCP_VEC_SIZE_M1);//this is ceil _quadrupoles_num, VCP_VEC_SIZE

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
		const vcp_double_vec m1_r_x = vcp_simd_broadcast(soa1_mol_pos_x + i);
		const vcp_double_vec m1_r_y = vcp_simd_broadcast(soa1_mol_pos_y + i);
		const vcp_double_vec m1_r_z = vcp_simd_broadcast(soa1_mol_pos_z + i);

		// Iterate over centers of second cell
		const vcp_doublesizedmask_vec compute_molecule_ljc = calcDistLookup<ForcePolicy>(soa1, i, i_ljc_idx, soa2._ljc_num, _LJCutoffRadiusSquare,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_doublesizedmask_vec compute_molecule_charges = calcDistLookup<ForcePolicy>(soa1, i, i_charge_idx, soa2._charges_num, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_doublesizedmask_vec compute_molecule_dipoles = calcDistLookup<ForcePolicy>(soa1, i, i_dipole_idx, soa2._dipoles_num, _cutoffRadiusSquare,
				soa2_dipoles_dist_lookup, soa2_dipoles_m_r_x, soa2_dipoles_m_r_y, soa2_dipoles_m_r_z,
				cutoffRadiusSquare,	end_dipoles_j, m1_r_x, m1_r_y, m1_r_z);
		const vcp_doublesizedmask_vec compute_molecule_quadrupoles = calcDistLookup<ForcePolicy>(soa1, i, i_quadrupole_idx, soa2._quadrupoles_num, _cutoffRadiusSquare,
				soa2_quadrupoles_dist_lookup, soa2_quadrupoles_m_r_x, soa2_quadrupoles_m_r_y, soa2_quadrupoles_m_r_z,
				cutoffRadiusSquare, end_quadrupoles_j, m1_r_x, m1_r_y, m1_r_z);

		if (!vcp_simd_movemask(compute_molecule_ljc)) {
			i_ljc_idx += soa1_mol_ljc_num[i];
		}
		else {
			// LJ force computation
			for (int local_i = 0; local_i < soa1_mol_ljc_num[i]; local_i++) {
				vcp_double_vec sum_fx1 = vcp_simd_zerov();
				vcp_double_vec sum_fy1 = vcp_simd_zerov();
				vcp_double_vec sum_fz1 = vcp_simd_zerov();
				const vcp_double_vec c_r_x1 = vcp_simd_broadcast(soa1_ljc_r_x + i_ljc_idx);
				const vcp_double_vec c_r_y1 = vcp_simd_broadcast(soa1_ljc_r_y + i_ljc_idx);
				const vcp_double_vec c_r_z1 = vcp_simd_broadcast(soa1_ljc_r_z + i_ljc_idx);
				// Iterate over each pair of centers in the second cell.
				size_t j = ForcePolicy::InitJ(i_ljc_idx);
				for (; j < end_ljc_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_ljc_dist_lookup + j);
					// Only go on if at least 1 of the forces has to be calculated.
					if (vcp_simd_movemask(forceMask) > 0) {
						const vcp_double_vec c_r_x2 = vcp_simd_load(soa2_ljc_r_x + j);
						const vcp_double_vec c_r_y2 = vcp_simd_load(soa2_ljc_r_y + j);
						const vcp_double_vec c_r_z2 = vcp_simd_load(soa2_ljc_r_z + j);

						const vcp_double_vec m_r_x2 = vcp_simd_load(soa2_ljc_m_r_x + j);
						const vcp_double_vec m_r_y2 = vcp_simd_load(soa2_ljc_m_r_y + j);
						const vcp_double_vec m_r_z2 = vcp_simd_load(soa2_ljc_m_r_z + j);

						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						vcp_double_vec fx, fy, fz;

						vcp_double_vec eps_24;
						vcp_double_vec sig2;
						unpackEps24Sig2(eps_24, sig2, _eps_sig[id_i], soa2_ljc_id + j);

						vcp_double_vec shift6;
						unpackShift6(shift6, _shift6[id_i], soa2_ljc_id + j);

						_loopBodyLJ<MacroPolicy>(
							m1_r_x, m1_r_y, m1_r_z, c_r_x1, c_r_y1, c_r_z1,
							m_r_x2, m_r_y2, m_r_z2, c_r_x2, c_r_y2, c_r_z2,
							fx, fy, fz,
							sum_upot6lj, sum_virial,
							forceMask,
							eps_24, sig2,
							shift6);
						const vcp_double_vec old_fx2 = vcp_simd_load(soa2_ljc_f_x + j);
						const vcp_double_vec new_fx2 = vcp_simd_sub(old_fx2, fx);
						vcp_simd_store(soa2_ljc_f_x + j, new_fx2);
						const vcp_double_vec old_fy2 = vcp_simd_load(soa2_ljc_f_y + j);
						const vcp_double_vec new_fy2 = vcp_simd_sub(old_fy2, fy);
						vcp_simd_store(soa2_ljc_f_y + j, new_fy2);
						const vcp_double_vec old_fz2 = vcp_simd_load(soa2_ljc_f_z + j);
						const vcp_double_vec new_fz2 = vcp_simd_sub(old_fz2, fz);
						vcp_simd_store(soa2_ljc_f_z + j, new_fz2);
						sum_fx1 = vcp_simd_add(sum_fx1, fx);
						sum_fy1 = vcp_simd_add(sum_fy1, fy);
						sum_fz1 = vcp_simd_add(sum_fz1, fz);
					}
				}

				hSum_Add_Store(soa1_ljc_f_x + i_ljc_idx, sum_fx1);
				hSum_Add_Store(soa1_ljc_f_y + i_ljc_idx, sum_fy1);
				hSum_Add_Store(soa1_ljc_f_z + i_ljc_idx, sum_fz1);

				i_ljc_idx++;
			}
		}

		// Computation of site interactions with charges

		if (!vcp_simd_movemask(compute_molecule_charges)) {
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

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_charge_idx + local_i);
				for (; j < end_charges_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {
						const vcp_double_vec q2 = vcp_simd_load(soa2_charges_q + j);

						const vcp_double_vec r2_x = vcp_simd_load(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = vcp_simd_load(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = vcp_simd_load(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z;
						_loopBodyCharge<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								sum_upotXpoles, sum_virial,
								forceMask);



						sum_f1_x = vcp_simd_add(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_add(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_add(sum_f1_z, f_z);

						vcp_double_vec f2_x = vcp_simd_load(soa2_charges_f_x + j);
						vcp_double_vec f2_y = vcp_simd_load(soa2_charges_f_y + j);
						vcp_double_vec f2_z = vcp_simd_load(soa2_charges_f_z + j);

						f2_x = vcp_simd_sub(f2_x, f_x);
						f2_y = vcp_simd_sub(f2_y, f_y);
						f2_z = vcp_simd_sub(f2_z, f_z);

						vcp_simd_store(soa2_charges_f_x + j, f2_x);
						vcp_simd_store(soa2_charges_f_y + j, f2_y);
						vcp_simd_store(soa2_charges_f_z + j, f2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx + local_i, sum_f1_z);

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

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M_x = vcp_simd_zerov();
				vcp_double_vec sum_M_y = vcp_simd_zerov();
				vcp_double_vec sum_M_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {

						const vcp_double_vec q = vcp_simd_load(soa2_charges_q + j);

						const vcp_double_vec r2_x = vcp_simd_load(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = vcp_simd_load(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = vcp_simd_load(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeDipole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = vcp_simd_sub(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_sub(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_sub(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_charges_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_charges_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_charges_f_z + j);

						sum_f2_x = vcp_simd_add(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_add(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_add(sum_f2_z, f_z);

						vcp_simd_store(soa2_charges_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_charges_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_charges_f_z + j, sum_f2_z);

						// Store torque

						sum_M_x = vcp_simd_add(sum_M_x, M_x);
						sum_M_y = vcp_simd_add(sum_M_y, M_y);
						sum_M_z = vcp_simd_add(sum_M_z, M_z);
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

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {

						const vcp_double_vec q = vcp_simd_load(soa2_charges_q + j);

						const vcp_double_vec r2_x = vcp_simd_load(soa2_charges_r_x + j);
						const vcp_double_vec r2_y = vcp_simd_load(soa2_charges_r_y + j);
						const vcp_double_vec r2_z = vcp_simd_load(soa2_charges_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_charges_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_charges_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_charges_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeQuadrupole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = vcp_simd_sub(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_sub(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_sub(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_charges_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_charges_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_charges_f_z + j);

						sum_f2_x = vcp_simd_add(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_add(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_add(sum_f2_z, f_z);

						vcp_simd_store(soa2_charges_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_charges_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_charges_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M_z);
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

				i_quadrupole_charge_idx++;
			}

			i_charge_idx += soa1_mol_charges_num[i];
		}

		// Computation of site interactions with dipoles

		// Continue with next molecule if no force has to be calculated
		if (!vcp_simd_movemask(compute_molecule_dipoles)) {
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

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_dipole_idx + local_i);
				for (; j < end_dipoles_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {
						const vcp_double_vec p2 = vcp_simd_load(soa2_dipoles_p + j);
						const vcp_double_vec e2_x = vcp_simd_load(soa2_dipoles_e_x + j);
						const vcp_double_vec e2_y = vcp_simd_load(soa2_dipoles_e_y + j);
						const vcp_double_vec e2_z = vcp_simd_load(soa2_dipoles_e_z + j);
						const vcp_double_vec r2_x = vcp_simd_load(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = vcp_simd_load(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = vcp_simd_load(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_dipoles_m_r_z + j);

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

						sum_f1_x = vcp_simd_add(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_add(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_add(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_dipoles_f_z + j);

						sum_f2_x = vcp_simd_sub(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_sub(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_sub(sum_f2_z, f_z);

						vcp_simd_store(soa2_dipoles_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_dipoles_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = vcp_simd_load(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = vcp_simd_load(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = vcp_simd_load(soa2_dipoles_M_z + j);

						sum_M2_x = vcp_simd_add(sum_M2_x, M2_x);
						sum_M2_y = vcp_simd_add(sum_M2_y, M2_y);
						sum_M2_z = vcp_simd_add(sum_M2_z, M2_z);

						vcp_simd_store(soa2_dipoles_M_x + j, sum_M2_x);
						vcp_simd_store(soa2_dipoles_M_y + j, sum_M2_y);
						vcp_simd_store(soa2_dipoles_M_z + j, sum_M2_z);
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

			}

			// Computation of charge-dipole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{

				const vcp_double_vec q = vcp_simd_broadcast(soa1_charges_q + i_charge_dipole_idx);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_charges_r_x + i_charge_dipole_idx);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_charges_r_y + i_charge_dipole_idx);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_charges_r_z + i_charge_dipole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_dipole_idx);
				for (; j < end_dipoles_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {

						const vcp_double_vec p = vcp_simd_load(soa2_dipoles_p + j);

						const vcp_double_vec e_x = vcp_simd_load(soa2_dipoles_e_x + j);
						const vcp_double_vec e_y = vcp_simd_load(soa2_dipoles_e_y + j);
						const vcp_double_vec e_z = vcp_simd_load(soa2_dipoles_e_z + j);

						const vcp_double_vec r2_x = vcp_simd_load(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = vcp_simd_load(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = vcp_simd_load(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeDipole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e_x, e_y, e_z, p,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = vcp_simd_add(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_add(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_add(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_dipoles_f_z + j);

						sum_f2_x = vcp_simd_sub(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_sub(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_sub(sum_f2_z, f_z);

						vcp_simd_store(soa2_dipoles_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_dipoles_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						vcp_double_vec sum_M2_x = vcp_simd_load(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = vcp_simd_load(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = vcp_simd_load(soa2_dipoles_M_z + j);

						sum_M2_x = vcp_simd_add(sum_M2_x, M_x);
						sum_M2_y = vcp_simd_add(sum_M2_y, M_y);
						sum_M2_z = vcp_simd_add(sum_M2_z, M_z);

						vcp_simd_store(soa2_dipoles_M_x + j, sum_M2_x);
						vcp_simd_store(soa2_dipoles_M_y + j, sum_M2_y);
						vcp_simd_store(soa2_dipoles_M_z + j, sum_M2_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_dipole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_dipole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_dipole_idx, sum_f1_z);

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

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_dipole_idx);
				for (; j < end_dipoles_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_dipoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {
						const vcp_double_vec p = vcp_simd_load(soa2_dipoles_p + j);
						const vcp_double_vec e2_x = vcp_simd_load(soa2_dipoles_e_x + j);
						const vcp_double_vec e2_y = vcp_simd_load(soa2_dipoles_e_y + j);
						const vcp_double_vec e2_z = vcp_simd_load(soa2_dipoles_e_z + j);
						const vcp_double_vec r2_x = vcp_simd_load(soa2_dipoles_r_x + j);
						const vcp_double_vec r2_y = vcp_simd_load(soa2_dipoles_r_y + j);
						const vcp_double_vec r2_z = vcp_simd_load(soa2_dipoles_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_dipoles_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_dipoles_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_dipoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipoleQuadrupole<MacroPolicy>(
								m2_r_x, m2_r_y, m2_r_z,	r2_x, r2_y, r2_z, e2_x, e2_y, e2_z, p,
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, e1_x, e1_y, e1_z, m,
								f_x, f_y, f_z, M2_x, M2_y, M2_z, M1_x, M1_y, M1_z,
								sum_upotXpoles, sum_virial,
								forceMask, ones);

						// Store forces

						sum_f1_x = vcp_simd_sub(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_sub(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_sub(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_dipoles_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_dipoles_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_dipoles_f_z + j);

						sum_f2_x = vcp_simd_add(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_add(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_add(sum_f2_z, f_z);

						vcp_simd_store(soa2_dipoles_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_dipoles_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_dipoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = vcp_simd_load(soa2_dipoles_M_x + j);
						vcp_double_vec sum_M2_y = vcp_simd_load(soa2_dipoles_M_y + j);
						vcp_double_vec sum_M2_z = vcp_simd_load(soa2_dipoles_M_z + j);

						sum_M2_x = vcp_simd_add(sum_M2_x, M2_x);
						sum_M2_y = vcp_simd_add(sum_M2_y, M2_y);
						sum_M2_z = vcp_simd_add(sum_M2_z, M2_z);

						vcp_simd_store(soa2_dipoles_M_x + j, sum_M2_x);
						vcp_simd_store(soa2_dipoles_M_y + j, sum_M2_y);
						vcp_simd_store(soa2_dipoles_M_z + j, sum_M2_z);
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

				i_quadrupole_dipole_idx++;

			}

			i_dipole_idx += soa1_mol_dipoles_num[i];
		}

		// Computation of site interactions with quadrupoles

		if (!vcp_simd_movemask(compute_molecule_quadrupoles)) {
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

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_quadrupole_idx + local_i);
				for (; j < end_quadrupoles_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {

						const vcp_double_vec mjj = vcp_simd_load(soa2_quadrupoles_m + j);
						const vcp_double_vec ejj_x = vcp_simd_load(soa2_quadrupoles_e_x + j);
						const vcp_double_vec ejj_y = vcp_simd_load(soa2_quadrupoles_e_y + j);
						const vcp_double_vec ejj_z = vcp_simd_load(soa2_quadrupoles_e_z + j);
						const vcp_double_vec rjj_x = vcp_simd_load(soa2_quadrupoles_r_x + j);
						const vcp_double_vec rjj_y = vcp_simd_load(soa2_quadrupoles_r_y + j);
						const vcp_double_vec rjj_z = vcp_simd_load(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_quadrupoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyQuadrupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, mii,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, mjj,
								f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								forceMask);

						// Store forces

						sum_f1_x = vcp_simd_add(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_add(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_add(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_quadrupoles_f_z + j);

						sum_f2_x = vcp_simd_sub(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_sub(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_sub(sum_f2_z, f_z);

						vcp_simd_store(soa2_quadrupoles_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_quadrupoles_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = vcp_simd_load(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M2_y = vcp_simd_load(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M2_z = vcp_simd_load(soa2_quadrupoles_M_z + j);

						sum_M2_x = vcp_simd_add(sum_M2_x, M2_x);
						sum_M2_y = vcp_simd_add(sum_M2_y, M2_y);
						sum_M2_z = vcp_simd_add(sum_M2_z, M2_z);

						vcp_simd_store(soa2_quadrupoles_M_x + j, sum_M2_x);
						vcp_simd_store(soa2_quadrupoles_M_y + j, sum_M2_y);
						vcp_simd_store(soa2_quadrupoles_M_z + j, sum_M2_z);
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

			}

			// Computation of charge-quadrupole interactions

			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++)
			{
				const vcp_double_vec q = vcp_simd_broadcast(soa1_charges_q + i_charge_quadrupole_idx);
				const vcp_double_vec r1_x = vcp_simd_broadcast(soa1_charges_r_x + i_charge_quadrupole_idx);
				const vcp_double_vec r1_y = vcp_simd_broadcast(soa1_charges_r_y + i_charge_quadrupole_idx);
				const vcp_double_vec r1_z = vcp_simd_broadcast(soa1_charges_r_z + i_charge_quadrupole_idx);

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				size_t j = ForcePolicy::InitJ(i_quadrupole_idx);
				for (; j < end_quadrupoles_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {

						const vcp_double_vec m = vcp_simd_load(soa2_quadrupoles_m + j);
						const vcp_double_vec e_x = vcp_simd_load(soa2_quadrupoles_e_x + j);
						const vcp_double_vec e_y = vcp_simd_load(soa2_quadrupoles_e_y + j);
						const vcp_double_vec e_z = vcp_simd_load(soa2_quadrupoles_e_z + j);
						const vcp_double_vec r2_x = vcp_simd_load(soa2_quadrupoles_r_x + j);
						const vcp_double_vec r2_y = vcp_simd_load(soa2_quadrupoles_r_y + j);
						const vcp_double_vec r2_z = vcp_simd_load(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_quadrupoles_m_r_z + j);


						vcp_double_vec f_x, f_y, f_z, M_x, M_y, M_z;

						_loopBodyChargeQuadrupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, e_x, e_y, e_z, m,
								f_x, f_y, f_z, M_x, M_y, M_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = vcp_simd_add(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_add(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_add(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_quadrupoles_f_z + j);

						sum_f2_x = vcp_simd_sub(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_sub(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_sub(sum_f2_z, f_z);

						vcp_simd_store(soa2_quadrupoles_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_quadrupoles_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque
						vcp_double_vec sum_M_x = vcp_simd_load(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M_y = vcp_simd_load(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M_z = vcp_simd_load(soa2_quadrupoles_M_z + j);

						sum_M_x = vcp_simd_add(sum_M_x, M_x);
						sum_M_y = vcp_simd_add(sum_M_y, M_y);
						sum_M_z = vcp_simd_add(sum_M_z, M_z);

						vcp_simd_store(soa2_quadrupoles_M_x + j, sum_M_x);
						vcp_simd_store(soa2_quadrupoles_M_y + j, sum_M_y);
						vcp_simd_store(soa2_quadrupoles_M_z + j, sum_M_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_quadrupole_idx, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_quadrupole_idx, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_quadrupole_idx, sum_f1_z);

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

				vcp_double_vec sum_f1_x = vcp_simd_zerov();
				vcp_double_vec sum_f1_y = vcp_simd_zerov();
				vcp_double_vec sum_f1_z = vcp_simd_zerov();

				vcp_double_vec sum_M1_x = vcp_simd_zerov();
				vcp_double_vec sum_M1_y = vcp_simd_zerov();
				vcp_double_vec sum_M1_z = vcp_simd_zerov();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_quadrupole_idx);
				for (; j < end_quadrupoles_j_longloop; j += VCP_VEC_SIZE) {
					const vcp_double_vec forceMask = vcp_simd_load(soa2_quadrupoles_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (vcp_simd_movemask(forceMask) > 0) {

						const vcp_double_vec m = vcp_simd_load(soa2_quadrupoles_m + j);
						const vcp_double_vec ejj_x = vcp_simd_load(soa2_quadrupoles_e_x + j);
						const vcp_double_vec ejj_y = vcp_simd_load(soa2_quadrupoles_e_y + j);
						const vcp_double_vec ejj_z = vcp_simd_load(soa2_quadrupoles_e_z + j);
						const vcp_double_vec rjj_x = vcp_simd_load(soa2_quadrupoles_r_x + j);
						const vcp_double_vec rjj_y = vcp_simd_load(soa2_quadrupoles_r_y + j);
						const vcp_double_vec rjj_z = vcp_simd_load(soa2_quadrupoles_r_z + j);

						const vcp_double_vec m2_r_x = vcp_simd_load(soa2_quadrupoles_m_r_x + j);
						const vcp_double_vec m2_r_y = vcp_simd_load(soa2_quadrupoles_m_r_y + j);
						const vcp_double_vec m2_r_z = vcp_simd_load(soa2_quadrupoles_m_r_z + j);

						vcp_double_vec f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z;

						_loopBodyDipoleQuadrupole<MacroPolicy>(
								m1_r_x, m1_r_y, m1_r_z, rii_x, rii_y, rii_z, eii_x, eii_y, eii_z, p,
								m2_r_x, m2_r_y, m2_r_z,	rjj_x, rjj_y, rjj_z, ejj_x, ejj_y, ejj_z, m,
								f_x, f_y, f_z, M1_x, M1_y, M1_z, M2_x, M2_y, M2_z,
								sum_upotXpoles, sum_virial,
								forceMask, zero);

						// Store forces

						sum_f1_x = vcp_simd_add(sum_f1_x, f_x);
						sum_f1_y = vcp_simd_add(sum_f1_y, f_y);
						sum_f1_z = vcp_simd_add(sum_f1_z, f_z);

						vcp_double_vec sum_f2_x = vcp_simd_load(soa2_quadrupoles_f_x + j);
						vcp_double_vec sum_f2_y = vcp_simd_load(soa2_quadrupoles_f_y + j);
						vcp_double_vec sum_f2_z = vcp_simd_load(soa2_quadrupoles_f_z + j);

						sum_f2_x = vcp_simd_sub(sum_f2_x, f_x);
						sum_f2_y = vcp_simd_sub(sum_f2_y, f_y);
						sum_f2_z = vcp_simd_sub(sum_f2_z, f_z);

						vcp_simd_store(soa2_quadrupoles_f_x + j, sum_f2_x);
						vcp_simd_store(soa2_quadrupoles_f_y + j, sum_f2_y);
						vcp_simd_store(soa2_quadrupoles_f_z + j, sum_f2_z);

						// Store torque

						sum_M1_x = vcp_simd_add(sum_M1_x, M1_x);
						sum_M1_y = vcp_simd_add(sum_M1_y, M1_y);
						sum_M1_z = vcp_simd_add(sum_M1_z, M1_z);

						vcp_double_vec sum_M2_x = vcp_simd_load(soa2_quadrupoles_M_x + j);
						vcp_double_vec sum_M2_y = vcp_simd_load(soa2_quadrupoles_M_y + j);
						vcp_double_vec sum_M2_z = vcp_simd_load(soa2_quadrupoles_M_z + j);

						sum_M2_x = vcp_simd_add(sum_M2_x, M2_x);
						sum_M2_y = vcp_simd_add(sum_M2_y, M2_y);
						sum_M2_z = vcp_simd_add(sum_M2_z, M2_z);

						vcp_simd_store(soa2_quadrupoles_M_x + j, sum_M2_x);
						vcp_simd_store(soa2_quadrupoles_M_y + j, sum_M2_y);
						vcp_simd_store(soa2_quadrupoles_M_z + j, sum_M2_z);
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

				i_dipole_quadrupole_idx++;

			}

			i_quadrupole_idx += soa1_mol_quadrupoles_num[i];
		}
	}

	hSum_Add_Store(&_upot6lj, sum_upot6lj);
	hSum_Add_Store(&_upotXpoles, sum_upotXpoles);
	hSum_Add_Store(&_virial, sum_virial);
	hSum_Add_Store(&_myRF, vcp_simd_sub(zero, sum_myRF));

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
