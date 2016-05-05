#include "particleContainer/ParticleCell.h"
#include "molecules/Molecule.h"
#include "utils/UnorderedVector.h"

#include <cassert>
#include <vector>



using namespace std;

ParticleCell::ParticleCell() :
		_molecules(), _cellDataSoA(0, 0, 0, 0, 0), _cellIndex(0) {
	for (int d = 0; d < 3; ++d) {
		_boxMin[d] = 0.0;
		_boxMax[d] = 0.0;
	}
}

ParticleCell::~ParticleCell() {
//	if(!isEmpty()) {
//		deallocateAllParticles();
//	}
}

void ParticleCell::removeAllParticles() {
	_molecules.clear();
}

void ParticleCell::deallocateAllParticles() {
	std::vector<Molecule * >::iterator it;
	for(it = _molecules.begin(); it != _molecules.end(); ++it) {
		delete *it;
	}
	removeAllParticles();
}

bool ParticleCell::addParticle(Molecule* particle_ptr, bool checkWhetherDuplicate) {
	bool wasInserted = false;

#ifndef NDEBUG
	bool isIn = particle_ptr->inBox(_boxMin, _boxMax);
//	assert(isIn);
#endif
	if (checkWhetherDuplicate == false) {
		_molecules.push_back(particle_ptr);
		wasInserted = true;
	} else {

		// perform a check whether this molecule exists (has been received) already

		vector<Molecule*>::const_iterator it;
		unsigned long pID = particle_ptr->id();

		bool found = false;

		for (it = _molecules.begin(); it != _molecules.end(); ++it) {
			if (pID == (*it)->id()) {
				found = true;
				break;
			}
		}

		if (not found) {
			_molecules.push_back(particle_ptr);
			wasInserted = true;
		}
	}
	return wasInserted;
}

vector<Molecule*>& ParticleCell::getParticlePointers() {
	return _molecules;
}

const vector<Molecule*>& ParticleCell::getParticlePointers() const {
	return _molecules;
}

bool ParticleCell::isEmpty() const {
	return _molecules.empty();
}

int ParticleCell::getMoleculeCount() const {
	return _molecules.size();
}

bool ParticleCell::deleteMolecule(unsigned long molecule_id) {
	bool found = false;
	vector<Molecule*>::iterator molecule_iter;

	for (molecule_iter = _molecules.begin(); molecule_iter != _molecules.end(); molecule_iter++) {
		Molecule *molecule = *molecule_iter;
		if (molecule->id() == molecule_id) {
			found = true;
			UnorderedVector::fastRemove(_molecules, molecule_iter);
			break;
		}
	}
	return found;
}

void ParticleCell::fastRemoveMolecule(std::vector<Molecule *>::iterator& it) {
	UnorderedVector::fastRemove(_molecules, it);
}

std::vector<Molecule*>& ParticleCell::filterLeavingMolecules() {
	std::vector<Molecule*>::iterator it;

	for (it = _molecules.begin(); it != _molecules.end();) {
		(*it)->setSoA(nullptr);
		const Molecule * const mol = *it;

		bool isStaying = mol->inBox(_boxMin, _boxMax);

		if (isStaying) {
			++it; // next iteration
		} else {
			_leavingMolecules.push_back(*it);

			// erase from _molecules:
			*it = _molecules.back(); // iterator points to a new molecule
			_molecules.back() = NULL;
			_molecules.pop_back();
			// don't advance iterator, it now points to a new molecule
		}
	}

	return _leavingMolecules;
}

void ParticleCell::getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer) {
	std::vector<Molecule *>::iterator particleIter;

	for (particleIter = _molecules.begin(); particleIter != _molecules.end();) {
		if ((*particleIter)->inBox(lowCorner, highCorner)) {
			particlePtrs.push_back(*particleIter);
			if (removeFromContainer) {
				UnorderedVector::fastRemove(_molecules, particleIter);
				// particleIter already points at next molecule, so continue without incrementing
				continue;
			}
		}
		++particleIter;
	}
}

void ParticleCell::buildSoACaches(const std::vector<size_t>& compIDs) {

	// Determine the total number of centers.
	size_t numMolecules = _molecules.size();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;

	for (size_t m = 0;  m < numMolecules; ++m) {
		nLJCenters += _molecules[m]->numLJcenters();
		nCharges += _molecules[m]->numCharges();
		nDipoles += _molecules[m]->numDipoles();
		nQuadrupoles += _molecules[m]->numQuadrupoles();
	}

	// Construct the SoA.
	CellDataSoA & soa = _cellDataSoA;
	soa.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);

	double* const soa_ljc_m_r_x = soa.ljc_m_r_xBegin();
	double* const soa_ljc_m_r_y = soa.ljc_m_r_yBegin();
	double* const soa_ljc_m_r_z = soa.ljc_m_r_zBegin();
	double* const soa_ljc_r_x = soa.ljc_r_xBegin();
	double* const soa_ljc_r_y = soa.ljc_r_yBegin();
	double* const soa_ljc_r_z = soa.ljc_r_zBegin();

	double* const soa_charges_m_r_x = soa.charges_m_r_xBegin();
	double* const soa_charges_m_r_y = soa.charges_m_r_yBegin();
	double* const soa_charges_m_r_z = soa.charges_m_r_zBegin();
	double* const soa_charges_r_x = soa.charges_r_xBegin();
	double* const soa_charges_r_y = soa.charges_r_yBegin();
	double* const soa_charges_r_z = soa.charges_r_zBegin();

	double* const soa_dipoles_m_r_x = soa.dipoles_m_r_xBegin();
	double* const soa_dipoles_m_r_y = soa.dipoles_m_r_yBegin();
	double* const soa_dipoles_m_r_z = soa.dipoles_m_r_zBegin();
	double* const soa_dipoles_r_x = soa.dipoles_r_xBegin();
	double* const soa_dipoles_r_y = soa.dipoles_r_yBegin();
	double* const soa_dipoles_r_z = soa.dipoles_r_zBegin();

	double* const soa_quadrupoles_m_r_x = soa.quadrupoles_m_r_xBegin();
	double* const soa_quadrupoles_m_r_y = soa.quadrupoles_m_r_yBegin();
	double* const soa_quadrupoles_m_r_z = soa.quadrupoles_m_r_zBegin();
	double* const soa_quadrupoles_r_x = soa.quadrupoles_r_xBegin();
	double* const soa_quadrupoles_r_y = soa.quadrupoles_r_yBegin();
	double* const soa_quadrupoles_r_z = soa.quadrupoles_r_zBegin();

	size_t iLJCenters = 0;
	size_t iCharges = 0;
	size_t iDipoles = 0;
	size_t iQuadrupoles = 0;
	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < _molecules.size(); ++i) {
		Molecule & M = *_molecules[i];
		const size_t mol_ljc_num = M.numLJcenters();
		const size_t mol_charges_num = M.numCharges();
		const size_t mol_dipoles_num = M.numDipoles();
		const size_t mol_quadrupoles_num = M.numQuadrupoles();
		const double mol_pos_x = M.r(0);
		const double mol_pos_y = M.r(1);
		const double mol_pos_z = M.r(2);

		soa._mol_pos.x(i) = mol_pos_x;
		soa._mol_pos.y(i) = mol_pos_y;
		soa._mol_pos.z(i) = mol_pos_z;

		soa._mol_ljc_num[i] = mol_ljc_num;
		soa._mol_charges_num[i] = mol_charges_num;
		soa._mol_dipoles_num[i] = mol_dipoles_num;
		soa._mol_quadrupoles_num[i] = mol_quadrupoles_num;

		M.setSoA(&soa);
		M.setStartIndexSoA_LJ(iLJCenters);
		M.setStartIndexSoA_C(iCharges);
		M.setStartIndexSoA_D(iDipoles);
		M.setStartIndexSoA_Q(iQuadrupoles);

		M.clearFM();

		M.normalizeQuaternion();

		for (size_t j = 0; j < mol_ljc_num; ++j, ++iLJCenters) {
			double centerPos[3];
			M.computeLJcenter_d(j,centerPos);
			soa_ljc_m_r_x[iLJCenters] = mol_pos_x;
			soa_ljc_m_r_y[iLJCenters] = mol_pos_y;
			soa_ljc_m_r_z[iLJCenters] = mol_pos_z;
			soa_ljc_r_x[iLJCenters] = centerPos[0] + mol_pos_x;
			soa_ljc_r_y[iLJCenters] = centerPos[1] + mol_pos_y;
			soa_ljc_r_z[iLJCenters] = centerPos[2] + mol_pos_z;
			soa._ljc_id[iLJCenters] = compIDs[M.componentid()] + j;
		}
		for (size_t j = 0; j < mol_charges_num; ++j, ++iCharges) {
			double centerPos[3];
			M.computeCharge_d(j, centerPos);
			soa_charges_m_r_x[iCharges] = mol_pos_x;
			soa_charges_m_r_y[iCharges] = mol_pos_y;
			soa_charges_m_r_z[iCharges] = mol_pos_z;
			soa_charges_r_x[iCharges] = centerPos[0] + mol_pos_x;
			soa_charges_r_y[iCharges] = centerPos[1] + mol_pos_y;
			soa_charges_r_z[iCharges] = centerPos[2] + mol_pos_z;
			soa._charges_q[iCharges] = M.component()->charge(j).q();
		}

		for (size_t j = 0; j < mol_dipoles_num; ++j, ++iDipoles) {
			double centerPos[3];
			M.computeDipole_d(j, centerPos);
			double orientation[3];
			M.computeDipole_e(j, orientation);
			soa_dipoles_m_r_x[iDipoles] = mol_pos_x;
			soa_dipoles_m_r_y[iDipoles] = mol_pos_y;
			soa_dipoles_m_r_z[iDipoles] = mol_pos_z;
			soa_dipoles_r_x[iDipoles] = centerPos[0] + mol_pos_x;
			soa_dipoles_r_y[iDipoles] = centerPos[1] + mol_pos_y;
			soa_dipoles_r_z[iDipoles] = centerPos[2] + mol_pos_z;
			soa._dipoles_p[iDipoles] = M.component()->dipole(j).absMy();
			soa._dipoles_e.x(iDipoles) = orientation[0];
			soa._dipoles_e.y(iDipoles) = orientation[1];
			soa._dipoles_e.z(iDipoles) = orientation[2];
		}

		for (size_t j = 0; j < mol_quadrupoles_num; ++j, ++iQuadrupoles) {
			double centerPos[3];
			M.computeQuadrupole_d(j, centerPos);
			double orientation[3];
			M.computeQuadrupole_e(j, orientation);
			soa_quadrupoles_m_r_x[iQuadrupoles] = mol_pos_x;
			soa_quadrupoles_m_r_y[iQuadrupoles] = mol_pos_y;
			soa_quadrupoles_m_r_z[iQuadrupoles] = mol_pos_z;
			soa_quadrupoles_r_x[iQuadrupoles] = centerPos[0] + mol_pos_x;
			soa_quadrupoles_r_y[iQuadrupoles] = centerPos[1] + mol_pos_y;
			soa_quadrupoles_r_z[iQuadrupoles] = centerPos[2] + mol_pos_z;
			soa._quadrupoles_m[iQuadrupoles] = M.component()->quadrupole(j).absQ();
			soa._quadrupoles_e.x(iQuadrupoles) = orientation[0];
			soa._quadrupoles_e.y(iQuadrupoles) = orientation[1];
			soa._quadrupoles_e.z(iQuadrupoles) = orientation[2];
		}
	}
}
