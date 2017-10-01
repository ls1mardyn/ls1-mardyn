/*
 * ParticleCellWR.cpp
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#include "particleContainer/ParticleCellWR.h"
#include "particleContainer/ParticleCell.h"

ParticleCell_WR::ParticleCell_WR() : _cellDataSoA_WR(0) {
	// TODO Auto-generated constructor stub

}

ParticleCell_WR::~ParticleCell_WR() {
	// TODO Auto-generated destructor stub
}

void ParticleCell_WR::deallocateAllParticles() {
	_cellDataSoA_WR.resize(0);
}

bool ParticleCell_WR::addParticle(Molecule& particle, bool checkWhetherDuplicate) {

#ifndef NDEBUG
	mardyn_assert(particle.inBox(_boxMin,_boxMax ));
#endif

	bool wasInserted;
	bool found = false;

	if (checkWhetherDuplicate) {
		size_t index;
		findMoleculeByID(found, index, particle.id());
	}

	if (not found) {
		_cellDataSoA_WR.appendMolecule(particle);
		wasInserted = true;
	} else {
		wasInserted = false;
	}

	return wasInserted;
}

void ParticleCell_WR::moleculesAtNew(size_t i, Molecule*& multipurposePointer) {
	mardyn_assert((int)i < getMoleculeCount());
	_cellDataSoA_WR.readMutableMolecule(i, *multipurposePointer);
}

void ParticleCell_WR::moleculesAtConstNew(size_t i, Molecule*& multipurposePointer) const {
	mardyn_assert((int)i < getMoleculeCount());
	_cellDataSoA_WR.readImmutableMolecule(i, *multipurposePointer);
}


bool ParticleCell_WR::isEmpty() const {
	return getMoleculeCount() == 0;
}

bool ParticleCell_WR::deleteMoleculeByIndex(size_t index) {
	bool found = true;
	mardyn_assert(index < static_cast<size_t>(getMoleculeCount()));
	_cellDataSoA_WR.deleteMolecule(index);
	return found;
}

int ParticleCell_WR::getMoleculeCount() const {
	return _cellDataSoA_WR.getMolNum();
}

void ParticleCell_WR::updateLeavingMoleculesBase(ParticleCellBase& otherCell) {
	ParticleCell_WR& oCell = downcastCellReferenceWR(otherCell);
	//const int oNumMols = oCell.getMoleculeCount();

	// TODO: use getRegion and put in a vector? will reduce number of calls to inBox by factor 2?

	// how many molecules travel from my cell to the other one
	const int numMolsMeToOther = countInRegion(oCell._boxMin, oCell._boxMax);

	// how many molecules travel from the other cell to me
	int numMolsOtherToMe = oCell.countInRegion(_boxMin, _boxMax);

	if (numMolsMeToOther == 0 and numMolsOtherToMe == 0)
		return;

	if (numMolsMeToOther >= numMolsOtherToMe) {
		swapAndAppendToCell(oCell);
	} else {
		oCell.swapAndAppendToCell(*this);
	}
}

void ParticleCell_WR::swapAndAppendToCell(ParticleCell_WR& other) {
	// holds: the number of molecules that move from "this" to "other" is >= 0
	Molecule dummy, otherDummy;
	Molecule *dummy_p = &dummy, *otherDummy_p = &otherDummy;
	int j = 0;
	for (int i = 0; i < getMoleculeCount(); ++i) {
		// find next molecule to move from this to other
		moleculesAtNew(i, dummy_p);
		if (not dummy.inBox(other._boxMin, other._boxMax))
			continue;

		// find next swap/insertion position
		for (;j <other.getMoleculeCount(); ++j) {
			other.moleculesAtNew(j, otherDummy_p);
			if (otherDummy.inBox(_boxMin,_boxMax))
				break;
		}

		if (j < other.getMoleculeCount()) {
			swapMolecules(i, other, j);
		} else {
			other.addParticle(dummy);
			deleteMoleculeByIndex(i);
			--i;
		}
	}
}

void ParticleCell_WR::swapMolecules(int i, ParticleCell_WR& other, int j) {
	Molecule read_i, read_j;
	      _cellDataSoA_WR.readImmutableMolecule(i, read_i);
	other._cellDataSoA_WR.readImmutableMolecule(j, read_j);

	      _cellDataSoA_WR.writeMolecule(i, read_j);
	other._cellDataSoA_WR.writeMolecule(j, read_i);


}

int ParticleCell_WR::countInRegion(double lowCorner[3], double highCorner[3]) const {
	int ret = 0;
	const int totalNumMols = getMoleculeCount();

	const vcp_real_calc * const soa1_mol_pos_x = _cellDataSoA_WR.r_xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = _cellDataSoA_WR.r_yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = _cellDataSoA_WR.r_zBegin();


	#if defined(_OPENMP)
	#pragma omp simd reduction(+:ret) aligned(soa1_mol_pos_x, soa1_mol_pos_y, soa1_mol_pos_z: 64)
	#endif
	for (int i = 0; i < totalNumMols; ++i) {
		bool isIn = soa1_mol_pos_x[i] >= lowCorner[0] and soa1_mol_pos_x[i] < highCorner[0] and
					soa1_mol_pos_y[i] >= lowCorner[1] and soa1_mol_pos_y[i] < highCorner[1] and
					soa1_mol_pos_z[i] >= lowCorner[2] and soa1_mol_pos_z[i] < highCorner[2];
		ret += static_cast<int>(isIn);
	}

	return ret;
}

void ParticleCell_WR::getRegion(double lowCorner[3], double highCorner[3],
		std::vector<Molecule*>& particlePtrs, bool removeFromContainer) {
}

void ParticleCell_WR::increaseMoleculeStorage(size_t numExtraMols) {
	_cellDataSoA_WR.increaseStorage(numExtraMols);
}

void ParticleCell_WR::prefetch() const {
	_cellDataSoA_WR.prefetch();
}
