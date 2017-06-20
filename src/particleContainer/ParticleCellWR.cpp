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

Molecule& ParticleCell_WR::moleculesAt(size_t i) {
	mardyn_assert((int)i < getMoleculeCount());
	_cellDataSoA_WR.readMutableMolecule(i, _dummy);
	return _dummy;
}

const Molecule& ParticleCell_WR::moleculesAtConst(size_t i) const {
	mardyn_assert((int)i < getMoleculeCount());
	_cellDataSoA_WR.readImmutableMolecule(i, const_cast<Molecule &>(_dummy));
	return _dummy;
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
	return _cellDataSoA_WR._mol_num;
}

void ParticleCell_WR::updateLeavingMoleculesBase(ParticleCellBase& otherCell) {
	ParticleCell_WR& oCell = downcastCellReferenceWR(otherCell);
	const int oNumMols = oCell.getMoleculeCount();

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
	int j = 0;
	for (int i = 0; i < getMoleculeCount(); ++i) {
		// find next molecule to move from this to other
		if (not moleculesAt(i).inBox(other._boxMin, other._boxMax))
			continue;

		// find next swap/insertion position
		for (;j <other.getMoleculeCount(); ++j) {
			if (not other.moleculesAt(j).inBox(_boxMin,_boxMax))
				continue;
		}

		if (j < other.getMoleculeCount()) {
			swapMolecules(i, other, j);
		} else {
			other.addParticle(moleculesAt(i));
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

	const vcp_real_calc * const soa1_mol_pos_x = _cellDataSoA_WR._mol_r.xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = _cellDataSoA_WR._mol_r.yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = _cellDataSoA_WR._mol_r.zBegin();
	vcp_real_calc low[3], high[3];
	for(int d = 0; d < 3; ++d) {
		low[d] = static_cast<vcp_real_calc>(lowCorner[d]);
		high[d] = static_cast<vcp_real_calc>(highCorner[d]);
	}

	#if defined(_OPENMP)
	#pragma omp simd reduction(+:ret) aligned(soa1_mol_pos_x, soa1_mol_pos_y, soa1_mol_pos_z: 64)
	#endif
	for (int i = 0; i < totalNumMols; ++i) {
		bool isIn = soa1_mol_pos_x[i] >= low[0] and soa1_mol_pos_x[i] < high[0] and
					soa1_mol_pos_y[i] >= low[1] and soa1_mol_pos_y[i] < high[1] and
					soa1_mol_pos_z[i] >= low[2] and soa1_mol_pos_z[i] < high[2];
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
