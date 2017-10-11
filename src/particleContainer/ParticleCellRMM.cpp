#include "particleContainer/ParticleCellRMM.h"
#include "particleContainer/ParticleCell.h"

CellBorderAndFlagManager ParticleCellRMM::_cellBorderAndFlagManager;

ParticleCellRMM::ParticleCellRMM() : _cellDataSoARMM(0) {
	// TODO Auto-generated constructor stub

}

ParticleCellRMM::~ParticleCellRMM() {
	// TODO Auto-generated destructor stub
}

void ParticleCellRMM::deallocateAllParticles() {
	_cellDataSoARMM.resize(0);
}

bool ParticleCellRMM::addParticle(Molecule& particle, bool checkWhetherDuplicate) {

#ifndef NDEBUG
	double boxMin[3] = {getBoxMin(0), getBoxMin(1), getBoxMin(2)};
	double boxMax[3] = {getBoxMax(0), getBoxMax(1), getBoxMax(2)};
	mardyn_assert(particle.inBox(boxMin, boxMax));
#endif

	bool wasInserted = false;
	bool found = false;

	if (checkWhetherDuplicate) {
		size_t index;
		found = findMoleculeByID(index, particle.id());
	}

	if (not found) {
		_cellDataSoARMM.appendMolecule(particle);
		wasInserted = true;
	} else {
		wasInserted = false;
	}
	return wasInserted;
}

Molecule ParticleCellRMM::buildAoSMolecule(size_t i) const {
	return _cellDataSoARMM.buildAoSMolecule(i);
}

void ParticleCellRMM::moleculesAtNew(size_t i, Molecule*& multipurposePointer) {
	mardyn_assert((int)i < getMoleculeCount());
	_cellDataSoARMM.readMutableMolecule(i, *multipurposePointer);
}

void ParticleCellRMM::moleculesAtConstNew(size_t i, Molecule*& multipurposePointer) const {
	mardyn_assert((int)i < getMoleculeCount());
	_cellDataSoARMM.readImmutableMolecule(i, *multipurposePointer);
}


bool ParticleCellRMM::isEmpty() const {
	return getMoleculeCount() == 0;
}

bool ParticleCellRMM::deleteMoleculeByIndex(size_t index) {
	bool found = true;
	mardyn_assert(index < static_cast<size_t>(getMoleculeCount()));
	_cellDataSoARMM.deleteMolecule(index);
	return found;
}

int ParticleCellRMM::getMoleculeCount() const {
	return _cellDataSoARMM.getMolNum();
}

void ParticleCellRMM::updateLeavingMoleculesBase(ParticleCellBase& otherCell) {
	ParticleCellRMM& oCell = downcastCellReferenceRMM(otherCell);
	//const int oNumMols = oCell.getMoleculeCount();

	// TODO: use getRegion and put in a vector? will reduce number of calls to inBox by factor 2?
	double boxMin[3] = {getBoxMin(0), getBoxMin(1), getBoxMin(2)};
	double boxMax[3] = {getBoxMax(0), getBoxMax(1), getBoxMax(2)};
	double otherCellBoxMin[3] = {otherCell.getBoxMin(0), otherCell.getBoxMin(1), otherCell.getBoxMin(2)};
	double otherCellBoxMax[3] = {otherCell.getBoxMax(0), otherCell.getBoxMax(1), otherCell.getBoxMax(2)};

	// how many molecules travel from my cell to the other one
	const int numMolsMeToOther = countInRegion(otherCellBoxMin, otherCellBoxMax);

	// how many molecules travel from the other cell to me
	int numMolsOtherToMe = oCell.countInRegion(boxMin, boxMax);

	if (numMolsMeToOther == 0 and numMolsOtherToMe == 0)
		return;

	if (numMolsMeToOther >= numMolsOtherToMe) {
		swapAndAppendToCell(oCell);
	} else {
		oCell.swapAndAppendToCell(*this);
	}
}

void ParticleCellRMM::swapAndAppendToCell(ParticleCellRMM& other) {
	double boxMin[3] = {getBoxMin(0), getBoxMin(1), getBoxMin(2)};
	double boxMax[3] = {getBoxMax(0), getBoxMax(1), getBoxMax(2)};
	double otherCellBoxMin[3] = {other.getBoxMin(0), other.getBoxMin(1), other.getBoxMin(2)};
	double otherCellBoxMax[3] = {other.getBoxMax(0), other.getBoxMax(1), other.getBoxMax(2)};

	// holds: the number of molecules that move from "this" to "other" is >= 0
	Molecule dummy, otherDummy;
	Molecule *dummy_p = &dummy, *otherDummy_p = &otherDummy;
	int j = 0;
	for (int i = 0; i < getMoleculeCount(); ++i) {
		// find next molecule to move from this to other
		moleculesAtNew(i, dummy_p);
		if (not dummy.inBox(otherCellBoxMin, otherCellBoxMax))
			continue;

		// find next swap/insertion position
		for (;j <other.getMoleculeCount(); ++j) {
			other.moleculesAtNew(j, otherDummy_p);
			if (otherDummy.inBox(boxMin, boxMax))
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

void ParticleCellRMM::swapMolecules(int i, ParticleCellRMM& other, int j) {
	Molecule read_i, read_j;
	      _cellDataSoARMM.readImmutableMolecule(i, read_i);
	other._cellDataSoARMM.readImmutableMolecule(j, read_j);

	      _cellDataSoARMM.writeMolecule(i, read_j);
	other._cellDataSoARMM.writeMolecule(j, read_i);


}

int ParticleCellRMM::countInRegion(double lowCorner[3], double highCorner[3]) const {
	int ret = 0;
	const int totalNumMols = getMoleculeCount();

	const vcp_real_calc * const soa1_mol_pos_x = _cellDataSoARMM.r_xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = _cellDataSoARMM.r_yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = _cellDataSoARMM.r_zBegin();


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

void ParticleCellRMM::getRegion(double lowCorner[3], double highCorner[3],
		std::vector<Molecule*>& particlePtrs, bool removeFromContainer) {
}

void ParticleCellRMM::increaseMoleculeStorage(size_t numExtraMols) {
	_cellDataSoARMM.increaseStorage(numExtraMols);
}

void ParticleCellRMM::prefetchForForce() const {
	_cellDataSoARMM.prefetchForForce();
}

void ParticleCellRMM::getLeavingMolecules(std::vector<Molecule>& appendBuffer) {
	const vcp_real_calc * const soa1_mol_pos_x = _cellDataSoARMM.r_xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = _cellDataSoARMM.r_yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = _cellDataSoARMM.r_zBegin();
	double boxMin[3] = {getBoxMin(0), getBoxMin(1), getBoxMin(2)};
	double boxMax[3] = {getBoxMax(0), getBoxMax(1), getBoxMax(2)};

	for (int i = 0; i < getMoleculeCount(); ) {

		double x = static_cast<double>(soa1_mol_pos_x[i]);
		double y = static_cast<double>(soa1_mol_pos_y[i]);
		double z = static_cast<double>(soa1_mol_pos_z[i]);
		bool xOut = x < boxMin[0] or x >= boxMax[0];
		bool yOut = y < boxMin[1] or y >= boxMax[1];
		bool zOut = z < boxMin[2] or z >= boxMax[2];
		bool satisfy = xOut or yOut or zOut;

		if (not satisfy) {
			// do nothing, advance i
			++i;
		} else {
			appendBuffer.push_back(buildAoSMolecule(i));
			deleteMoleculeByIndex(i);
			// don't advance i
		}
	}
}
