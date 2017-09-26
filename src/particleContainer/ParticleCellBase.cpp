/*
 * ParticleCellBase.cpp
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#include "ParticleCellBase.h"
#include "ensemble/EnsembleBase.h"
#include "utils/Random.h"
#include "Simulation.h"

ParticleCellBase::ParticleCellBase() {
	// TODO Auto-generated constructor stub

}

ParticleCellBase::~ParticleCellBase() {
	// TODO Auto-generated destructor stub
}

bool ParticleCellBase::deleteMoleculeByID(unsigned long molid) {
	bool found = false;

	size_t index;
	findMoleculeByID(found, index, molid);
	if (found) {
		deleteMoleculeByIndex(index);
	}

	return found;
}

void ParticleCellBase::findMoleculeByID(bool& wasFound, size_t& index, unsigned long molid) const {
	wasFound = false;
	int numMolecules = getMoleculeCount();

	for (int i = 0; i < numMolecules; ++i) {
		if (moleculesAtConst(i).id() == molid) {
			index = i;
			wasFound = true;
			break;
		}
	}
}

template <typename T>
std::array<T, 3> getRandomVelocity(T temperature, Random& RNG) {
	std::array<T,3> ret;

	// Velocity
	for (int dim = 0; dim < 3; dim++) {
		ret[dim] = RNG.uniformRandInRange(-0.5f, 0.5f);
	}
	T dotprod_v = 0;
	for (unsigned int i = 0; i < ret.size(); i++) {
		dotprod_v += ret[i] * ret[i];
	}
	// Velocity Correction
	const T three = static_cast<T>(3.0);
	T vCorr = sqrt(three * temperature / dotprod_v);
	for (unsigned int i = 0; i < ret.size(); i++) {
		ret[i] *= vCorr;
	}

	return ret;
}

template <typename T>
bool PositionIsInBox1D(const T l, const T u, const T r) {
#ifdef __INTEL_COMPILER
		#pragma float_control(precise, on)
		#pragma fenv_access(on)
#endif
	return r >= l and r < u;
}


template <typename T>
bool PositionIsInBox3D(const T l[3], const T u[3], const T r[3]) {
	bool in = true;
	for (int d = 0; d < 3; ++d) {
		in &= PositionIsInBox1D(l[d], u[d], r[d]);
	}
	return in;
}

unsigned long ParticleCellBase::initCubicGrid(int numMoleculesPerDimension,
		double simBoxLength, Random & RNG) {

	double spacing = simBoxLength / numMoleculesPerDimension;
	double origin1 = spacing * 0.25; // origin of the first DrawableMolecule
	double origin2 = spacing * 0.25 * 3.; // origin of the first DrawableMolecule

	vcp_real_calc T = global_simulation->getEnsemble()->T();

	int start_i = floor((getBoxMin(0) / simBoxLength) * numMoleculesPerDimension) - 1;
	int start_j = floor((getBoxMin(1) / simBoxLength) * numMoleculesPerDimension) - 1;
	int start_k = floor((getBoxMin(2) / simBoxLength) * numMoleculesPerDimension) - 1;

	int end_i = ceil((getBoxMax(0) / simBoxLength) * numMoleculesPerDimension) + 1;
	int end_j = ceil((getBoxMax(1) / simBoxLength) * numMoleculesPerDimension) + 1;
	int end_k = ceil((getBoxMax(2) / simBoxLength) * numMoleculesPerDimension) + 1;

	std::vector<Molecule> buffer;
	int numMolsUpperBound = (end_k - start_k) * (end_j - start_j) * (end_i - start_i) * 2;
	buffer.reserve(numMolsUpperBound);

	// step 1: count!
	unsigned long numInserted = 0;

	const bool isHalo = isHaloCell();
	const bool notHalo = not isHalo;


	for (int i = start_i; i < end_i; ++i) {

		//// positions should be initially created in vcp_real_calc, but comparisons to _boxMin and _boxMax should be in double precision!

		double i_rc = static_cast<double>(i);
		vcp_real_calc x1 = origin1 + i_rc * spacing;
		vcp_real_calc x2 = origin2 + i_rc * spacing;
		bool x1In = PositionIsInBox1D(_boxMin[0], _boxMax[0], static_cast<double>(x1));
		bool x2In = PositionIsInBox1D(_boxMin[0], _boxMax[0], static_cast<double>(x2));

		if(not (x1In or x2In)) {
			continue;
		}

		for (int j = start_j; j < end_j; ++j) {
			double j_rc = static_cast<double>(j);
			vcp_real_calc y1 = origin1 + j_rc * spacing;
			vcp_real_calc y2 = origin2 + j_rc * spacing;
			bool y1In = PositionIsInBox1D(_boxMin[1], _boxMax[1], static_cast<double>(y1));
			bool y2In = PositionIsInBox1D(_boxMin[1], _boxMax[1], static_cast<double>(y2));

			if(not (y1In or y2In)) {
				continue;
			}

			for (int k = start_k; k < end_k; ++k) {

				double k_rc = static_cast<double>(k);
				vcp_real_calc z1 = origin1 + k_rc * spacing;
				vcp_real_calc z2 = origin2 + k_rc * spacing;
				bool z1In = PositionIsInBox1D(_boxMin[2], _boxMax[2], static_cast<double>(z1));
				bool z2In = PositionIsInBox1D(_boxMin[2], _boxMax[2], static_cast<double>(z2));

				if (x1In and y1In and z1In) {
					++numInserted;
					std::array<vcp_real_calc,3> v = getRandomVelocity<vcp_real_calc>(T, RNG);
					Molecule dummy(0, &(global_simulation->getEnsemble()->getComponents()->at(0)),
						x1, y1, z1, v[0], -v[1], v[2]);
					buffer.push_back(dummy);
				}

				if (x2In and y2In and z2In) {
					++numInserted;
					std::array<vcp_real_calc,3> v = getRandomVelocity<vcp_real_calc>(T, RNG);
					Molecule dummy(0, &(global_simulation->getEnsemble()->getComponents()->at(0)),
						x2, y2, z2, v[0], -v[1], v[2]);
					buffer.push_back(dummy);
				}
			}
		}
	}

	// step 2: preallocate!
	increaseMoleculeStorage(numInserted);

	// step 3: add. Also for Halo-cells, we want the current thread to also "touch" the memory, not just allocate it.
	// lookup OpenMP "first touch" policy
	for (unsigned long i = 0; i < numInserted; ++i) {
		addParticle(buffer[i]);
	}

	if (isHalo) {
		deallocateAllParticles();

		// and return 0, instead of a number.
		numInserted = 0;
	}

	return numInserted;
}
