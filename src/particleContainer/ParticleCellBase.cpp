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
#include "utils/generator/EqualVelocityAssigner.h"

CellBorderAndFlagManager ParticleCellBase::_cellBorderAndFlagManager;

ParticleCellBase::ParticleCellBase() {
	// TODO Auto-generated constructor stub
}

ParticleCellBase::~ParticleCellBase() {
	// TODO Auto-generated destructor stub
}

bool ParticleCellBase::deleteMoleculeByID(unsigned long molid) {
	size_t index;
	bool found = findMoleculeByID(index, molid);
	if (found) {
		deleteMoleculeByIndex(index);
	}
	return found;
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

unsigned long ParticleCellBase::initCubicGrid(const std::array<unsigned long, 3> &numMoleculesPerDimension,
											  const std::array<double, 3> &simBoxLength, Random &RNG) {

	std::array<double, 3> spacing{};
	std::array<double, 3> origin1{}; // origin of the first DrawableMolecule
	std::array<double, 3> origin2{}; // origin of the first DrawableMolecule
	for (int d = 0; d < 3; ++d) {
		spacing[d] = simBoxLength[d] / static_cast<double>(numMoleculesPerDimension[d]);
		origin1[d] = spacing[d] * 0.25;
		origin2[d] = spacing[d] * 0.25 * 3.;
	}

	vcp_real_calc T = global_simulation->getEnsemble()->T();
	EqualVelocityAssigner eqVeloAssigner(T);

	double boxMin[3] = {getBoxMin(0), getBoxMin(1), getBoxMin(2)};
	double boxMax[3] = {getBoxMax(0), getBoxMax(1), getBoxMax(2)};

	int start_i = floor((boxMin[0] / simBoxLength[0]) * static_cast<double>(numMoleculesPerDimension[0])) - 1;
	int start_j = floor((boxMin[1] / simBoxLength[1]) * static_cast<double>(numMoleculesPerDimension[1])) - 1;
	int start_k = floor((boxMin[2] / simBoxLength[2]) * static_cast<double>(numMoleculesPerDimension[2])) - 1;

	int end_i = ceil((boxMax[0] / simBoxLength[0]) * static_cast<double>(numMoleculesPerDimension[0])) + 1;
	int end_j = ceil((boxMax[1] / simBoxLength[1]) * static_cast<double>(numMoleculesPerDimension[1])) + 1;
	int end_k = ceil((boxMax[2] / simBoxLength[2]) * static_cast<double>(numMoleculesPerDimension[2])) + 1;

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
		vcp_real_calc x1 = origin1[0] + i_rc * spacing[0];
		vcp_real_calc x2 = origin2[0] + i_rc * spacing[0];
		bool x1In = PositionIsInBox1D(boxMin[0], boxMax[0], static_cast<double>(x1));
		bool x2In = PositionIsInBox1D(boxMin[0], boxMax[0], static_cast<double>(x2));

		if(not (x1In or x2In)) {
			continue;
		}

		for (int j = start_j; j < end_j; ++j) {
			double j_rc = static_cast<double>(j);
			vcp_real_calc y1 = origin1[1] + j_rc * spacing[1];
			vcp_real_calc y2 = origin2[1] + j_rc * spacing[1];
			bool y1In = PositionIsInBox1D(boxMin[1], boxMax[1], static_cast<double>(y1));
			bool y2In = PositionIsInBox1D(boxMin[1], boxMax[1], static_cast<double>(y2));

			if(not (y1In or y2In)) {
				continue;
			}

			for (int k = start_k; k < end_k; ++k) {

				double k_rc = static_cast<double>(k);
				vcp_real_calc z1 = origin1[2] + k_rc * spacing[2];
				vcp_real_calc z2 = origin2[2] + k_rc * spacing[2];
				bool z1In = PositionIsInBox1D(boxMin[2], boxMax[2], static_cast<double>(z1));
				bool z2In = PositionIsInBox1D(boxMin[2], boxMax[2], static_cast<double>(z2));

				if (x1In and y1In and z1In) {
					++numInserted;
					// Init molecule with zero velocity and use the EqualVelocityAssigner in the next step
					Molecule dummy(0, &(global_simulation->getEnsemble()->getComponents()->at(0)),
						x1, y1, z1, 0.0, 0.0, 0.0);
					eqVeloAssigner.assignVelocity(&dummy);
					buffer.push_back(dummy);
				}

				if (x2In and y2In and z2In) {
					++numInserted;
					// Init molecule with zero velocity and use the EqualVelocityAssigner in the next step
					Molecule dummy(0, &(global_simulation->getEnsemble()->getComponents()->at(0)),
						x2, y2, z2, 0.0, 0.0, 0.0);
					
					eqVeloAssigner.assignVelocity(&dummy);
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
