#pragma once

#include <variant>

/**
 * Adds the force values (force (F), torque (M) and virial (Vi)) of haloMolecule to the particle at the given position.
 * The function will reuse the previousIterator if after an increment it points to a particle at the given position.
 * @param moleculeContainer The molecule container.
 * @param position The position of the particle to which the force values will be added.
 * @param previousIterator Can be used to speed up the force exchange. The iterator returned from this function should
 * be used for this.
 * @param haloMolecule The molecule whose force values will be added to the particle at the given position.
 * @return Iterator to the current particle. Can be used in the next iteration as previousIterator.
 */
std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> addValuesAndGetIterator(
	ParticleContainer* moleculeContainer, const double* position,
	std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> previousIterator, Molecule& haloMolecule);
