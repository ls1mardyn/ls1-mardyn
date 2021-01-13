#pragma once

#include <variant>

/**
 *
 * @param moleculeContainer
 * @param position
 * @param previousIterator
 * @param haloMolecule
 * @return
 */
std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> addValuesAndGetIterator(
	ParticleContainer* moleculeContainer, const double* position,
	std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> previousIterator,
	Molecule& haloMolecule);