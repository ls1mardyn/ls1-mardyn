#include <particleContainer/ParticleContainer.h>
#include <particleContainer/ParticleIterator.h>
#include <variant>

std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> addValuesAndGetIterator(
	ParticleContainer* moleculeContainer, const double* position,
	std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> previousIterator,
	Molecule& haloMolecule) {
	// Reuse the previous iterator in case the particles match.
	bool found = false;
	std::visit(
		[&](auto& previousIter) {
			if (previousIter.isValid()) {
				++previousIter;
				if (previousIter.isValid()) {
					const double epsi = moleculeContainer->getCutoff() * 1e-6;
					if (fabs(previousIter->r(0) - position[0]) <= epsi and
						fabs(previousIter->r(1) - position[1]) <= epsi and
						fabs(previousIter->r(2) - position[2]) <= epsi) {
						found = true;
					}
				}
			}
		},
		previousIterator);
	// If the previousIterator matches, use it!
	auto originalIter = found ? previousIterator : moleculeContainer->getMoleculeAtPosition(position);

	std::visit(
		[&](auto originalIter) {
			if (not originalIter.isValid()) {
				// This should not happen
				std::cout << "Original molecule not found";
				mardyn_exit(1);
			}

			mardyn_assert(originalIter->getID() == haloMolecule.getID());

			originalIter->Fadd(haloMolecule.F_arr().data());
			originalIter->Madd(haloMolecule.M_arr().data());
			originalIter->Viadd(haloMolecule.Vi_arr().data());
		},
		originalIter);
	return originalIter;
}