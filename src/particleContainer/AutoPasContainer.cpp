/**
 * @file AutoPasContainer.cpp
 * @author seckler
 * @date 19.09.18
 */
#include <exception>
#include "AutoPasContainer.h"

AutoPasContainer::AutoPasContainer() :
		_cutoff(0.), _verletSkin(0.3), _verletRebuildFrequency(10u), _tuningFrequency(100u), _autopasContainer() {
}

void AutoPasContainer::readXML(XMLfileUnits &xmlconfig) {}

bool AutoPasContainer::rebuild(double *bBoxMin, double *bBoxMax) {
	mardyn_assert(_cutoff > 0.);
	std::array<double, 3> boxMin{bBoxMin[0], bBoxMin[1], bBoxMin[2]};
	std::array<double, 3> boxMax{bBoxMax[0], bBoxMax[1], bBoxMax[2]};
	_autopasContainer.init(boxMin, boxMax, _cutoff, _verletSkin, _verletRebuildFrequency, autopas::allContainerOptions,
	                       autopas::allTraversalOptions, _tuningFrequency);

	/// @todo return sendHaloAndLeavingTogether, (always false) for simplicity.
	return false;
}

void AutoPasContainer::update() {
	_autopasContainer.updateContainer();
}

bool AutoPasContainer::addParticle(Molecule &particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate,
                                   const bool &rebuildCaches) {
	_autopasContainer.addParticle(particle);
}

void AutoPasContainer::addParticles(std::vector<Molecule> &particles, bool checkWhetherDuplicate) {
	for (auto &particle: particles) {
		addParticle(particle, true, checkWhetherDuplicate);
	}
}

void AutoPasContainer::traverseCells(CellProcessor &cellProcessor) {
	autopas::LJFunctor<Molecule, CellType> functor;
	functor.setGlobals(_cutoff, 1., 1., 0.);
	_autopasContainer.iteratePairwise(&functor, autopas::DataLayoutOption::aos);
}

void AutoPasContainer::traverseNonInnermostCells(CellProcessor &cellProcessor) {
	throw std::runtime_error("not yet implemented");
}

void AutoPasContainer::traversePartialInnermostCells(CellProcessor &cellProcessor, unsigned int stage, int stageCount) {
	throw std::runtime_error("not yet implemented");
}

unsigned long AutoPasContainer::getNumberOfParticles() {
	return 0;
}

void AutoPasContainer::clear() {
	_autopasContainer.deleteAllParticles();
}

void AutoPasContainer::deleteOuterParticles() {
	_autopasContainer.deleteHaloParticles();
}

double AutoPasContainer::get_halo_L(int /*index*/) const {
	return _cutoff;
}

double AutoPasContainer::getCutoff() {
	return _cutoff;
}

void AutoPasContainer::deleteMolecule(Molecule &molecule, const bool &rebuildCaches) {
	throw std::runtime_error("not yet implemented");
}

double AutoPasContainer::getEnergy(ParticlePairsHandler *particlePairsHandler, Molecule *m1,
								   CellProcessor &cellProcessor) {
	throw std::runtime_error("not yet implemented");
}

void AutoPasContainer::updateInnerMoleculeCaches() {
	throw std::runtime_error("not yet implemented");
}

void AutoPasContainer::updateBoundaryAndHaloMoleculeCaches() {
	throw std::runtime_error("not yet implemented");
}

void AutoPasContainer::updateMoleculeCaches() {
	// nothing needed
}

bool AutoPasContainer::getMoleculeAtPosition(const double *pos, Molecule **result) {
	return false;
}

unsigned long AutoPasContainer::initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
                                              std::array<double, 3> simBoxLength) {
	throw std::runtime_error("not yet implemented");
}

double *AutoPasContainer::getCellLength() {
	throw std::runtime_error("not yet implemented");
}

