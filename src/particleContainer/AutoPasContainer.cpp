/**
 * @file AutoPasContainer.cpp
 * @author seckler
 * @date 19.09.18
 */
#include <exception>
#include <particleContainer/adapter/VectorizedCellProcessor.h>
#include <particleContainer/adapter/LegacyCellProcessor.h>
#include "AutoPasContainer.h"
#include "Simulation.h"
#include "Domain.h"

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
	memcpy(_boundingBoxMin, bBoxMin, 3*sizeof(double));
	memcpy(_boundingBoxMax, bBoxMax, 3*sizeof(double));
	/// @todo return sendHaloAndLeavingTogether, (always false) for simplicity.
	return false;
}

void AutoPasContainer::update() {
	_autopasContainer.updateContainer();
}

bool AutoPasContainer::addParticle(Molecule &particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate,
                                   const bool &rebuildCaches) {
	if(particle.inBox(_boundingBoxMin, _boundingBoxMax)){
		_autopasContainer.addParticle(particle);
	} else {
		_autopasContainer.addHaloParticle(particle);
	}

}

bool AutoPasContainer::addHaloParticle(Molecule &particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate,
                                       const bool &rebuildCaches) {
	_autopasContainer.addHaloParticle(particle);
	return true;
}

void AutoPasContainer::addParticles(std::vector<Molecule> &particles, bool checkWhetherDuplicate) {
	for (auto &particle: particles) {
		addParticle(particle, true, checkWhetherDuplicate);
	}
}

void AutoPasContainer::traverseCells(CellProcessor &cellProcessor) {
	if(dynamic_cast<VectorizedCellProcessor*> (&cellProcessor) or dynamic_cast<LegacyCellProcessor*> (&cellProcessor)){

		double epsilon, sigma, shift;
		{
			auto iter = iterator();
			if (not iter.isValid()) {
				return;
			}
			auto ljcenter = iter->component()->ljcenter(0);
			epsilon = ljcenter.eps();
			sigma = ljcenter.sigma();
			shift = ljcenter.shift6()/6.;

		}
		// lower and upper corner of the local domain needed to correctly calculate the global values
		std::array<double, 3> lowCorner = {_boundingBoxMin[0],_boundingBoxMin[1],_boundingBoxMin[2]};
		std::array<double, 3> highCorner = {_boundingBoxMax[0],_boundingBoxMax[1],_boundingBoxMax[2]};

		// generate the functor
		autopas::LJFunctor<Molecule, CellType, /*calculateGlobals*/ true> functor(_cutoff, epsilon, sigma, shift,
																				  lowCorner, highCorner,
																				  /*duplicatedCalculation*/ true);

		functor.resetGlobalValues();
		_autopasContainer.iteratePairwise(&functor, autopas::DataLayoutOption::aos);
		functor.postProcessGlobalValues(/*newton3*/ true);
		double upot = functor.getUpot();
		double virial = functor.getVirial();

		// _myRF is always zero for lj only!
		global_simulation->getDomain()->setLocalVirial(virial /*+ 3.0 * _myRF*/);
		// _upotXpoles is zero as we do not have any dipoles or quadrupoles
		global_simulation->getDomain()->setLocalUpot(upot /* _upotXpoles + _myRF*/);
	} else {
		global_log->warning() << "only lj functors are supported for traversals." << std::endl;
	}
}

void AutoPasContainer::traverseNonInnermostCells(CellProcessor &cellProcessor) {
	throw std::runtime_error("not yet implemented");
}

void AutoPasContainer::traversePartialInnermostCells(CellProcessor &cellProcessor, unsigned int stage, int stageCount) {
	throw std::runtime_error("not yet implemented");
}

unsigned long AutoPasContainer::getNumberOfParticles() {
	return _autopasContainer.getNumberOfParticles();
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
	std::array<double, 3> pos_arr{pos[0], pos[1], pos[2]};
	for (auto iter = iterator(); iter.isValid(); ++iter) {
		if (iter->getR() == pos_arr) {
			*result = &(*iter);
			return true;
		}
	}
	return false;
}

unsigned long AutoPasContainer::initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
                                              std::array<double, 3> simBoxLength) {
	throw std::runtime_error("not yet implemented");
}

double *AutoPasContainer::getCellLength() {
	throw std::runtime_error("not yet implemented");
}

autopas::IteratorBehavior convertBehaviorToAutoPas(ParticleIterator::Type t){
	switch(t){
		case ParticleIterator::Type::ALL_CELLS:
			return autopas::IteratorBehavior::haloAndOwned;
		case ParticleIterator::Type::ONLY_INNER_AND_BOUNDARY:
			return autopas::IteratorBehavior::ownedOnly;
	}
}

ParticleIterator AutoPasContainer::iterator(ParticleIterator::Type t) {

	return _autopasContainer.begin(convertBehaviorToAutoPas(t));
}

RegionParticleIterator
AutoPasContainer::regionIterator(const double *startCorner, const double *endCorner, ParticleIterator::Type t) {
	std::array<double, 3> lowCorner{startCorner[0],startCorner[1],startCorner[2]};
	std::array<double, 3> highCorner{endCorner[0],endCorner[1],endCorner[2]};
	return _autopasContainer.getRegionIterator(lowCorner, highCorner, convertBehaviorToAutoPas(t));
}



