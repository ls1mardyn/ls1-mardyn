/**
 * @file AutoPasContainer.cpp
 * @author seckler
 * @date 19.09.18
 */
#include "AutoPasContainer.h"
#include <particleContainer/adapter/LegacyCellProcessor.h>
#include <particleContainer/adapter/VectorizedCellProcessor.h>
#include <exception>
#include "Domain.h"
#include "Simulation.h"
#include "autopas/utils/Logger.h"
#include "autopas/utils/StringUtils.h"
#include "parallel/DomainDecompBase.h"

AutoPasContainer::AutoPasContainer()
	: _cutoff(0.),
	  _verletSkin(0.3),
	  _verletRebuildFrequency(10u),
	  _tuningFrequency(1000u),
	  _tuningSamples(3u),
	  _autopasContainer(),
	  _traversalChoices(autopas::allTraversalOptions),
	  _containerChoices(autopas::allContainerOptions),
	  _selectorStrategy(autopas::SelectorStrategyOption::fastestMedian),
	  _dataLayoutChoices{autopas::DataLayoutOption::soa},
	  _newton3Choices{autopas::Newton3Option::enabled} {
#ifdef ENABLE_MPI
	std::stringstream logFileName;

	auto timeNow = chrono::system_clock::now();
	auto time_tNow = std::chrono::system_clock::to_time_t(timeNow);

	auto maxRank = global_simulation->domainDecomposition().getNumProcs();
	auto numDigitsMaxRank = std::to_string(maxRank).length();

	logFileName << "AutoPas_Rank" << setfill('0') << setw(numDigitsMaxRank)
				<< global_simulation->domainDecomposition().getRank() << "_"
				<< std::put_time(std::localtime(&time_tNow), "%Y-%m-%d_%H-%M-%S") << ".log";

	_logFile.open(logFileName.str());
	_autopasContainer = decltype(_autopasContainer)(_logFile);
#endif
}

void AutoPasContainer::readXML(XMLfileUnits &xmlconfig) {
	string oldPath(xmlconfig.getcurrentnodepath());

	// set default values here!

	_traversalChoices = autopas::utils::StringUtils::parseTraversalOptions(
		string_utils::toLowercase(xmlconfig.getNodeValue_string("allowedTraversals", "c08")));
	_containerChoices = autopas::utils::StringUtils::parseContainerOptions(
		string_utils::toLowercase(xmlconfig.getNodeValue_string("allowedContainers", "linked-cell")));

	_selectorStrategy = autopas::utils::StringUtils::parseSelectorStrategy(
		string_utils::toLowercase(xmlconfig.getNodeValue_string("selectorStrategy", "median")));

	_dataLayoutChoices = autopas::utils::StringUtils::parseDataLayout(
		string_utils::toLowercase(xmlconfig.getNodeValue_string("dataLayouts", "soa")));

	_newton3Choices = autopas::utils::StringUtils::parseNewton3Options(
		string_utils::toLowercase(xmlconfig.getNodeValue_string("newton3", "enabled")));

	_tuningSamples = (unsigned int)xmlconfig.getNodeValue_int("tuningSamples", 3);
	_tuningFrequency = (unsigned int)xmlconfig.getNodeValue_int("tuningInterval", 500);

	xmlconfig.getNodeValue("rebuildFrequency", _verletRebuildFrequency);
	xmlconfig.getNodeValue("skin", _verletSkin);

	std::stringstream dataLayoutChoicesStream;
	for_each(_dataLayoutChoices.begin(), _dataLayoutChoices.end(),
			 [&](auto &choice) { dataLayoutChoicesStream << autopas::utils::StringUtils::to_string(choice) << " "; });
	std::stringstream containerChoicesStream;
	for_each(_containerChoices.begin(), _containerChoices.end(),
			 [&](auto &choice) { containerChoicesStream << autopas::utils::StringUtils::to_string(choice) << " "; });
	std::stringstream traversalChoicesStream;
	for_each(_traversalChoices.begin(), _traversalChoices.end(),
			 [&](auto &choice) { traversalChoicesStream << autopas::utils::StringUtils::to_string(choice) << " "; });
	std::stringstream newton3ChoicesStream;
	for_each(_newton3Choices.begin(), _newton3Choices.end(),
			 [&](auto &choice) { newton3ChoicesStream << autopas::utils::StringUtils::to_string(choice) << " "; });

	int valueOffset = 20;
	global_log->info() << "AutoPas configuration:" << endl
					   << setw(valueOffset) << left << "Data Layout "
					   << ": " << dataLayoutChoicesStream.str() << endl
					   << setw(valueOffset) << left << "Container "
					   << ": " << containerChoicesStream.str() << endl
					   << setw(valueOffset) << left << "Traversals "
					   << ": " << traversalChoicesStream.str() << endl
					   << setw(valueOffset) << left << "Selector strategy "
					   << ": " << autopas::utils::StringUtils::to_string(_selectorStrategy) << endl
					   << setw(valueOffset) << left << "Tuning frequency"
					   << ": " << _tuningFrequency << endl
					   << setw(valueOffset) << left << "Number of samples "
					   << ": " << _tuningSamples << endl
					   << setw(valueOffset) << left << "Newton3"
					   << ": " << newton3ChoicesStream.str() << endl;
	xmlconfig.changecurrentnode(oldPath);
}

bool AutoPasContainer::rebuild(double *bBoxMin, double *bBoxMax) {
	mardyn_assert(_cutoff > 0.);
	std::array<double, 3> boxMin{bBoxMin[0], bBoxMin[1], bBoxMin[2]};
	std::array<double, 3> boxMax{bBoxMax[0], bBoxMax[1], bBoxMax[2]};

	_autopasContainer.setBoxMin(boxMin);
	_autopasContainer.setBoxMax(boxMax);
	_autopasContainer.setCutoff(_cutoff);
	_autopasContainer.setVerletSkin(_verletSkin);
	_autopasContainer.setVerletRebuildFrequency(_verletRebuildFrequency);
	_autopasContainer.setTuningInterval(_tuningFrequency);
	_autopasContainer.setNumSamples(_tuningSamples);
	_autopasContainer.setSelectorStrategy(_selectorStrategy);
	_autopasContainer.setAllowedContainers(_containerChoices);
	_autopasContainer.setAllowedTraversals(_traversalChoices);
	_autopasContainer.setAllowedDataLayouts(_dataLayoutChoices);
	_autopasContainer.setAllowedNewton3Options(_newton3Choices);
	_autopasContainer.init();
	autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);

	memcpy(_boundingBoxMin, bBoxMin, 3 * sizeof(double));
	memcpy(_boundingBoxMax, bBoxMax, 3 * sizeof(double));
	/// @todo return sendHaloAndLeavingTogether, (always false) for simplicity.
	return false;
}

void AutoPasContainer::update() {
	if (not _invalidParticles.empty()) {
		global_log->error() << "AutoPasContainer: trying to update container, even though invalidParticles still "
							   "exist. This would lead to lost particle => ERROR!"
							<< std::endl;
		Simulation::exit(434);
	}

	std::tie(_invalidParticles, _hasInvalidParticles) = _autopasContainer.updateContainer();
}

void AutoPasContainer::forcedUpdate() {
	if (not _invalidParticles.empty()) {
		global_log->error() << "AutoPasContainer: trying to force update container, even though invalidParticles still "
							   "exist. This would lead to lost particle => ERROR!"
							<< std::endl;
		Simulation::exit(435);
	}
	_hasInvalidParticles = true;
	_invalidParticles = _autopasContainer.updateContainerForced();
}

bool AutoPasContainer::addParticle(Molecule &particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate,
								   const bool &rebuildCaches) {
	if (particle.inBox(_boundingBoxMin, _boundingBoxMax)) {
		_autopasContainer.addParticle(particle);
	} else {
		_autopasContainer.addOrUpdateHaloParticle(particle);
	}
	return true;
}

bool AutoPasContainer::addHaloParticle(Molecule &particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate,
									   const bool &rebuildCaches) {
	_autopasContainer.addOrUpdateHaloParticle(particle);
	return true;
}

void AutoPasContainer::addParticles(std::vector<Molecule> &particles, bool checkWhetherDuplicate) {
	for (auto &particle : particles) {
		addParticle(particle, true, checkWhetherDuplicate);
	}
}

void AutoPasContainer::traverseCells(CellProcessor &cellProcessor) {
	if (dynamic_cast<VectorizedCellProcessor *>(&cellProcessor) or
		dynamic_cast<LegacyCellProcessor *>(&cellProcessor)) {
		global_log->info() << "AutoPasContainer: traverseCells" << std::endl;
		double epsilon=0., sigma=0., shift=0.;
		{
			auto iter = iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			if (iter.isValid()) {
				auto ljcenter = iter->component()->ljcenter(0);
				epsilon = ljcenter.eps();
				sigma = ljcenter.sigma();
				shift = ljcenter.shift6() / 6.;
			}
		}
		// lower and upper corner of the local domain needed to correctly calculate the global values
		std::array<double, 3> lowCorner = {_boundingBoxMin[0], _boundingBoxMin[1], _boundingBoxMin[2]};
		std::array<double, 3> highCorner = {_boundingBoxMax[0], _boundingBoxMax[1], _boundingBoxMax[2]};

		// generate the functor
		autopas::LJFunctor<Molecule, CellType, autopas::FunctorN3Modes::Both,
						   /*calculateGlobals*/ true>
			functor(_cutoff, epsilon, sigma, shift, /*duplicatedCalculation*/ true);
#if defined(_OPENMP)
#pragma omp parallel
#endif
		for (auto iter = iterator(ParticleIterator::ALL_CELLS); iter.isValid(); ++iter) {
			iter->clearFM();
		}

		_autopasContainer.iteratePairwise(&functor);
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
	throw std::runtime_error("AutoPasContainer::traverseNonInnermostCells() not yet implemented");
}

void AutoPasContainer::traversePartialInnermostCells(CellProcessor &cellProcessor, unsigned int stage, int stageCount) {
	throw std::runtime_error("AutoPasContainer::traversePartialInnermostCells() not yet implemented");
}

unsigned long AutoPasContainer::getNumberOfParticles() {
	unsigned long count = 0;
	for(auto iter = iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter){
		++count;
	}
	return count;
	//return _autopasContainer.getNumberOfParticles(); // this is currently buggy!, so we use iterators instead.
}

void AutoPasContainer::clear() { _autopasContainer.deleteAllParticles(); }

void AutoPasContainer::deleteOuterParticles() {
	global_log->info() << "deleting outer particles by using forced update" << std::endl;
	auto invalidParticles = _autopasContainer.updateContainerForced();
	if (not invalidParticles.empty()) {
		throw std::runtime_error(
			"AutoPasContainer: Invalid particles ignored in deleteOuterParticles, check that your rebalance rate is a "
			"multiple of the rebuild rate!");
	}
}

double AutoPasContainer::get_halo_L(int /*index*/) const { return _cutoff; }

double AutoPasContainer::getCutoff() const { return _cutoff; }

double AutoPasContainer::getInteractionLength() const { return _cutoff + _verletSkin; }

double AutoPasContainer::getSkin() const { return _verletSkin; }

void AutoPasContainer::deleteMolecule(Molecule &molecule, const bool &rebuildCaches) {
	throw std::runtime_error("AutoPasContainer::deleteMolecule() not yet implemented");
}

double AutoPasContainer::getEnergy(ParticlePairsHandler *particlePairsHandler, Molecule *m1,
								   CellProcessor &cellProcessor) {
	throw std::runtime_error("AutoPasContainer::getEnergy() not yet implemented");
}

void AutoPasContainer::updateInnerMoleculeCaches() {
	throw std::runtime_error("AutoPasContainer::updateInnerMoleculeCaches() not yet implemented");
}

void AutoPasContainer::updateBoundaryAndHaloMoleculeCaches() {
	throw std::runtime_error("AutoPasContainer::updateBoundaryAndHaloMoleculeCaches() not yet implemented");
}

void AutoPasContainer::updateMoleculeCaches() {
	// nothing needed
}

bool AutoPasContainer::getMoleculeAtPosition(const double *pos, Molecule **result) {
	std::array<double, 3> pos_arr{pos[0], pos[1], pos[2]};
	for (auto iter = this->iterator(ParticleIterator::ALL_CELLS); iter.isValid(); ++iter) {
		if (iter->getR() == pos_arr) {
			*result = &(*iter);
			return true;
		}
	}
	return false;
}

unsigned long AutoPasContainer::initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
											  std::array<double, 3> simBoxLength) {
	throw std::runtime_error("AutoPasContainer::initCubicGrid() not yet implemented");
}

double *AutoPasContainer::getCellLength() {
	throw std::runtime_error("AutoPasContainer::getCellLength() not yet implemented");
}

double *AutoPasContainer::getHaloSize() {
	static std::array<double, 3> haloLength{_verletSkin + _cutoff};
	return haloLength.data();
}

autopas::IteratorBehavior convertBehaviorToAutoPas(ParticleIterator::Type t) {
	switch (t) {
		case ParticleIterator::Type::ALL_CELLS:
			return autopas::IteratorBehavior::haloAndOwned;
		case ParticleIterator::Type::ONLY_INNER_AND_BOUNDARY:
			return autopas::IteratorBehavior::ownedOnly;
	}
	throw std::runtime_error("Unknown iterator type.");
}

ParticleIterator AutoPasContainer::iterator(ParticleIterator::Type t) {
	return _autopasContainer.begin(convertBehaviorToAutoPas(t));
}

RegionParticleIterator AutoPasContainer::regionIterator(const double *startCorner, const double *endCorner,
														ParticleIterator::Type t) {
	std::array<double, 3> lowCorner{startCorner[0], startCorner[1], startCorner[2]};
	std::array<double, 3> highCorner{endCorner[0], endCorner[1], endCorner[2]};
	return RegionParticleIterator{
		_autopasContainer.getRegionIterator(lowCorner, highCorner, convertBehaviorToAutoPas(t))};
}
