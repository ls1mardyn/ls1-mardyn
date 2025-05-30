/**
 * @file AutoPasContainer.cpp
 * @author seckler
 * @date 19.09.18
 */
#include "AutoPasContainer.h"
#include <particleContainer/adapter/LegacyCellProcessor.h>
#include <particleContainer/adapter/VectorizedCellProcessor.h>
#include <exception>
#include <sstream>
#include "Domain.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/logging/Logger.h"
#include "parallel/DomainDecompBase.h"
#include "utils/generator/EqualVelocityAssigner.h"

// Declare the main AutoPas class and the computeInteractions() methods with all used functors as extern template
// instantiation. They are instantiated in the respective cpp file inside the templateInstantiations folder.
//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctor<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctor<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctor<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctor<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
#ifdef __AVX__
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorAVX<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorAVX<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorAVX<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorAVX<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
#endif
#ifdef __ARM_FEATURE_SVE
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorSVE<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorSVE<
				Molecule,
				/*applyShift*/ true,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorSVE<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ true,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
extern template bool autopas::AutoPas<Molecule>::computeInteractions(
		mdLib::LJFunctorSVE<
				Molecule,
				/*applyShift*/ false,
				/*mixing*/ false,
				autopas::FunctorN3Modes::Both,
				/*calculateGlobals*/ true
		> *);
#endif
//! @endcond

AutoPasContainer::AutoPasContainer(double cutoff) : _cutoff(cutoff), _particlePropertiesLibrary(cutoff) {
	// use autopas defaults. This block is important when we do not read from an XML like in the unit tests
	_verletSkin = _autopasContainer.getVerletSkin();
	_verletRebuildFrequency = _autopasContainer.getVerletRebuildFrequency();
	_tuningFrequency = _autopasContainer.getTuningInterval();
	_tuningSamples = _autopasContainer.getNumSamples();
	_maxEvidence = _autopasContainer.getMaxEvidence();
	_traversalChoices = _autopasContainer.getAllowedTraversals();
	_containerChoices = _autopasContainer.getAllowedContainers();
	_selectorStrategy = _autopasContainer.getSelectorStrategy();
	_tuningStrategyOptions = _autopasContainer.getTuningStrategyOptions();
	_tuningMetricOption = _autopasContainer.getTuningMetricOption();
	_tuningAcquisitionFunction = _autopasContainer.getAcquisitionFunction();
	_dataLayoutChoices = _autopasContainer.getAllowedDataLayouts();
	_newton3Choices = _autopasContainer.getAllowedNewton3Options();
	_verletClusterSize = _autopasContainer.getVerletClusterSize();
	_relativeOptimumRange = _autopasContainer.getRelativeOptimumRange();
	_relativeBlacklistRange = _autopasContainer.getRelativeBlacklistRange();
	_maxTuningPhasesWithoutTest = _autopasContainer.getMaxTuningPhasesWithoutTest();
	_evidenceForPrediction = _autopasContainer.getEvidenceFirstPrediction();
	_extrapolationMethod = _autopasContainer.getExtrapolationMethodOption();

#ifdef ENABLE_MPI
	std::stringstream logFileName, outputSuffix;

	auto timeNow = std::chrono::system_clock::now();
	auto time_tNow = std::chrono::system_clock::to_time_t(timeNow);

	auto maxRank = global_simulation->domainDecomposition().getNumProcs();
	auto numDigitsMaxRank = std::to_string(maxRank).length();
	auto myRank = global_simulation->domainDecomposition().getRank();

	logFileName << "AutoPas_Rank" << std::setfill('0') << std::setw(numDigitsMaxRank) << myRank << "_"
				<< std::put_time(std::localtime(&time_tNow), "%Y-%m-%d_%H-%M-%S") << ".log";

	_logFile.open(logFileName.str());
	_autopasContainer = decltype(_autopasContainer)(_logFile);

	outputSuffix << "Rank" << std::setfill('0') << std::setw(numDigitsMaxRank) << myRank << "_";
	_autopasContainer.setOutputSuffix(outputSuffix.str());
#endif
}

/**
 * Safe method to parse autopas options from the xml.
 * It checks for exceptions, prints the exceptions and all possible options.
 * If the option does not exist in the xml the given default is returned
 * @tparam OptionType Type of the option to parse
 * @param xmlconfig
 * @param xmlString XML whose content we like to parse.
 * @param defaultValue Set of options that is returned if nothing was found in the xml
 * @return
 */
template <class OptionType, class OutputContainer = std::set<OptionType>>
auto parseAutoPasOption(XMLfileUnits &xmlconfig, const std::string &xmlString,
						const OutputContainer &defaultValue) {
	auto stringInXml = string_utils::toLowercase(xmlconfig.getNodeValue_string(xmlString));
	if (stringInXml.empty()) {
		return defaultValue;
	}
	try {
		return OptionType::template parseOptions<OutputContainer>(stringInXml);
	} catch (const std::exception &e) {
		std::ostringstream error_message;
		error_message << "AutoPasContainer: error when parsing " << xmlString << ":" << std::endl;
		error_message << e.what() << std::endl;
		error_message << "Possible options: "
							<< autopas::utils::ArrayUtils::to_string(OptionType::getAllOptions()) << std::endl;
		MARDYN_EXIT(error_message.str());
		// dummy return
		return decltype(OptionType::template parseOptions<OutputContainer>(""))();
	}
}

void AutoPasContainer::readXML(XMLfileUnits &xmlconfig) {
	std::string oldPath(xmlconfig.getcurrentnodepath());

	// if any option is not specified in the XML use the autopas defaults
	// get option values from xml
	_traversalChoices = parseAutoPasOption<autopas::TraversalOption>(xmlconfig, "allowedTraversals", _traversalChoices);
	_containerChoices = parseAutoPasOption<autopas::ContainerOption>(xmlconfig, "allowedContainers", _containerChoices);
	_selectorStrategy =
		*parseAutoPasOption<autopas::SelectorStrategyOption>(xmlconfig, "selectorStrategy", {_selectorStrategy})
			 .begin();
	_tuningStrategyOptions =
		parseAutoPasOption<autopas::TuningStrategyOption, std::vector<autopas::TuningStrategyOption>>(xmlconfig, "tuningStrategies", {_tuningStrategyOptions});
	_tuningMetricOption = *parseAutoPasOption<autopas::TuningMetricOption>(xmlconfig, "tuningMetric", {_tuningMetricOption}).begin();
	_extrapolationMethod = *parseAutoPasOption<autopas::ExtrapolationMethodOption>(xmlconfig, "extrapolationMethod",
																				   {_extrapolationMethod})
								.begin();
	_dataLayoutChoices = parseAutoPasOption<autopas::DataLayoutOption>(xmlconfig, "dataLayouts", _dataLayoutChoices);
	_newton3Choices = parseAutoPasOption<autopas::Newton3Option>(xmlconfig, "newton3", _newton3Choices);
	_tuningAcquisitionFunction = *parseAutoPasOption<autopas::AcquisitionFunctionOption>(
									  xmlconfig, "tuningAcquisitionFunction", {_tuningAcquisitionFunction})
									  .begin();
	// get numeric options from xml
	// int
	_maxEvidence = static_cast<unsigned int>(xmlconfig.getNodeValue_int("maxEvidence", static_cast<int>(_maxEvidence)));
	_maxTuningPhasesWithoutTest = static_cast<unsigned int>(
		xmlconfig.getNodeValue_int("tuningPhasesWithoutTest", static_cast<int>(_maxTuningPhasesWithoutTest)));
	_evidenceForPrediction = static_cast<unsigned int>(
		xmlconfig.getNodeValue_int("evidenceForPrediction", static_cast<int>(_evidenceForPrediction)));
	_tuningSamples =
		static_cast<unsigned int>(xmlconfig.getNodeValue_int("tuningSamples", static_cast<int>(_tuningSamples)));
	_tuningFrequency =
		static_cast<unsigned int>(xmlconfig.getNodeValue_int("tuningInterval", static_cast<int>(_tuningFrequency)));
	_verletRebuildFrequency = static_cast<unsigned int>(
		xmlconfig.getNodeValue_int("rebuildFrequency", static_cast<int>(_verletRebuildFrequency)));
	_verletClusterSize = static_cast<unsigned int>(
		xmlconfig.getNodeValue_int("verletClusterSize", static_cast<int>(_verletClusterSize)));

	// double
	const auto vlSkin = xmlconfig.getNodeValue_double("skin", -1);
	const auto vlSkinPerTimestep = xmlconfig.getNodeValue_double("skinPerTimestep", -1);
	if (vlSkin == -1 and vlSkinPerTimestep == -1) {
		// stick with the default value
	} else if (vlSkin == -1) {
		_verletSkin = vlSkinPerTimestep * _verletRebuildFrequency;
	} else if (vlSkinPerTimestep == -1 ){
		_verletSkin = vlSkin;
	} else {
		std::ostringstream error_message;
		error_message << "Input XML specifies skin AND skinPerTimestep. Please choose only one." << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	_relativeOptimumRange = xmlconfig.getNodeValue_double("optimumRange", _relativeOptimumRange);
	_relativeBlacklistRange = xmlconfig.getNodeValue_double("blacklistRange", _relativeBlacklistRange);

	// string
	xmlconfig.getNodeValue("ruleFile", _ruleFileName);
	std::string functorChoiceStr{};
	xmlconfig.getNodeValue("functor", functorChoiceStr);
	if (functorChoiceStr.empty()) {
#ifdef __ARM_FEATURE_SVE
		functorChoiceStr = "sve";
#elif defined(__AVX__)
		functorChoiceStr = "avx";
#endif
	}
	if (functorChoiceStr.find("avx") != std::string::npos) {
		functorOption = FunctorOption::AVX;
		Log::global_log->info() << "Selected AVX Functor." << std::endl;
#ifndef __AVX__
		Log::global_log->warning() << "Selected AVX Functor but AVX is not supported! Switching to autoVec functor." << std::endl;
		functorOption = FunctorOption::autoVec;
#endif
	} else if (functorChoiceStr.find("sve") != std::string::npos) {
		functorOption = FunctorOption::SVE;
		Log::global_log->info() << "Selected SVE Functor." << std::endl;
#ifndef __ARM_FEATURE_SVE
		Log::global_log->warning() << "Selected SVE Functor but SVE is not supported! Switching to autoVec functor." << std::endl;
		functorOption = FunctorOption::autoVec;
#endif
	} else {
		functorOption = FunctorOption::autoVec;
		Log::global_log->info() << "Selected autoVec Functor." << std::endl;
	}


	// AutoPas log level
	auto logLevelStr = xmlconfig.getNodeValue_string("logLevel", "");
    // if anything was found try to parse it
	if (not logLevelStr.empty()) {
		// if this is not parsable it defaults to LogLevel::off
		_logLevel = spdlog::level::from_str(logLevelStr);
	}

	xmlconfig.changecurrentnode(oldPath);
}

bool AutoPasContainer::rebuild(double *bBoxMin, double *bBoxMax) {
	mardyn_assert(_cutoff > 0.);
	std::array<double, 3> boxMin{bBoxMin[0], bBoxMin[1], bBoxMin[2]};
	std::array<double, 3> boxMax{bBoxMax[0], bBoxMax[1], bBoxMax[2]};

	memcpy(_boundingBoxMin, bBoxMin, 3 * sizeof(double));
	memcpy(_boundingBoxMax, bBoxMax, 3 * sizeof(double));

	// check if autopas is already initialized
	if (_autopasContainerIsInitialized) {
		const auto emigrants = _autopasContainer.resizeBox(boxMin, boxMax);
		if (not emigrants.empty()) {
			throw std::runtime_error("AutoPasContainer::rebuild(): After resizing the container some particles"
				" were dropped and no effort is made to reinsert them. They should have been collected earlier!" +
				autopas::utils::ArrayUtils::to_string(emigrants, "\n", {"",""}));
		}
		// TODO: maybe only force this if the box and num particles changed too much?
		_autopasContainer.forceRetune();
		return false;
	}

	// The following code is only executed if _autopasContainer has not been initialized, yet.
	_autopasContainer.setBoxMin(boxMin);
	_autopasContainer.setBoxMax(boxMax);
	_autopasContainer.setCutoff(_cutoff);
	_autopasContainer.setVerletSkin(_verletSkin);
	_autopasContainer.setVerletRebuildFrequency(_verletRebuildFrequency);
	_autopasContainer.setVerletClusterSize(_verletClusterSize);
	_autopasContainer.setTuningInterval(_tuningFrequency);
	_autopasContainer.setNumSamples(_tuningSamples);
	_autopasContainer.setSelectorStrategy(_selectorStrategy);
	_autopasContainer.setAllowedContainers(_containerChoices);
	_autopasContainer.setAllowedTraversals(_traversalChoices);
	_autopasContainer.setAllowedDataLayouts(_dataLayoutChoices);
	_autopasContainer.setAllowedNewton3Options(_newton3Choices);
	_autopasContainer.setTuningStrategyOption(_tuningStrategyOptions);
	_autopasContainer.setTuningMetricOption(_tuningMetricOption);
	_autopasContainer.setRuleFileName(_ruleFileName);
	_autopasContainer.setAcquisitionFunction(_tuningAcquisitionFunction);
	_autopasContainer.setMaxEvidence(_maxEvidence);
	_autopasContainer.setRelativeOptimumRange(_relativeOptimumRange);
	_autopasContainer.setMaxTuningPhasesWithoutTest(_maxTuningPhasesWithoutTest);
	_autopasContainer.setRelativeBlacklistRange(_relativeBlacklistRange);
	_autopasContainer.setEvidenceFirstPrediction(_evidenceForPrediction);
	_autopasContainer.setExtrapolationMethodOption(_extrapolationMethod);
	autopas::Logger::get()->set_level(_logLevel);
	_autopasContainer.init();
	_autopasContainerIsInitialized = true;

	// print full configuration to the command line
	constexpr int valueOffset = 28;
	Log::global_log->info() << "AutoPas configuration:\n"
					   << std::setw(valueOffset) << std::left << "Data Layout "
					   << ": " << autopas::utils::ArrayUtils::to_string(_autopasContainer.getAllowedDataLayouts())
					   << "\n"
					   << std::setw(valueOffset) << std::left << "Container "
					   << ": " << autopas::utils::ArrayUtils::to_string(_autopasContainer.getAllowedContainers())
					   << "\n"
					   << std::setw(valueOffset) << std::left << "Cell size Factor "
					   << ": " << _autopasContainer.getAllowedCellSizeFactors() << "\n"
					   << std::setw(valueOffset) << std::left << "Traversals "
					   << ": " << autopas::utils::ArrayUtils::to_string(_autopasContainer.getAllowedTraversals())
					   << "\n"
					   << std::setw(valueOffset) << std::left << "Newton3"
					   << ": " << autopas::utils::ArrayUtils::to_string(_autopasContainer.getAllowedNewton3Options())
					   << "\n"
					   << std::setw(valueOffset) << std::left << "Tuning Metric"
					   << ": " << _autopasContainer.getTuningMetricOption() << "\n"
					   << std::setw(valueOffset) << std::left << "Tuning strategies "
					   << ": " << autopas::utils::ArrayUtils::to_string(_autopasContainer.getTuningStrategyOptions())
					   << "\n"
					   << std::setw(valueOffset) << std::left << "Selector strategy "
					   << ": " << _autopasContainer.getSelectorStrategy() << "\n"
					   << std::setw(valueOffset) << std::left << "Tuning frequency"
					   << ": " << _autopasContainer.getTuningInterval() << "\n"
					   << std::setw(valueOffset) << std::left << "Number of samples "
					   << ": " << _autopasContainer.getNumSamples() << "\n"
					   << std::setw(valueOffset) << std::left << "Tuning Acquisition Function"
					   << ": " << _autopasContainer.getAcquisitionFunction() << "\n"
					   << std::setw(valueOffset) << std::left << "Number of evidence "
					   << ": " << _autopasContainer.getMaxEvidence() << "\n"
					   << std::setw(valueOffset) << std::left << "Verlet Cluster size "
					   << ": " << _autopasContainer.getVerletClusterSize() << "\n"
					   << std::setw(valueOffset) << std::left << "Rebuild frequency "
					   << ": " << _autopasContainer.getVerletRebuildFrequency() << "\n"
					   << std::setw(valueOffset) << std::left << "Verlet Skin "
					   << ": " << _autopasContainer.getVerletSkin() << "\n"
					   << std::setw(valueOffset) << std::left << "Optimum Range "
					   << ": " << _autopasContainer.getRelativeOptimumRange() << "\n"
					   << std::setw(valueOffset) << std::left << "Tuning Phases without test "
					   << ": " << _autopasContainer.getMaxTuningPhasesWithoutTest() << "\n"
					   << std::setw(valueOffset) << std::left << "Blacklist Range "
					   << ": " << _autopasContainer.getRelativeBlacklistRange() << "\n"
					   << std::setw(valueOffset) << std::left << "Evidence for prediction "
					   << ": " << _autopasContainer.getEvidenceFirstPrediction() << "\n"
					   << std::setw(valueOffset) << std::left << "Extrapolation method "
					   << ": " << _autopasContainer.getExtrapolationMethodOption() << "\n"
					   << std::setw(valueOffset) << std::left << "Rule File "
					   << ": " << _autopasContainer.getRuleFileName() << std::endl;

	/// @todo return sendHaloAndLeavingTogether, (always false) for simplicity.
	return false;
}

void AutoPasContainer::update() {
	// in case we update the container before handling the invalid particles, this might lead to lost particles.
	if (not _invalidParticles.empty()) {
		std::ostringstream error_message;
		error_message << "AutoPasContainer: trying to update container, even though invalidParticles still "
							   "exist. This would lead to lost particles => ERROR!\n"
							   "Remaining invalid particles:\n"
							<< autopas::utils::ArrayUtils::to_string(_invalidParticles, "\n", {"", ""})
							<< std::endl;
		MARDYN_EXIT(error_message.str());
	}

	_invalidParticles = _autopasContainer.updateContainer();
}

bool AutoPasContainer::addParticle(Molecule &particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate,
								   const bool &rebuildCaches) {
	if (particle.inBox(_boundingBoxMin, _boundingBoxMax)) {
		_autopasContainer.addParticle(particle);
	} else {
		_autopasContainer.addHaloParticle(particle);
	}
	return true;
}

bool AutoPasContainer::addHaloParticle(Molecule &particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate,
									   const bool &rebuildCaches) {
	_autopasContainer.addHaloParticle(particle);
	return true;
}

void AutoPasContainer::addParticles(std::vector<Molecule> &particles, bool checkWhetherDuplicate) {
	for (auto &particle : particles) {
		addParticle(particle, true, checkWhetherDuplicate);
	}
}

template <typename F>
std::pair<double, double> AutoPasContainer::iterateWithFunctor(F &&functor) {
	// here we call the actual autopas' computeInteractions method to compute the forces.
	_autopasContainer.computeInteractions(&functor);
	const double upot = functor.getPotentialEnergy();
	const double virial = functor.getVirial();
	return std::make_pair(upot, virial);
}

template <bool shifting>
void AutoPasContainer::traverseTemplateHelper() {
#if defined(_OPENMP)
#pragma omp parallel
#endif
	for (auto iter = iterator(ParticleIterator::ALL_CELLS); iter.isValid(); ++iter) {
		iter->clearFM();
	}

	double upot, virial;

	// Check if all components have the same eps24 and sigma. If that is the case, we can skip the mixing rules, which
	// is faster!
	const auto numComponents = _particlePropertiesLibrary.getNumberRegisteredSiteTypes();
	const double epsilonFirstComponent = _particlePropertiesLibrary.getEpsilon(0);
	const double sigmaFirstComponent = _particlePropertiesLibrary.getSigma(0);
	bool allSame = true;
	for (auto i = 1ul; i < numComponents; ++i) {
		allSame &= _particlePropertiesLibrary.getEpsilon(i) == epsilonFirstComponent;
		allSame &= _particlePropertiesLibrary.getSigma(i) == sigmaFirstComponent;
	}

	if (not allSame) {
		Log::global_log->debug() << "AutoPasContainer: Using mixing." << std::endl;
		switch (functorOption) {
			case FunctorOption::SVE: {
#ifdef __ARM_FEATURE_SVE
				// Generate the functor. Should be regenerated every iteration to wipe internally saved globals.
				mdLib::LJFunctorSVE<Molecule, /*applyShift*/ shifting, /*mixing*/ true, autopas::FunctorN3Modes::Both,
						/*calculateGlobals*/ true>
						functor(_cutoff, _particlePropertiesLibrary);

				std::tie(upot, virial) = iterateWithFunctor(functor);
#else
				throw std::runtime_error("SVE Functor not compiled due to lack of compiler support!");
#endif
				break;
			}
			case FunctorOption::AVX: {
#ifdef __AVX__
				// Generate the functor. Should be regenerated every iteration to wipe internally saved globals.
				mdLib::LJFunctorAVX<Molecule, /*applyShift*/ shifting, /*mixing*/ true, autopas::FunctorN3Modes::Both,
						/*calculateGlobals*/ true>
						functor(_cutoff, _particlePropertiesLibrary);

				std::tie(upot, virial) = iterateWithFunctor(functor);
#else
				throw std::runtime_error("AVX Functor not compiled due to lack of compiler support!");
#endif
				break;
			}
			case FunctorOption::autoVec: {
				// Generate the functor. Should be regenerated every iteration to wipe internally saved globals.
				mdLib::LJFunctor<Molecule, /*applyShift*/ shifting, /*mixing*/ true, autopas::FunctorN3Modes::Both,
						/*calculateGlobals*/ true>
						functor(_cutoff, _particlePropertiesLibrary);

				std::tie(upot, virial) = iterateWithFunctor(functor);
				break;
			}
			default: {
				throw std::runtime_error("Unknown functor choice!");
			}
		}
	} else {
		Log::global_log->debug() << "AutoPasContainer: Not using mixing." << std::endl;
		switch (functorOption) {
			case FunctorOption::SVE: {
#ifdef __ARM_FEATURE_SVE
				// Generate the functor. Should be regenerated every iteration to wipe internally saved globals.
				mdLib::LJFunctorSVE<Molecule, /*applyShift*/ shifting, /*mixing*/ false, autopas::FunctorN3Modes::Both,
						/*calculateGlobals*/ true>
						functor(_cutoff);
				functor.setParticleProperties(epsilonFirstComponent * 24, sigmaFirstComponent * sigmaFirstComponent);
				std::tie(upot, virial) = iterateWithFunctor(functor);
#else
				throw std::runtime_error("SVE Functor not compiled due to lack of compiler support!");
#endif
				break;
			}
			case FunctorOption::AVX: {
#ifdef __AVX__
				// Generate the functor. Should be regenerated every iteration to wipe internally saved globals.
				mdLib::LJFunctorAVX<Molecule, /*applyShift*/ shifting, /*mixing*/ false, autopas::FunctorN3Modes::Both,
						/*calculateGlobals*/ true>
						functor(_cutoff);
				functor.setParticleProperties(epsilonFirstComponent * 24, sigmaFirstComponent * sigmaFirstComponent);
				std::tie(upot, virial) = iterateWithFunctor(functor);
#else
				throw std::runtime_error("AVX Functor not compiled due to lack of compiler support!");
#endif
				break;
			}
			case FunctorOption::autoVec: {
				// Generate the functor. Should be regenerated every iteration to wipe internally saved globals.
				mdLib::LJFunctor<Molecule, /*applyShift*/ shifting, /*mixing*/ false, autopas::FunctorN3Modes::Both,
						/*calculateGlobals*/ true>
						functor(_cutoff);
				functor.setParticleProperties(epsilonFirstComponent * 24, sigmaFirstComponent * sigmaFirstComponent);
				std::tie(upot, virial) = iterateWithFunctor(functor);
				break;
			}
			default: {
				throw std::runtime_error("Unknown functor choice!");
			}
		}
	}

	// _myRF is always zero for lj only!
	global_simulation->getDomain()->setLocalVirial(virial /*+ 3.0 * _myRF*/);
	// _upotXpoles is zero as we do not have any dipoles or quadrupoles
	global_simulation->getDomain()->setLocalUpot(upot /* _upotXpoles + _myRF*/);
}

void AutoPasContainer::traverseCells(CellProcessor &cellProcessor) {
	if (dynamic_cast<VectorizedCellProcessor *>(&cellProcessor) or
		dynamic_cast<LegacyCellProcessor *>(&cellProcessor)) {
		// only initialize ppl if it is empty
		bool hasShift = false;
		bool hasNoShift = false;

		if (_particlePropertiesLibrary.getNumberRegisteredSiteTypes() == 0) {
			const auto components = global_simulation->getEnsemble()->getComponents();
			for (const auto &c : *components) {
				_particlePropertiesLibrary.addSiteType(c.getLookUpId(),c.ljcenter(0).m());
				_particlePropertiesLibrary.addLJParametersToSite(c.getLookUpId(), c.ljcenter(0).eps(), c.ljcenter(0).sigma());
			}
			_particlePropertiesLibrary.calculateMixingCoefficients();
			size_t numComponentsAdded = 0;
			for (const auto &c : *components) {
				if (c.ljcenter(0).shift6() != 0.) {
					hasShift = true;
					const double autoPasShift6 =
						_particlePropertiesLibrary.getMixingShift6(numComponentsAdded, numComponentsAdded);
					const double ls1Shift6 = c.ljcenter(0).shift6();
					if (std::fabs((autoPasShift6 - ls1Shift6) / ls1Shift6) > 1.e-10) {
						// warn if shift differs relatively by more than 1.e-10
						Log::global_log->warning() << "Dangerous shift6 detected: AutoPas will use: " << autoPasShift6
											  << ", while normal ls1 mode uses: " << ls1Shift6 << std::endl
											  << "Please check that your shifts are calculated correctly." << std::endl;
					}
					++numComponentsAdded;
				} else {
					hasNoShift = true;
				}
			}
		}
		if (hasShift and hasNoShift) {
			// if some particles require shifting and some don't:
			// throw an error, as AutoPas does not support this, yet.
			throw std::runtime_error("AutoPas does not support mixed shifting state!");
		}
		if (hasShift) {
			traverseTemplateHelper</*shifting*/ true>();
		} else {
			traverseTemplateHelper</*shifting*/ false>();
		}

	} else {
		Log::global_log->warning() << "only lj functors are supported for traversals." << std::endl;
	}
}

void AutoPasContainer::traverseNonInnermostCells(CellProcessor &cellProcessor) {
	throw std::runtime_error("AutoPasContainer::traverseNonInnermostCells() not yet implemented");
}

void AutoPasContainer::traversePartialInnermostCells(CellProcessor &cellProcessor, unsigned int stage, int stageCount) {
	throw std::runtime_error("AutoPasContainer::traversePartialInnermostCells() not yet implemented");
}

unsigned long AutoPasContainer::getNumberOfParticles(ParticleIterator::Type t /* = ParticleIterator::ONLY_INNER_AND_BOUNDARY */) {
	unsigned long count = 0;
	for (auto iter = iterator(t); iter.isValid(); ++iter) {
		++count;
	}
	return count;
	// return _autopasContainer.getNumberOfParticles(); // todo: this is currently buggy!, so we use iterators instead.
}

void AutoPasContainer::clear() { _autopasContainer.deleteAllParticles(); }

void AutoPasContainer::deleteOuterParticles() {
	Log::global_log->info() << "deleting outer particles by using forced update" << std::endl;
	auto invalidParticles = _autopasContainer.updateContainer();
	if (not invalidParticles.empty()) {
		throw std::runtime_error(
			"AutoPasContainer: deleteOuterParticles(): Invalid particles are returned! Please ensure that you are "
			"properly updating your container!");
	}
}

double AutoPasContainer::getHaloWidthForDimension(int /*index*/) const { return _cutoff; }

double AutoPasContainer::getCutoff() const { return _cutoff; }

double AutoPasContainer::getSkin() const { return _verletSkin; }

size_t AutoPasContainer::getRebuildFrequency() const {return _verletRebuildFrequency;};

void AutoPasContainer::deleteMolecule(ParticleIterator &moleculeIter, const bool & /*rebuildCaches*/) {
	_autopasContainer.deleteParticle(moleculeIter);
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

std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> AutoPasContainer::getMoleculeAtPosition(
	const double *pos) {
	std::array<double, 3> pos_arr{pos[0], pos[1], pos[2]};
	for (auto iter = this->iterator(ParticleIterator::ALL_CELLS); iter.isValid(); ++iter) {
		if (iter->getR() == pos_arr) {
			return iter;
		}
	}
	return {};  // default initialized iter is invalid.
}

unsigned long AutoPasContainer::initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
											  std::array<double, 3> simBoxLength) {

	vcp_real_calc T = global_simulation->getEnsemble()->T();
	EqualVelocityAssigner eqVeloAssigner(T);

	const std::array<double, 3> spacing = autopas::utils::ArrayMath::div(simBoxLength,
																		 autopas::utils::ArrayUtils::static_cast_copy_array<double>(
																				 numMoleculesPerDimension));
	size_t numMolecules = 0;
	for (int grid = 0; grid < 2; ++grid) {
		// grid starts 1/4 away from the corner. Grids are offset by 1/2 spacing
		const std::array<double, 3> offset = autopas::utils::ArrayMath::mulScalar(spacing, 0.25 + .5 * grid);
		for (int z = 0; z < numMoleculesPerDimension[2]; ++z) {
			const double posZ = offset[2] + spacing[2] * z;
			for (int y = 0; y < numMoleculesPerDimension[1]; ++y) {
				const double posY = offset[1] + spacing[1] * y;
				for (int x = 0; x < numMoleculesPerDimension[0]; ++x) {
					const double posX = offset[0] + spacing[0] * x;
					// Init molecule with zero velocity and use the EqualVelocityAssigner in the next step
					Molecule m(numMolecules++, &(global_simulation->getEnsemble()->getComponents()->at(0)), posX, posY,
							   posZ, 0.0, 0.0, 0.0);
					eqVeloAssigner.assignVelocity(&m);
					addParticle(m, true, false, false);
				}
			}
		}
	}
	return numMolecules;
}

double *AutoPasContainer::getCellLength() {
	throw std::runtime_error("AutoPasContainer::getCellLength() not yet implemented");
}

double *AutoPasContainer::getHaloSize() {
	static std::array<double, 3> haloLength{_verletSkin + _cutoff, _verletSkin + _cutoff, _verletSkin + _cutoff};
	return haloLength.data();
}

autopas::IteratorBehavior convertBehaviorToAutoPas(ParticleIterator::Type t) {
	switch (t) {
		case ParticleIterator::Type::ALL_CELLS:
			return autopas::IteratorBehavior::ownedOrHalo;
		case ParticleIterator::Type::ONLY_INNER_AND_BOUNDARY:
			return autopas::IteratorBehavior::owned;
	}
	throw std::runtime_error("Unknown iterator type.");
}

ParticleIterator AutoPasContainer::iterator(ParticleIterator::Type t) {
	return ParticleIterator{_autopasContainer.begin(convertBehaviorToAutoPas(t))};
}

RegionParticleIterator AutoPasContainer::regionIterator(const double *startCorner, const double *endCorner,
														ParticleIterator::Type t) {
	const std::array<double, 3> lowCorner{startCorner[0], startCorner[1], startCorner[2]};
	const std::array<double, 3> highCorner{endCorner[0], endCorner[1], endCorner[2]};
	return RegionParticleIterator{
		_autopasContainer.getRegionIterator(lowCorner, highCorner, convertBehaviorToAutoPas(t))};
}
std::string AutoPasContainer::getConfigurationAsString() { return _autopasContainer.getCurrentConfigs().at(autopas::InteractionTypeOption::pairwise).get().toString(); }
