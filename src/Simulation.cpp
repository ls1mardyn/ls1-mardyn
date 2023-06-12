#define SIMULATION_SRC
#include "Simulation.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "Common.h"
#include "Domain.h"
#include "particleContainer/LinkedCells.h"

#ifdef MARDYN_AUTOPAS
#include "particleContainer/AutoPasContainer.h"
#endif

#include "parallel/DomainDecompBase.h"
#include "parallel/NonBlockingMPIHandlerBase.h"
#include "parallel/NonBlockingMPIMultiStepHandler.h"
#include "molecules/Molecule.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#include "parallel/GeneralDomainDecomposition.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/adapter/VCP1CLJRMM.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "integrators/LeapfrogRMM.h"

#include "plugins/PluginBase.h"
#include "plugins/PluginFactory.h"

#include "io/MmpldWriter.h"
#include "io/RDF.h"
#include "io/FlopRateWriter.h"

#include "io/ASCIIReader.h"
#include "io/BinaryReader.h"
#include "io/CubicGridGeneratorInternal.h"
#include "io/MemoryProfiler.h"
#include "io/Mkesfera.h"
#include "io/MultiObjectGenerator.h"
#include "io/PerCellGenerator.h"
#include "io/ReplicaGenerator.h"
#include "io/TcTS.h"
#include "io/TimerProfiler.h"
#ifdef ENABLE_ADIOS2
#include "io/Adios2Reader.h"
#endif

#include "ensemble/GrandCanonicalEnsemble.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/CavityEnsemble.h"

#include "thermostats/VelocityScalingThermostat.h"
#include "thermostats/TemperatureControl.h"

#include "utils/FileUtils.h"
#include "utils/Logger.h"

#include "longRange/LongRangeCorrection.h"
#include "longRange/Homogeneous.h"
#include "longRange/Planar.h"
#include "longRange/NoLRC.h"

#include "bhfmm/FastMultipoleMethod.h"
#include "bhfmm/cellProcessors/VectorizedLJP2PCellProcessor.h"


Simulation* global_simulation;

Simulation::Simulation()
	:
	_simulationTime(0.0),
	_maxMoleculeId(0),
	_cutoffRadius(0.0),
	_LJCutoffRadius(0.0),
	_collectThermostatDirectedVelocity(100),
	// ANDERSEN DEPRECATED
	_thermostatType(VELSCALE_THERMOSTAT),
	_numberOfTimesteps(1),
	_simstep(0),
	_initSimulation(0),
	_initCanonical(_initSimulation + 1),  // simulate canonical (with Tfactor == 1.) from begin on!
	_initGrandCanonical(10000000),
	_initStatistics(20000),
	_ensemble(nullptr),
	_moleculeContainer(nullptr),
	_particlePairsHandler(nullptr),
	_cellProcessor(nullptr),
	_domainDecomposition(nullptr),
	_integrator(nullptr),
	_domain(nullptr),
	_inputReader(nullptr),
	_outputPrefix("mardyn"),
	_momentumInterval(1000),
	_rand(8624),
	_longRangeCorrection(nullptr),
	_temperatureControl(nullptr),
	_FMM(nullptr),
	_timerProfiler(),
#ifdef TASKTIMINGPROFILE
	_taskTimingProfiler(new TaskTimingProfiler),
#endif /* TASKTIMINGPROFILE */
	_forced_checkpoint_time(0),
	_loopCompTime(0.0),
	_loopCompTimeSteps(0)
{
	_timeFromStart.start();
	_ensemble = new CanonicalEnsemble();
	initialize();
}

Simulation::~Simulation() {
	delete _ensemble;
	_ensemble = nullptr;
	delete _moleculeContainer;
	_moleculeContainer = nullptr;
	delete _particlePairsHandler;
	_particlePairsHandler = nullptr;
	delete _cellProcessor;
	_cellProcessor = nullptr;
	delete _domainDecomposition;
	_domainDecomposition = nullptr;
	delete _integrator;
	_integrator = nullptr;
	delete _domain;
	_domain = nullptr;
	delete _inputReader;
	_inputReader = nullptr;
	delete _longRangeCorrection;
	_longRangeCorrection = nullptr;
	delete _temperatureControl;
	_temperatureControl = nullptr;
	delete _FMM;
	_FMM = nullptr;

	/* destruct plugins and remove from plugin list */
	_plugins.remove_if([](PluginBase *pluginPtr) {delete pluginPtr; return true;} );
}

void Simulation::exit(int exitcode) {
	// .. to avoid code duplication ..
	mardyn_exit(exitcode);
}

void Simulation::readXML(XMLfileUnits& xmlconfig) {
	/* timers */
	if(xmlconfig.changecurrentnode("programtimers")) {
		_timerProfiler.readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}

	/* integrator */
	if(xmlconfig.changecurrentnode("integrator")) {
		std::string integratorType;
		xmlconfig.getNodeValue("@type", integratorType);
		Log::global_log->info() << "Integrator type: " << integratorType << std::endl;
		if(integratorType == "Leapfrog") {
#ifdef ENABLE_REDUCED_MEMORY_MODE
			Log::global_log->error() << "The reduced memory mode (RMM) requires the LeapfrogRMM integrator." << std::endl;
			Simulation::exit(-1);
#endif
			_integrator = new Leapfrog();
		} else if (integratorType == "LeapfrogRMM") {
			_integrator = new LeapfrogRMM();
		} else {
			Log::global_log-> error() << "Unknown integrator " << integratorType << std::endl;
			Simulation::exit(1);
		}
		_integrator->readXML(xmlconfig);
		_integrator->init();
		xmlconfig.changecurrentnode("..");
	} else {
		Log::global_log->error() << "Integrator section missing." << std::endl;
	}

	/* run section */
	if(xmlconfig.changecurrentnode("run")) {
		xmlconfig.getNodeValueReduced("currenttime", _simulationTime);
		Log::global_log->info() << "Simulation start time: " << _simulationTime << std::endl;
		/* steps */
		xmlconfig.getNodeValue("equilibration/steps", _initStatistics);
		Log::global_log->info() << "Number of equilibration steps: " << _initStatistics << std::endl;
		xmlconfig.getNodeValue("production/steps", _numberOfTimesteps);
		Log::global_log->info() << "Number of timesteps: " << _numberOfTimesteps << std::endl;
		xmlconfig.getNodeValue("production/loop-abort-time", _maxWallTime);
		if(_maxWallTime != -1) {
			Log::global_log->info() << "Max loop-abort-time set: " << _maxWallTime << std::endl;
			_wallTimeEnabled = true;
		}
		xmlconfig.changecurrentnode("..");
	} else {
		Log::global_log->error() << "Run section missing." << std::endl;
	}

	/* ensemble */
	if(xmlconfig.changecurrentnode("ensemble")) {
		std::string ensembletype;
		xmlconfig.getNodeValue("@type", ensembletype);
		Log::global_log->info() << "Ensemble: " << ensembletype<< std::endl;
		if (ensembletype == "NVT") {
			delete _ensemble;
			_ensemble = new CanonicalEnsemble();
		} else if (ensembletype == "muVT") {
			Log::global_log->error() << "muVT ensemble not completely implemented via XML input." << std::endl;
			/* TODO: GrandCanonical terminates as readXML is not implemented and it is not tested yet. */
			_ensemble = new GrandCanonicalEnsemble();
		} else {
			Log::global_log->error() << "Unknown ensemble type: " << ensembletype << std::endl;
			Simulation::exit(1);
		}
		_ensemble->readXML(xmlconfig);
		/** @todo Here we store data in the _domain member as long as we do not use the ensemble everywhere */
		for (int d = 0; d < 3; d++) {
			_domain->setGlobalLength(d, _ensemble->domain()->length(d));
		}
		_domain->setGlobalTemperature(_ensemble->T());

		xmlconfig.changecurrentnode("..");
	}
	else {
		Log::global_log->error() << "Ensemble section missing." << std::endl;
		Simulation::exit(1);
	}

	//The mixing coefficents have to be read in the ensemble part
	//if this didn't happen fill the mixing coeffs with default values
	auto& dmixcoeff = _domain->getmixcoeff();
	const std::size_t compNum = _ensemble->getComponents()->size();
	//1 Comps: 0 coeffs; 2 Comps: 2 coeffs; 3 Comps: 6 coeffs; 4 Comps 12 coeffs
	const std::size_t neededCoeffs = compNum*(compNum-1);
	if(dmixcoeff.size() < neededCoeffs){
		Log::global_log->warning() << "Not enough mixing coefficients were given! (Filling the missing ones with 1)" << '\n';
		Log::global_log->warning() << "This can happen because the xml-input doesn't support these yet!" << std::endl;
		unsigned int numcomponents = _simulation.getEnsemble()->getComponents()->size();
		for (unsigned int i = dmixcoeff.size(); i < neededCoeffs; i++) {
			dmixcoeff.push_back(1);
		}
	}

	/* algorithm */
	if(xmlconfig.changecurrentnode("algorithm")) {
		/* cutoffs */
		if(xmlconfig.changecurrentnode("cutoffs")) {
			if(xmlconfig.getNodeValueReduced("defaultCutoff", _cutoffRadius)) {
				Log::global_log->info() << "dimensionless default cutoff radius:\t" << _cutoffRadius << std::endl;
			}
			if(xmlconfig.getNodeValueReduced("radiusLJ", _LJCutoffRadius)) {
				Log::global_log->info() << "dimensionless LJ cutoff radius:\t" << _LJCutoffRadius << std::endl;
				for(auto &component: *(_ensemble->getComponents())) {
					component.updateAllLJcentersShift(_LJCutoffRadius);
				}
			}
			/** @todo introduce maxCutoffRadius here for datastructures, ...
			 *        maybe use map/list to store cutoffs for different potentials? */
			_cutoffRadius = std::max(_cutoffRadius, _LJCutoffRadius);
			if(_cutoffRadius <= 0) {
				Log::global_log->error() << "cutoff radius <= 0." << std::endl;
				Simulation::exit(1);
			}
			Log::global_log->info() << "dimensionless cutoff radius:\t" << _cutoffRadius << std::endl;
			xmlconfig.changecurrentnode("..");
		} else {
			Log::global_log->error() << "Cutoff section missing." << std::endl;
			Simulation::exit(1);
		}

		/* electrostatics */
		/** @todo This may be better go into a physical section for constants? */
		if(xmlconfig.changecurrentnode("electrostatic[@type='ReactionField']")) {
			double epsilonRF = 0;
			xmlconfig.getNodeValueReduced("epsilon", epsilonRF);
			Log::global_log->info() << "Epsilon Reaction Field: " << epsilonRF << std::endl;
			_domain->setepsilonRF(epsilonRF);
			xmlconfig.changecurrentnode("..");
		} else {
			Log::global_log->error() << "Electrostatics section for reaction field setup missing." << std::endl;
			Simulation::exit(1);
		}

		if (xmlconfig.changecurrentnode("electrostatic[@type='FastMultipoleMethod']")) {
#ifdef MARDYN_AUTOPAS
			Log::global_log->fatal()
				<< "The fast multipole method is not compatible with AutoPas. Please disable the AutoPas mode (ENABLE_AUTOPAS)!"
				<< std::endl;
			Simulation::exit(1);
#endif
			_FMM = new bhfmm::FastMultipoleMethod();
			_FMM->readXML(xmlconfig);
			xmlconfig.changecurrentnode("..");
		}

		/* parallelisation */
		if(xmlconfig.changecurrentnode("parallelisation")) {
			std::string parallelisationtype("DomainDecomposition");
			xmlconfig.getNodeValue("@type", parallelisationtype);
			Log::global_log->info() << "Parallelisation type: " << parallelisationtype << std::endl;
			xmlconfig.getNodeValue("overlappingP2P", _overlappingP2P);
			Log::global_log->info() << "Using overlapping p2p communication: " << _overlappingP2P << std::endl;
		#ifdef ENABLE_MPI
			/// @todo Dummy Decomposition now included in DecompBase - still keep this name?
			if(parallelisationtype == "DummyDecomposition") {
				Log::global_log->error() << "DummyDecomposition not available in parallel mode." << std::endl;
			}
			else if(parallelisationtype == "DomainDecomposition") {
				delete _domainDecomposition;
				_domainDecomposition = new DomainDecomposition();
			}
			else if(parallelisationtype == "KDDecomposition") {
				delete _domainDecomposition;
				_domainDecomposition = new KDDecomposition(getcutoffRadius(), _ensemble->getComponents()->size());
			} else if (parallelisationtype == "GeneralDomainDecomposition") {
				double skin = 0.;
				bool forceLatchingToLinkedCellsGrid = false;
				// We need the skin here (to specify the smallest possible partition), so we extract it from the AutoPas
				// container's xml, because the ParticleContainer needs to be instantiated later. :/
				xmlconfig.changecurrentnode("..");
				if (xmlconfig.changecurrentnode("datastructure")) {
					std::string datastructuretype;
					xmlconfig.getNodeValue("@type", datastructuretype);
					if (datastructuretype == "AutoPas" or datastructuretype == "AutoPasContainer") {
						// check if skin is specified
						if (xmlconfig.getNodeValue("skin", skin) == 0) {
							Log::global_log->error() << "Skin not set in datastructure/AutoPas. "
												   "This will lead to a different interaction length in the container "
												   "vs the GeneralDomainDecomposition which can lead ALL to shrink the "
												   "domain too small."
												  << std::endl;
							this->exit(512435340);
						}
					} else {
						Log::global_log->warning() << "Using the GeneralDomainDecomposition without AutoPas is not "
												 "thoroughly tested and considered BETA."
											  << std::endl;
						// Force grid! This is needed, as the linked cells container assumes a grid and the calculation
						// of global values will be faulty without one!
						Log::global_log->info() << "Forcing a grid for the GeneralDomainDecomposition! This is required "
												 "to get correct global values!"
											  << std::endl;
						forceLatchingToLinkedCellsGrid = true;
					}
					Log::global_log->info() << "Using skin = " << skin << " for the GeneralDomainDecomposition." << std::endl;
				} else {
					Log::global_log->error() << "Datastructure section missing" << std::endl;
					Simulation::exit(1);
				}
				if(not xmlconfig.changecurrentnode("../parallelisation")){
					Log::global_log->error() << "Could not go back to parallelisation path. Aborting." << std::endl;
					Simulation::exit(1);
				}
				delete _domainDecomposition;
				_domainDecomposition = new GeneralDomainDecomposition(getcutoffRadius() + skin, _domain, forceLatchingToLinkedCellsGrid);
			} else {
				Log::global_log->error() << "Unknown parallelisation type: " << parallelisationtype << std::endl;
				Simulation::exit(1);
			}
		#else /* serial */
			if(parallelisationtype != "DummyDecomposition") {
				Log::global_log->warning()
						<< "Executable was compiled without support for parallel execution: "
						<< parallelisationtype
						<< " not available. Using serial mode." << std::endl;
				//Simulation::exit(1);
			}
			//_domainDecomposition = new DomainDecompBase();  // already set in initialize()
		#endif
			_domainDecomposition->readXML(xmlconfig);
        #ifdef ENABLE_MPI
            if(auto kdd = dynamic_cast<KDDecomposition*>(_domainDecomposition)) {
                kdd->init(_domain);
			}
        #endif

			std::string loadTimerStr("SIMULATION_FORCE_CALCULATION");
			xmlconfig.getNodeValue("timerForLoad", loadTimerStr);
			Log::global_log->info() << "Using timer " << loadTimerStr << " for the load calculation." << std::endl;
			_timerForLoad = timers()->getTimer(loadTimerStr);
			if (not _timerForLoad) {
				Log::global_log->error() << "'timerForLoad' set to a timer that does not exist('" << loadTimerStr
									<< "')! Aborting!" << std::endl;
				exit(1);
			}

			std::size_t timerForLoadAveragingLength{1ul};
			xmlconfig.getNodeValue("timerForLoad_AveragingLength", timerForLoadAveragingLength);
			Log::global_log->info() << "Averaging over " << timerForLoadAveragingLength
							   << "time steps for the load calculation." << std::endl;
			if(timerForLoadAveragingLength < 1ul) {
				Log::global_log->fatal() << "timerForLoadAveragingLength has to be at least 1" << std::endl;
				Simulation::exit(15843);
			}
			_lastTraversalTimeHistory.setCapacity(timerForLoadAveragingLength);

			xmlconfig.changecurrentnode("..");
		}
		else {
		#ifdef ENABLE_MPI
			Log::global_log->error() << "Parallelisation section missing." << std::endl;
			Simulation::exit(1);
		#else /* serial */
			// set _timerForLoad, s.t. it always exists.
			_timerForLoad = timers()->getTimer("SIMULATION_COMPUTATION");
			//_domainDecomposition = new DomainDecompBase(); // already set in initialize()
		#endif
		}

		/* datastructure */
		if(xmlconfig.changecurrentnode("datastructure")) {
			std::string datastructuretype;
			xmlconfig.getNodeValue("@type", datastructuretype);
			Log::global_log->info() << "Datastructure type: " << datastructuretype << std::endl;
			if(datastructuretype == "LinkedCells") {
#ifdef MARDYN_AUTOPAS
				Log::global_log->fatal()
					<< "LinkedCells not compiled (use AutoPas instead, or compile with disabled autopas mode)!"
					<< std::endl;
				Simulation::exit(33);
#else
				_moleculeContainer = new LinkedCells();
				/** @todo Review if we need to know the max cutoff radius usable with any datastructure. */
				Log::global_log->info() << "Setting cell cutoff radius for linked cell datastructure to " << _cutoffRadius << std::endl;
				_moleculeContainer->setCutoff(_cutoffRadius);
#endif
			} else if(datastructuretype == "AdaptiveSubCells") {
				Log::global_log->warning() << "AdaptiveSubCells no longer supported." << std::endl;
				Simulation::exit(-1);
			} else if(datastructuretype == "AutoPas" || datastructuretype == "AutoPasContainer") {
#ifdef MARDYN_AUTOPAS
				Log::global_log->info() << "Using AutoPas container." << std::endl;
				_moleculeContainer = new AutoPasContainer(_cutoffRadius);
				Log::global_log->info() << "Setting cell cutoff radius for AutoPas container to " << _cutoffRadius << std::endl;
#else
				Log::global_log->fatal() << "AutoPas not compiled (use LinkedCells instead, or compile with enabled autopas mode)!" << std::endl;
				Simulation::exit(33);
#endif
			}
			else {
				Log::global_log->error() << "Unknown data structure type: " << datastructuretype << std::endl;
				Simulation::exit(1);
			}
			_moleculeContainer->readXML(xmlconfig);

			double bBoxMin[3];
			double bBoxMax[3];
			_domainDecomposition->getBoundingBoxMinMax(_domain, bBoxMin, bBoxMax);
			bool sendTogether = _moleculeContainer->rebuild(bBoxMin, bBoxMax);
			_domainDecomposition->updateSendLeavingWithCopies(sendTogether);
			xmlconfig.changecurrentnode("..");
		} else {
			Log::global_log->error() << "Datastructure section missing" << std::endl;
			Simulation::exit(1);
		}

		// TODO: move parts to readXML in TemperatureControl?
		/* thermostats */
		if(xmlconfig.changecurrentnode("thermostats")) {
			long numThermostats = 0;
			XMLfile::Query query = xmlconfig.query("thermostat");
			numThermostats = query.card();
			Log::global_log->info() << "Number of thermostats: " << numThermostats << std::endl;
			if(numThermostats > 1) {
				Log::global_log->info() << "Enabling component wise thermostat" << std::endl;
				_velocityScalingThermostat.enableComponentwise();
			}
			std::string oldpath = xmlconfig.getcurrentnodepath();
			XMLfile::Query::const_iterator thermostatIter;
			for( thermostatIter = query.begin(); thermostatIter; thermostatIter++ ) {
				xmlconfig.changecurrentnode( thermostatIter );
				std::string thermostattype;
				xmlconfig.getNodeValue("@type", thermostattype);
				if(thermostattype == "VelocityScaling") {
					double temperature = _ensemble->T();
					xmlconfig.getNodeValue("temperature", temperature);
					std::string componentName("global");
					xmlconfig.getNodeValue("@componentId", componentName);
					if(componentName == "global"){
						_domain->setGlobalTemperature(temperature);
						Log::global_log->info() << "Adding global velocity scaling thermostat, T = " << temperature << std::endl;
					}
					else {
						int componentId = 0;
						componentId = getEnsemble()->getComponent(componentName)->ID();
						int thermostatID = _domain->getThermostat(componentId);
						_domain->setTargetTemperature(thermostatID, temperature);
						Log::global_log->info() << "Adding velocity scaling thermostat for component '" << componentName << "' (ID: " << componentId << "), T = " << temperature << std::endl;
					}
				}
				else if(thermostattype == "TemperatureControl") {
                    if (nullptr == _temperatureControl) {
                        _temperatureControl = new TemperatureControl();
                        _temperatureControl->readXML(xmlconfig);
                    } else {
                        Log::global_log->error() << "Instance of TemperatureControl allready exist! Programm exit ..."
                                            << std::endl;
                        Simulation::exit(-1);
                    }
                }
				else
				{
					Log::global_log->warning() << "Unknown thermostat " << thermostattype << std::endl;
					continue;
				}
			}
			xmlconfig.changecurrentnode(oldpath);
			xmlconfig.changecurrentnode("..");
		}
		else {
			Log::global_log->warning() << "Thermostats section missing." << std::endl;
		}

		/* long range correction */
		if(xmlconfig.changecurrentnode("longrange") )
		{
			std::string type;
			if( !xmlconfig.getNodeValue("@type", type) )
			{
				Log::global_log->error() << "LongRangeCorrection: Missing type specification. Program exit ..." << std::endl;
				Simulation::exit(-1);
			}
			if("planar" == type)
			{
				unsigned int nSlabs = 10;
				delete _longRangeCorrection;
				Log::global_log->info() << "Initializing planar LRC." << std::endl;
				_longRangeCorrection = new Planar(_cutoffRadius, _LJCutoffRadius, _domain, _domainDecomposition, _moleculeContainer, nSlabs, global_simulation);
				_longRangeCorrection->readXML(xmlconfig);
			}
			else if("homogeneous" == type)
			{
				delete _longRangeCorrection;
				Log::global_log->info() << "Initializing homogeneous LRC." << std::endl;
				_longRangeCorrection = new Homogeneous(_cutoffRadius, _LJCutoffRadius, _domain, global_simulation);
			}
			else if("none" == type)
			{
				delete _longRangeCorrection;
				_longRangeCorrection = new NoLRC(_cutoffRadius, _LJCutoffRadius, _domain, global_simulation);
			}
			else
			{
				Log::global_log->error() << "LongRangeCorrection: Wrong type. Expected type == homogeneous|planar|none. Program exit ..." << std::endl;
                Simulation::exit(-1);
			}
			xmlconfig.changecurrentnode("..");
		} else {
			delete _longRangeCorrection;
			Log::global_log->info() << "Initializing default homogeneous LRC, as no LRC was defined." << std::endl;
			_longRangeCorrection = new Homogeneous(_cutoffRadius, _LJCutoffRadius, _domain, global_simulation);
		}

		xmlconfig.changecurrentnode(".."); /* algorithm section */
	}
	else {
		Log::global_log->error() << "Algorithm section missing." << std::endl;
	}

	Log::global_log -> info() << "Registering default plugins..." << std::endl;
    // REGISTERING/ENABLING PLUGINS
	PluginFactory<PluginBase> pluginFactory;
    pluginFactory.registerDefaultPlugins();
	Log::global_log -> info() << "Successfully registered plugins." << std::endl;

	long numPlugs = 0;
	numPlugs += pluginFactory.enablePlugins(_plugins, xmlconfig, "plugin", _domain);
	numPlugs += pluginFactory.enablePlugins(_plugins, xmlconfig, "output/outputplugin", _domain);
	Log::global_log -> info() << "Number of enabled Plugins: " << numPlugs << std::endl;

	Log::global_log -> info() << "Registering callbacks." << std::endl;
	for(auto&& plugin : _plugins) {
		plugin->registerCallbacks(_callbacks);
	}
	Log::global_log -> info() << _callbacks.size() << " callbacks registered." << std::endl;

	Log::global_log -> info() << "Accessing callbacks." << std::endl;
	for(auto&& plugin : _plugins) {
		plugin->accessAllCallbacks(_callbacks);
	}
	Log::global_log -> info() << "Accessed callbacks." << std::endl;


	std::string oldpath = xmlconfig.getcurrentnodepath();

	if(xmlconfig.changecurrentnode("ensemble/phasespacepoint/file")) {
		Log::global_log->info() << "Reading phase space from file." << std::endl;
		std::string pspfiletype;
		xmlconfig.getNodeValue("@type", pspfiletype);
		Log::global_log->info() << "Phase space file type: " << pspfiletype << std::endl;

		if (pspfiletype == "ASCII") {
			_inputReader = new ASCIIReader();
			_inputReader->readXML(xmlconfig);
		}
		else if (pspfiletype == "binary") {
			_inputReader = new BinaryReader();
			_inputReader->readXML(xmlconfig);
			//!@todo read header should be either part of readPhaseSpace or readXML.
			double timestepLength = 0.005;  // <-- TODO: should be removed from parameter list
			_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
		}
#ifdef ENABLE_ADIOS2
        else if (pspfiletype == "adios2") {
			_inputReader = new Adios2Reader();
			_inputReader->readXML(xmlconfig);
		}
#endif
		else {
			Log::global_log->error() << "Unknown phase space file type" << std::endl;
			Simulation::exit(-1);
		}
	}
	xmlconfig.changecurrentnode(oldpath);

	oldpath = xmlconfig.getcurrentnodepath();
	if(xmlconfig.changecurrentnode("ensemble/phasespacepoint/generator")) {
		std::string generatorName;
		xmlconfig.getNodeValue("@name", generatorName);
		Log::global_log->info() << "Initializing phase space using generator: " << generatorName << std::endl;
		if(generatorName == "MultiObjectGenerator") {
			_inputReader = new MultiObjectGenerator();
		}
		else if(generatorName == "mkesfera") {
			_inputReader = new MkesferaGenerator();
		}
		else if(generatorName == "mkTcTS") {
			_inputReader = new MkTcTSGenerator();
		}
		else if (generatorName == "CubicGridGenerator") {
			_inputReader = new CubicGridGeneratorInternal();
		}
		else if (generatorName == "ReplicaGenerator") {
			_inputReader = new ReplicaGenerator();
		}
		else if (generatorName == "PerCellGenerator") {
			_inputReader = new PerCellGenerator();
		}
		else {
			Log::global_log->error() << "Unknown generator: " << generatorName << std::endl;
			Simulation::exit(1);
		}
		_inputReader->readXML(xmlconfig);
	}
	xmlconfig.changecurrentnode(oldpath);

	oldpath = xmlconfig.getcurrentnodepath();

	if(xmlconfig.changecurrentnode("ensemble/phasespacepoint")) {
		bool ignoreCheckpointTime = false;
		if(xmlconfig.getNodeValue("ignoreCheckpointTime", ignoreCheckpointTime)) {
			if(ignoreCheckpointTime)
				_simulationTime = 0;
		}
	}

	xmlconfig.changecurrentnode(oldpath);

	/** Prepare start options, affecting behavior of method prepare_start() */
	_prepare_start_opt.refreshIDs = false;

	oldpath = xmlconfig.getcurrentnodepath();
	if(xmlconfig.changecurrentnode("options")) {
		unsigned long numOptions = 0;
		XMLfile::Query query = xmlconfig.query("option");
		numOptions = query.card();
		Log::global_log->info() << "Number of prepare start options: " << numOptions << std::endl;

		XMLfile::Query::const_iterator optionIter;
		for( optionIter = query.begin(); optionIter; optionIter++ ) {
			xmlconfig.changecurrentnode(optionIter);
			std::string strOptionName;
			xmlconfig.getNodeValue("@name", strOptionName);
			if(strOptionName == "refreshIDs") {
				bool bVal = false;
				xmlconfig.getNodeValue(".", bVal);
				_prepare_start_opt.refreshIDs = bVal;
				if(_prepare_start_opt.refreshIDs)
					Log::global_log->info() << "Particle IDs will be refreshed before simulation start." << std::endl;
				else
					Log::global_log->info() << "Particle IDs will NOT be refreshed before simulation start." << std::endl;
			}
			else
			{
				Log::global_log->warning() << "Unknown option '" << strOptionName << "'" << std::endl;
			}
		}
	}
	xmlconfig.changecurrentnode(oldpath);
}


void Simulation::readConfigFile(std::string filename) {
	std::string extension(getFileExtension(filename.c_str()));
	Log::global_log->debug() << "Found config filename extension: " << extension << std::endl;
	if (extension == "xml") {
		initConfigXML(filename);
	}
	else {
		Log::global_log->error() << "Unknown config file extension '" << extension << "'." << std::endl;
		Simulation::exit(1);
	}
}

void Simulation::initConfigXML(const std::string& inputfilename) {
	Log::global_log->info() << "Initializing XML config file: " << inputfilename << std::endl;

	try{
		XMLfileUnits inp(inputfilename);

		Log::global_log->debug() << "Input XML:" << std::endl << std::string(inp) << std::endl;

		if(inp.changecurrentnode("/mardyn") < 0) {
			Log::global_log->error() << "Cound not find root node /mardyn in XML input file." << std::endl;
			Log::global_log->fatal() << "Not a valid MarDyn XML input file." << std::endl;
			Simulation::exit(1);
		}

		std::string version("unknown");
		inp.getNodeValue("@version", version);
		Log::global_log->info() << "MarDyn XML config file version: " << version << std::endl;

		Log::global_log->set_log_level(inp.getNodeValue_string("loglevel", "info"));

		if(inp.changecurrentnode("simulation")) {
			readXML(inp);
			inp.changecurrentnode("..");
		} // simulation-section
		else {
			Log::global_log->error() << "Simulation section missing" << std::endl;
			Simulation::exit(1);
		}
	} catch (const std::exception& e) {
		Log::global_log->error() << "Error in XML config. Please check your input file!" << std::endl;
		Log::global_log->error() << "Exception: " << e.what() << std::endl;
		Simulation::exit(7);
	}

#ifdef ENABLE_MPI
	// if we are using the DomainDecomposition, please complete its initialization:
	DomainDecomposition *temp = dynamic_cast<DomainDecomposition*>(_domainDecomposition);
	if (temp != nullptr) {
		temp->initCommunicationPartners(_cutoffRadius, _domain, _moleculeContainer);
	}
#endif

	// read particle data (or generate particles, if a generator is chosen)
	timers()->registerTimer("PHASESPACE_CREATION",  std::vector<std::string>{"SIMULATION_IO"}, new Timer());
	timers()->setOutputString("PHASESPACE_CREATION", "Phasespace creation took:");
	timers()->getTimer("PHASESPACE_CREATION")->start();

	_inputReader->readPhaseSpace(_moleculeContainer, _domain, _domainDecomposition);
	timers()->getTimer("PHASESPACE_CREATION")->stop();

	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	unsigned long globalNumMolecules = _domain->getglobalNumMolecules(true, _moleculeContainer, _domainDecomposition);
	double rho_global = globalNumMolecules / _ensemble->V();
	Log::global_log->info() << "Setting domain class parameters: N_global: " << globalNumMolecules
					   << ", rho_global: " << rho_global << ", T_global: " << _ensemble->T() << std::endl;
	_domain->setGlobalTemperature(_ensemble->T());
	_domain->setglobalRho(rho_global);

	_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
	//domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

    _ensemble->initConfigXML(_moleculeContainer);
}

void Simulation::updateForces() {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		// iterate over all particles in case we are not using the Full Shell method.
		for (auto i = _moleculeContainer->iterator(ParticleIterator::ALL_CELLS); i.isValid(); ++i){
			i->calcFM();
		}
	} // end pragma omp parallel
}

void Simulation::prepare_start() {
	Log::global_log->info() << "Initializing simulation" << std::endl;

	Log::global_log->info() << "Initialising cell processor" << std::endl;
	if (!_legacyCellProcessor) {
#ifndef ENABLE_REDUCED_MEMORY_MODE
		Log::global_log->info() << "Using vectorized cell processor." << std::endl;
		_cellProcessor = new VectorizedCellProcessor( *_domain, _cutoffRadius, _LJCutoffRadius);
#else
		Log::global_log->info() << "Using reduced memory mode (RMM) cell processor." << std::endl;
		_cellProcessor = new VCP1CLJRMM( *_domain, _cutoffRadius, _LJCutoffRadius);
#endif
	} else {
		Log::global_log->info() << "Using legacy cell processor." << std::endl;
		_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _particlePairsHandler);
	}

	if (_FMM != nullptr) {

		double globalLength[3];
		for (int i = 0; i < 3; i++) {
			globalLength[i] = _domain->getGlobalLength(i);
		}
		double bBoxMin[3];
		double bBoxMax[3];
		for (int i = 0; i < 3; i++) {
			bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
			bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
		}
		_FMM->init(globalLength, bBoxMin, bBoxMax, _moleculeContainer->getCellLength(), _moleculeContainer);

		delete _cellProcessor;
		_cellProcessor = new bhfmm::VectorizedLJP2PCellProcessor(*_domain, _LJCutoffRadius, _cutoffRadius);
	}

#ifdef ENABLE_MPI
	if(auto *kdd = dynamic_cast<KDDecomposition*>(_domainDecomposition); kdd != nullptr){
		kdd->fillTimeVecs(&_cellProcessor);
	}
#endif

	Log::global_log->info() << "Clearing halos" << std::endl;
	_moleculeContainer->deleteOuterParticles();
	Log::global_log->info() << "Updating domain decomposition" << std::endl;

	if(getMemoryProfiler()) {
		getMemoryProfiler()->doOutput("without halo copies");
	}

	updateParticleContainerAndDecomposition(1.0, false);

	if(getMemoryProfiler()) {
		getMemoryProfiler()->doOutput("with halo copies");
	}

#ifdef ENABLE_REDUCED_MEMORY_MODE
	// the leapfrog integration requires that we move the velocities by one half-timestep
	// so we halve vcp1clj_wr_cellProcessor::_dtInvm
	VCP1CLJRMM * vcp1clj_wr_cellProcessor = static_cast<VCP1CLJRMM * >(_cellProcessor);
	double dt_inv_m = vcp1clj_wr_cellProcessor->getDtInvm();
	vcp1clj_wr_cellProcessor->setDtInvm(dt_inv_m * 0.5);
#endif /* ENABLE_REDUCED_MEMORY_MODE */

	Log::global_log->info() << "Performing initial force calculation" << std::endl;
	global_simulation->timers()->start("SIMULATION_FORCE_CALCULATION");
	_moleculeContainer->traverseCells(*_cellProcessor);

	global_simulation->timers()->stop("SIMULATION_FORCE_CALCULATION");

	if (_longRangeCorrection != nullptr) {
		Log::global_log->info() << "Initializing LongRangeCorrection" << std::endl;
		_longRangeCorrection->init();
	} else {
		Log::global_log->fatal() << "No _longRangeCorrection set!" << std::endl;
		Simulation::exit(93742);
	}
	// longRangeCorrection is a site-wise force plugin, so we have to call it before updateForces()
	_longRangeCorrection->calculateLongRange();

	// Update forces in molecules so they can be exchanged - future
	updateForces();

	// Exchange forces if it's required by the cell container.
	if(_moleculeContainer->requiresForceExchange()){
		_domainDecomposition->exchangeForces(_moleculeContainer, _domain);
	}

#ifdef ENABLE_REDUCED_MEMORY_MODE
	// now set vcp1clj_wr_cellProcessor::_dtInvm back.
	vcp1clj_wr_cellProcessor->setDtInvm(dt_inv_m);
#endif /* ENABLE_REDUCED_MEMORY_MODE */

	_loopCompTime = global_simulation->timers()->getTime("SIMULATION_FORCE_CALCULATION");
	global_simulation->timers()->reset("SIMULATION_FORCE_CALCULATION");
	++_loopCompTimeSteps;

	if (_FMM != nullptr) {
		Log::global_log->info() << "Performing initial FMM force calculation" << std::endl;
		_FMM->computeElectrostatics(_moleculeContainer);
	}

	/** Init TemperatureControl beta_trans, beta_rot log-files, register as observer if plugin DistControl is in use. */
	if(nullptr != _temperatureControl)
		_temperatureControl->prepare_start();  // Has to be called before plugin initialization (see below): plugin->init(...)

	// initializing plugins and starting plugin timers
	for (auto& plugin : _plugins) {
		Log::global_log->info() << "Initializing plugin " << plugin->getPluginName() << std::endl;
		plugin->init(_moleculeContainer, _domainDecomposition, _domain);
		std::string timer_name = plugin->getPluginName();
		// TODO: real timer
		global_simulation->timers()->registerTimer(timer_name, std::vector<std::string>{"SIMULATION_PER_STEP_IO"}, new Timer());
		std::string timer_plugin_string = std::string("Plugin ") + timer_name + std::string(" took:");
		global_simulation->timers()->setOutputString(timer_name, timer_plugin_string);
	}

	//afterForces Plugin Call
	Log::global_log->debug() << "[AFTER FORCES] Performing AfterForces plugin call"
						<< std::endl;
	for (auto plugin : _plugins) {
		Log::global_log->debug() << "[AFTER FORCES] Plugin: "
							<< plugin->getPluginName() << std::endl;
		plugin->afterForces(_moleculeContainer, _domainDecomposition, _simstep);
	}

#ifndef MARDYN_AUTOPAS
	// clear halo
	Log::global_log->info() << "Clearing halos" << std::endl;
	_moleculeContainer->deleteOuterParticles();
#endif

	// here we have to call calcFM() manually, otherwise force and moment are not
	// updated inside the molecule (actually this is done in upd_postF)
	// integrator->eventForcesCalculated should not be called, since otherwise the velocities would already be updated.
	//updateForces();

	Log::global_log->info() << "Calculating global values" << std::endl;
	_domain->calculateThermostatDirectedVelocity(_moleculeContainer);

	_domain->calculateVelocitySums(_moleculeContainer);

	_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer,
			true, 1.0);
	Log::global_log->debug() << "Calculating global values finished." << std::endl;

	_ensemble->prepare_start();

	_simstep = _initSimulation = (unsigned long) round(_simulationTime / _integrator->getTimestepLength() );
	Log::global_log->info() << "Set initial time step to start from to " << _initSimulation << std::endl;
	Log::global_log->info() << "System initialised with " << _domain->getglobalNumMolecules(true, _moleculeContainer, _domainDecomposition) << " molecules." << std::endl;

	/** refresh particle IDs */
	if(_prepare_start_opt.refreshIDs)
		this->refreshParticleIDs();
}

void Simulation::simulate() {
	Log::global_log->info() << "Started simulation" << std::endl;

	_ensemble->updateGlobalVariable(_moleculeContainer, NUM_PARTICLES);
	Log::global_log->debug() << "Number of particles in the Ensemble: " << _ensemble->N() << std::endl;
// 	ensemble.updateGlobalVariable(ENERGY);
// 	Log::global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E() << std::endl;
// 	ensemble.updateGlobalVariable(TEMPERATURE);
// 	Log::global_log->debug() << "Temperature of the Ensemble: " << ensemble.T() << std::endl;*/

	/***************************************************************************/
	/* BEGIN MAIN LOOP                                                         */
	/***************************************************************************/

	global_simulation->timers()->setOutputString("SIMULATION_LOOP", "Computation in main loop took:");
	global_simulation->timers()->setOutputString("SIMULATION_DECOMPOSITION", "Decomposition took:");
	global_simulation->timers()->setOutputString("SIMULATION_COMPUTATION", "Computation took:");
	global_simulation->timers()->setOutputString("SIMULATION_PER_STEP_IO", "IO in main loop took:");
	global_simulation->timers()->setOutputString("SIMULATION_FORCE_CALCULATION", "Force calculation took:");
	global_simulation->timers()->setOutputString("SIMULATION_MPI_OMP_COMMUNICATION", "Communication took:");
	global_simulation->timers()->setOutputString("SIMULATION_UPDATE_CONTAINER", "Container update took:");
	global_simulation->timers()->setOutputString("SIMULATION_UPDATE_CACHES", "Cache update took:");
	global_simulation->timers()->setOutputString("COMMUNICATION_PARTNER_INIT_SEND", "initSend() took:");
	global_simulation->timers()->setOutputString("COMMUNICATION_PARTNER_TEST_RECV", "testRecv() took:");

	// all timers except the ioTimer measure inside the main loop
	Timer* loopTimer = global_simulation->timers()->getTimer("SIMULATION_LOOP"); ///< timer for the entire simulation loop (synced)
	Timer* decompositionTimer = global_simulation->timers()->getTimer("SIMULATION_DECOMPOSITION"); ///< timer for decomposition: sub-timer of loopTimer
	Timer* computationTimer = global_simulation->timers()->getTimer("SIMULATION_COMPUTATION"); ///< timer for computation: sub-timer of loopTimer
	Timer* perStepIoTimer = global_simulation->timers()->getTimer("SIMULATION_PER_STEP_IO"); ///< timer for io in simulation loop: sub-timer of loopTimer
	Timer* forceCalculationTimer = global_simulation->timers()->getTimer("SIMULATION_FORCE_CALCULATION"); ///< timer for force calculation: sub-timer of computationTimer
	Timer* mpiOMPCommunicationTimer = global_simulation->timers()->getTimer("SIMULATION_MPI_OMP_COMMUNICATION"); ///< timer for measuring MPI-OMP communication time: sub-timer of decompositionTimer

	//loopTimer->set_sync(true);
	//global_simulation->timers()->setSyncTimer("SIMULATION_LOOP", true);
#ifdef WITH_PAPI
	const char *papi_event_list[] = { "PAPI_TOT_CYC", "PAPI_TOT_INS" /*, "PAPI_VEC_DP", "PAPI_L2_DCM", "PAPI_L2_ICM", "PAPI_L1_ICM", "PAPI_DP_OPS", "PAPI_VEC_INS" }; */
	int num_papi_events = sizeof(papi_event_list) / sizeof(papi_event_list[0]);
	loopTimer->add_papi_counters(num_papi_events, (char**) papi_event_list);
#endif
	loopTimer->start();
#ifndef NDEBUG
#ifndef ENABLE_MPI
		unsigned particleNoTest;
#endif
#endif
	if(getMemoryProfiler()) {
		getMemoryProfiler()->doOutput();
	}

	pluginEndStepCall(_initSimulation);


	// keepRunning() increments the simstep counter before the first iteration
	_simstep = _initSimulation;

	// stores the timing info for the previous load. This is used for the load calculation and the rebalancing.
	double previousTimeForLoad = 0.;

	while (keepRunning()) {
		Log::global_log->debug() << "timestep: " << getSimulationStep() << std::endl;
		Log::global_log->debug() << "simulation time: " << getSimulationTime() << std::endl;
		global_simulation->timers()->incrementTimerTimestepCounter();

		computationTimer->start();

        // beforeEventNewTimestep Plugin Call
        Log::global_log -> debug() << "[BEFORE EVENT NEW TIMESTEP] Performing beforeEventNewTimestep plugin call" << std::endl;
        for (auto plugin : _plugins) {
            Log::global_log -> debug() << "[BEFORE EVENT NEW TIMESTEP] Plugin: " << plugin->getPluginName() << std::endl;
            plugin->beforeEventNewTimestep(_moleculeContainer, _domainDecomposition, _simstep);
        }

        _ensemble->beforeEventNewTimestep(_moleculeContainer, _domainDecomposition, _simstep);

		_integrator->eventNewTimestep(_moleculeContainer, _domain);

        // beforeForces Plugin Call
        Log::global_log -> debug() << "[BEFORE FORCES] Performing BeforeForces plugin call" << std::endl;
        for (auto plugin : _plugins) {
            Log::global_log -> debug() << "[BEFORE FORCES] Plugin: " << plugin->getPluginName() << std::endl;
            plugin->beforeForces(_moleculeContainer, _domainDecomposition, _simstep);
        }

		computationTimer->stop();



#if defined(ENABLE_MPI)
		bool overlapCommComp = _overlappingP2P;
#else
		bool overlapCommComp = false;
#endif

		double startEtime = computationTimer->get_etime();
		if (overlapCommComp) {
			double currentTime = _timerForLoad->get_etime();
			performOverlappingDecompositionAndCellTraversalStep(currentTime - previousTimeForLoad);
			previousTimeForLoad = currentTime;
			// Force timer and computation timer are running at this point!
		}
		else {
			decompositionTimer->start();
			// ensure that all Particles are in the right cells and exchange Particles
			Log::global_log->debug() << "Updating container and decomposition" << std::endl;

			double currentTime = _timerForLoad->get_etime();
			updateParticleContainerAndDecomposition(currentTime - previousTimeForLoad, true);
			previousTimeForLoad = currentTime;

			decompositionTimer->stop();

			// Force calculation and other pair interaction related computations
			Log::global_log->debug() << "Traversing pairs" << std::endl;
			computationTimer->start();
			forceCalculationTimer->start();

			_moleculeContainer->traverseCells(*_cellProcessor);
			// Force timer and computation timer are running at this point!
		}

		// siteWiseForces Plugin Call
		Log::global_log -> debug() << "[SITEWISE FORCES] Performing siteWiseForces plugin call" << std::endl;
		for (auto plugin : _plugins) {
			Log::global_log -> debug() << "[SITEWISE FORCES] Plugin: " << plugin->getPluginName() << std::endl;
			plugin->siteWiseForces(_moleculeContainer, _domainDecomposition, _simstep);
		}

		// longRangeCorrection is a site-wise force plugin, so we have to call it before updateForces()
		_longRangeCorrection->calculateLongRange();

		// Update forces in molecules so they can be exchanged
		updateForces();

		forceCalculationTimer->stop();
		computationTimer->stop();

		decompositionTimer->start();
		// Exchange forces if it's required by the cell container.
		if(_moleculeContainer->requiresForceExchange()){
			Log::global_log->debug() << "Exchanging Forces" << std::endl;
			_domainDecomposition->exchangeForces(_moleculeContainer, _domain);
		}
		decompositionTimer->stop();
		_loopCompTime += computationTimer->get_etime() - startEtime;
		_loopCompTimeSteps ++;

		computationTimer->start();


		if (_FMM != nullptr) {
			Log::global_log->debug() << "Performing FMM calculation" << std::endl;
			_FMM->computeElectrostatics(_moleculeContainer);
		}

		//afterForces Plugin Call
		Log::global_log -> debug() << "[AFTER FORCES] Performing AfterForces plugin call" << std::endl;
		for (auto plugin : _plugins) {
			Log::global_log -> debug() << "[AFTER FORCES] Plugin: " << plugin->getPluginName() << std::endl;
			plugin->afterForces(_moleculeContainer, _domainDecomposition, _simstep);
		}

		_ensemble->afterForces(_moleculeContainer, _domainDecomposition, _cellProcessor, _simstep);

		// TODO: test deletions and insertions
		Log::global_log->debug() << "Deleting outer particles / clearing halo." << std::endl;
#ifndef MARDYN_AUTOPAS
		_moleculeContainer->deleteOuterParticles();
#endif

		if (!(_simstep % _collectThermostatDirectedVelocity)) {
			_domain->calculateThermostatDirectedVelocity(_moleculeContainer);
		}

		_longRangeCorrection->writeProfiles(_domainDecomposition, _domain, _simstep);

		_ensemble->beforeThermostat(_simstep, _initStatistics);

		Log::global_log->debug() << "Inform the integrator (forces calculated)" << std::endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);

		// calculate the global macroscopic values from the local values
		Log::global_log->debug() << "Calculate macroscopic values" << std::endl;
		_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer,
				(!(_simstep % _collectThermostatDirectedVelocity)), Tfactor(_simstep));

		// scale velocity and angular momentum
        // TODO: integrate into Temperature Control
		if ( !_domain->NVE() && _temperatureControl == nullptr) {
			if (_thermostatType ==VELSCALE_THERMOSTAT) {
				Log::global_log->debug() << "Velocity scaling" << std::endl;
				if (_domain->severalThermostats()) {
					_velocityScalingThermostat.enableComponentwise();
					for(unsigned int cid = 0; cid < global_simulation->getEnsemble()->getComponents()->size(); cid++) {
						int thermostatId = _domain->getThermostat(cid);
						_velocityScalingThermostat.setBetaTrans(thermostatId, _domain->getGlobalBetaTrans(thermostatId));
						_velocityScalingThermostat.setBetaRot(thermostatId, _domain->getGlobalBetaRot(thermostatId));
						Log::global_log->debug() << "Thermostat for CID: " << cid << " thermID: " << thermostatId
								<< " B_trans: " << _velocityScalingThermostat.getBetaTrans(thermostatId)
								<< " B_rot: " << _velocityScalingThermostat.getBetaRot(thermostatId) << std::endl;
						double v[3];
						for(int d = 0; d < 3; d++) {
							v[d] = _domain->getThermostatDirectedVelocity(thermostatId, d);
						}
						_velocityScalingThermostat.setVelocity(thermostatId, v);
					}
				}
				else {
					_velocityScalingThermostat.setGlobalBetaTrans(_domain->getGlobalBetaTrans());
					_velocityScalingThermostat.setGlobalBetaRot(_domain->getGlobalBetaRot());
					/* TODO */
					// Undirected global thermostat not implemented!
				}
				_velocityScalingThermostat.apply(_moleculeContainer);


			}
		} else if ( _temperatureControl != nullptr) {
			// mheinen 2015-07-27 --> TEMPERATURE_CONTROL
           _temperatureControl->DoLoopsOverMolecules(_domainDecomposition, _moleculeContainer, _simstep);
        }
        // <-- TEMPERATURE_CONTROL



		advanceSimulationTime(_integrator->getTimestepLength());

		/* BEGIN PHYSICAL SECTION:
		 * the system is in a consistent state so we can extract global variables
		 */
		//! @todo the number of particles per component stored in components has to be
		//!       updated here in case we insert/remove particles
		// _ensemble->updateGlobalVariable(NUM_PARTICLES);
		// Log::global_log->debug() << "Number of particles in the Ensemble: " << _ensemble->N() << std::endl;
		/*
		ensemble.updateGlobalVariable(ENERGY);
		Log::global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E() << std::endl;
		ensemble.updateGlobalVariable(TEMPERATURE);
		Log::global_log->debug() << "Temperature of the Ensemble: " << ensemble.T() << std::endl;
		*/
		/* END PHYSICAL SECTION */


		computationTimer->stop();
		perStepIoTimer->start();

		// CALL ALL PLUGIN ENDSTEP METHODS
		pluginEndStepCall(_simstep);

		if( (_forced_checkpoint_time > 0) && (loopTimer->get_etime() >= _forced_checkpoint_time) ) {
			/* force checkpoint for specified time */
			std::string cpfile(_outputPrefix + ".timed.restart.dat");
			Log::global_log->info() << "Writing timed, forced checkpoint to file '" << cpfile << "'" << std::endl;
			_domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition, _simulationTime);
			_forced_checkpoint_time = -1; /* disable for further timesteps */
		}
		perStepIoTimer->stop();
	}
	loopTimer->stop();
	/***************************************************************************/
	/* END MAIN LOOP                                                           */
	/***************************************************************************/


	global_simulation->timers()->registerTimer("SIMULATION_FINAL_IO", std::vector<std::string>{"SIMULATION_IO"}, new Timer());
	global_simulation->timers()->setOutputString("SIMULATION_FINAL_IO", "Final IO took:");
	global_simulation->timers()->getTimer("SIMULATION_FINAL_IO")->start();
    if( _finalCheckpoint ) {
        /* write final checkpoint */
        std::string cpfile(_outputPrefix + ".restart.dat");
        Log::global_log->info() << "Writing final checkpoint to file '" << cpfile << "'" << std::endl;
        _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition, _simulationTime, false);
    }

	Log::global_log->info() << "Finish plugins" << std::endl;
	for (auto plugin : _plugins) {
		plugin->finish(_moleculeContainer, _domainDecomposition, _domain);
	}
	global_simulation->timers()->getTimer("SIMULATION_FINAL_IO")->stop();

	Log::global_log->info() << "Timing information:" << std::endl;
	global_simulation->timers()->printTimers();
	global_simulation->timers()->resetTimers();
	if(getMemoryProfiler()) {
		getMemoryProfiler()->doOutput();
	}
	Log::global_log->info() << std::endl;

#ifdef WITH_PAPI
	Log::global_log->info() << "PAPI counter values for loop timer:"  << std::endl;
	for(int i = 0; i < loopTimer->get_papi_num_counters(); i++) {
		Log::global_log->info() << "  " << papi_event_list[i] << ": " << loopTimer->get_global_papi_counter(i) << std::endl;
	}
#endif /* WITH_PAPI */
}

void Simulation::pluginEndStepCall(unsigned long simstep) {

	std::list<PluginBase*>::iterator pluginIter;
	for (pluginIter = _plugins.begin(); pluginIter != _plugins.end(); pluginIter++) {
		PluginBase* plugin = (*pluginIter);
		Log::global_log->debug() << "Plugin end of step: " << plugin->getPluginName() << std::endl;
		global_simulation->timers()->start(plugin->getPluginName());
		plugin->endStep(_moleculeContainer, _domainDecomposition, _domain, simstep);
		global_simulation->timers()->stop(plugin->getPluginName());
	}


	if (_domain->thermostatWarning())
		Log::global_log->warning() << "Thermostat!" << std::endl;
	/* TODO: thermostat */
	Log::global_log->info() << "Simstep = " << simstep << "\tT = "
					   << _domain->getGlobalCurrentTemperature() << "\tU_pot = "
					   << _domain->getGlobalUpot() << "\tp = "
					   << _domain->getGlobalPressure() << std::endl;
	if (std::isnan(_domain->getGlobalCurrentTemperature()) || std::isnan(_domain->getGlobalUpot()) || std::isnan(_domain->getGlobalPressure())) {
		Log::global_log->error() << "NaN detected, exiting." << std::endl;
		Simulation::exit(1);
	}
}

void Simulation::finalize() {
	if (_FMM != nullptr) {
		_FMM->printTimers();
		auto * temp = dynamic_cast<bhfmm::VectorizedLJP2PCellProcessor*>(_cellProcessor);
		temp->printTimers();
	}

#ifdef TASKTIMINGPROFILE
	std::string outputFileName = "taskTimings_"
								 + std::to_string(std::time(nullptr))
								 + ".csv";
	_taskTimingProfiler->dump(outputFileName);
#endif /* TASKTIMINGPROFILE */

	if (_domainDecomposition != nullptr) {
		delete _domainDecomposition;
		_domainDecomposition = nullptr;
	}

	_plugins.remove_if([](PluginBase * plugin) { delete plugin; return true; });
	global_simulation = nullptr;
}

void Simulation::updateParticleContainerAndDecomposition(double lastTraversalTime, bool useTimers) {
	// The particles have moved, so the leaving + halo particles have changed and have to be communicated.
	if(not useTimers){
		global_simulation->timers()->deactivateTimer("SIMULATION_UPDATE_CONTAINER");
		global_simulation->timers()->deactivateTimer("SIMULATION_MPI_OMP_COMMUNICATION");
		global_simulation->timers()->deactivateTimer("SIMULATION_UPDATE_CACHES");
	}

	global_simulation->timers()->start("SIMULATION_UPDATE_CONTAINER");
	_moleculeContainer->update();
	global_simulation->timers()->stop("SIMULATION_UPDATE_CONTAINER");

	_lastTraversalTimeHistory.insert(lastTraversalTime);
	double averageLastTraversalTime =
		std::accumulate(_lastTraversalTimeHistory.begin(), _lastTraversalTimeHistory.end(), 0.) /
		_lastTraversalTimeHistory.size();

	bool forceRebalancing = false;
	global_simulation->timers()->start("SIMULATION_MPI_OMP_COMMUNICATION");
	_domainDecomposition->balanceAndExchange(averageLastTraversalTime, forceRebalancing, _moleculeContainer,
											 _domain);
	global_simulation->timers()->stop("SIMULATION_MPI_OMP_COMMUNICATION");

	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	global_simulation->timers()->start("SIMULATION_UPDATE_CACHES");
	_moleculeContainer->updateMoleculeCaches();
	global_simulation->timers()->stop("SIMULATION_UPDATE_CACHES");

	if(not useTimers){
		global_simulation->timers()->activateTimer("SIMULATION_UPDATE_CONTAINER");
		global_simulation->timers()->activateTimer("SIMULATION_MPI_OMP_COMMUNICATION");
		global_simulation->timers()->activateTimer("SIMULATION_UPDATE_CACHES");
	}
}

void Simulation::performOverlappingDecompositionAndCellTraversalStep(double etime) {
	bool forceRebalancing = false;

#ifdef ENABLE_MPI
	auto* dd = dynamic_cast<DomainDecompMPIBase*>(_domainDecomposition);
	if (not dd) {
		Log::global_log->fatal() << "DomainDecompMPIBase* required for overlapping comm, but dynamic_cast failed." << std::endl;
		Simulation::exit(873456);
	}
	NonBlockingMPIMultiStepHandler nonBlockingMPIHandler {dd, _moleculeContainer, _domain, _cellProcessor};

	//NonBlockingMPIHandlerBase nonBlockingMPIHandler {dd, _moleculeContainer, _domain, _cellProcessor};

	nonBlockingMPIHandler.performOverlappingTasks(forceRebalancing, etime);
#else
	Log::global_log->fatal() << "performOverlappingDecompositionAndCellTraversalStep() called with disabled MPI." << std::endl;
	Simulation::exit(873457);
#endif
}

void Simulation::setDomainDecomposition(DomainDecompBase* domainDecomposition) {
	delete _domainDecomposition;
	_domainDecomposition = domainDecomposition;
}

unsigned long Simulation::getTotalNumberOfMolecules() const {
	return _domain->N();
}

/* FIXME: we should provide a more general way of doing this */
double Simulation::Tfactor(unsigned long simstep) {
	double xi = (double) (simstep - _initSimulation) / (double) (_initCanonical - _initSimulation);
	if ((xi < 0.1) || (xi > 0.9))
		return 1.0;
	else if (xi < 0.3)
		return 15.0 * xi - 0.5;
	else if (xi < 0.4)
		return 10.0 - 20.0 * xi;
	else if (xi < 0.6)
		return 2.0;
	else
		return 4 - 10.0 * xi / 3.0;
}

void Simulation::initialize() {
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif

	global_simulation = this;

	delete _domainDecomposition;
#ifndef ENABLE_MPI
	Log::global_log->info() << "Creating alibi domain decomposition ... " << std::endl;
	_domainDecomposition = new DomainDecompBase();
#else
	Log::global_log->info() << "Creating standard domain decomposition ... " << std::endl;
	_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
#endif

	_outputPrefix.append(gettimestring());

	Log::global_log->info() << "Creating domain ..." << std::endl;
	_domain = new Domain(ownrank);
	Log::global_log->info() << "Creating ParticlePairs2PotForceAdapter ..." << std::endl;
	_particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);

	Log::global_log->info() << "Initialization done" << std::endl;
}

bool Simulation::keepRunning() {

	// Simstep Criterion
	if (_simstep >= _numberOfTimesteps){
		Log::global_log->info() << "Maximum Simstep reached: " << _simstep << std::endl;
		return false;
	}
	// WallTime Criterion, elapsed time since Simulation constructor
	else if(_wallTimeEnabled && _timeFromStart.get_etime_running() > _maxWallTime){
		Log::global_log->info() << "Maximum Walltime reached (s): " << _maxWallTime << std::endl;
		return false;
	}
	else{
		_simstep++;
		return true;
	}
}


PluginBase* Simulation::getPlugin(const std::string& name)  {
	for(auto& plugin : _plugins) {
		if(name == plugin->getPluginName()) {
			return plugin;
		}
	}
	return nullptr;
}

unsigned long Simulation::getNumberOfTimesteps() const {
	return _numberOfTimesteps;
}

CellProcessor *Simulation::getCellProcessor() const {
	return _cellProcessor;
}

void Simulation::refreshParticleIDs()
{
    unsigned long start_ID = 0;
#ifdef ENABLE_MPI
    unsigned long num_molecules_local = _moleculeContainer->getNumberOfParticles();
    MPI_Exscan(&num_molecules_local, &start_ID, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

	for (auto pit = _moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		pit->setid(++start_ID);
	}
}
