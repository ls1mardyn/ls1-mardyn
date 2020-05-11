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
#include "io/MultiObjectGenerator.h"
#include "io/TcTS.h"
#include "io/Mkesfera.h"
#include "io/CubicGridGeneratorInternal.h"
#include "io/ReplicaGenerator.h"
#include "io/TimerProfiler.h"
#include "io/MemoryProfiler.h"

#include "ensemble/GrandCanonicalEnsemble.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/CavityEnsemble.h"

#include "thermostats/VelocityScalingThermostat.h"
#include "thermostats/TemperatureControl.h"

#include "utils/FileUtils.h"
#include "utils/OptionParser.h"
#include "utils/Logger.h"

#include "longRange/LongRangeCorrection.h"
#include "longRange/Homogeneous.h"
#include "longRange/Planar.h"

#include "bhfmm/FastMultipoleMethod.h"
#include "bhfmm/cellProcessors/VectorizedLJP2PCellProcessor.h"

#include "plugins/VectorizationTuner.h"

using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;


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
		string integratorType;
		xmlconfig.getNodeValue("@type", integratorType);
		global_log->info() << "Integrator type: " << integratorType << endl;
		if(integratorType == "Leapfrog") {
#ifdef ENABLE_REDUCED_MEMORY_MODE
			global_log->error() << "The reduced memory mode (RMM) requires the LeapfrogRMM integrator." << endl;
			Simulation::exit(-1);
#endif
			_integrator = new Leapfrog();
		} else if (integratorType == "LeapfrogRMM") {
			_integrator = new LeapfrogRMM();
		} else {
			global_log-> error() << "Unknown integrator " << integratorType << endl;
			Simulation::exit(1);
		}
		_integrator->readXML(xmlconfig);
		_integrator->init();
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->error() << "Integrator section missing." << endl;
	}

	/* run section */
	if(xmlconfig.changecurrentnode("run")) {
		xmlconfig.getNodeValueReduced("currenttime", _simulationTime);
		global_log->info() << "Simulation start time: " << _simulationTime << endl;
		/* steps */
		xmlconfig.getNodeValue("equilibration/steps", _initStatistics);
		global_log->info() << "Number of equilibration steps: " << _initStatistics << endl;
		xmlconfig.getNodeValue("production/steps", _numberOfTimesteps);
		global_log->info() << "Number of timesteps: " << _numberOfTimesteps << endl;
		xmlconfig.getNodeValue("production/loop-abort-time", _maxWallTime);
		if(_maxWallTime != -1) {
			global_log->info() << "Max loop-abort-time set: " << _maxWallTime << endl;
			_wallTimeEnabled = true;
		}
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->error() << "Run section missing." << endl;
	}

	/* ensemble */
	if(xmlconfig.changecurrentnode("ensemble")) {
		string ensembletype;
		xmlconfig.getNodeValue("@type", ensembletype);
		global_log->info() << "Ensemble: " << ensembletype<< endl;
		if (ensembletype == "NVT") {
			delete _ensemble;
			_ensemble = new CanonicalEnsemble();
		} else if (ensembletype == "muVT") {
			global_log->error() << "muVT ensemble not completely implemented via XML input." << endl;
			/* TODO: GrandCanonical terminates as readXML is not implemented and it is not tested yet. */
			_ensemble = new GrandCanonicalEnsemble();
		} else {
			global_log->error() << "Unknown ensemble type: " << ensembletype << endl;
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
		global_log->error() << "Ensemble section missing." << endl;
		Simulation::exit(1);
	}

	//The mixing coefficents have to be read in the ensemble part
	//if this didn't happen fill the mixing coeffs with default values
	auto& dmixcoeff = _domain->getmixcoeff();
	const std::size_t compNum = _ensemble->getComponents()->size();
	//1 Comps: 0 coeffs; 2 Comps: 2 coeffs; 3 Comps: 6 coeffs; 4 Comps 12 coeffs
	const std::size_t neededCoeffs = compNum*(compNum-1);
	if(dmixcoeff.size() < neededCoeffs){
		global_log->warning() << "Not enough mixing coefficients were given! (Filling the missing ones with 1)" << '\n';
		global_log->warning() << "This can happen because the xml-input doesn't support these yet!" << endl;
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
				global_log->info() << "dimensionless default cutoff radius:\t" << _cutoffRadius << endl;
			}
			if(xmlconfig.getNodeValueReduced("radiusLJ", _LJCutoffRadius)) {
				global_log->info() << "dimensionless LJ cutoff radius:\t" << _LJCutoffRadius << endl;
			}
			/** @todo introduce maxCutoffRadius here for datastructures, ...
			 *        maybe use map/list to store cutoffs for different potentials? */
			_cutoffRadius = max(_cutoffRadius, _LJCutoffRadius);
			if(_cutoffRadius <= 0) {
				global_log->error() << "cutoff radius <= 0." << endl;
				Simulation::exit(1);
			}
			global_log->info() << "dimensionless cutoff radius:\t" << _cutoffRadius << endl;
			xmlconfig.changecurrentnode("..");
		} else {
			global_log->error() << "Cutoff section missing." << endl;
			Simulation::exit(1);
		}

		/* electrostatics */
		/** @todo This may be better go into a physical section for constants? */
		if(xmlconfig.changecurrentnode("electrostatic[@type='ReactionField']")) {
			double epsilonRF = 0;
			xmlconfig.getNodeValueReduced("epsilon", epsilonRF);
			global_log->info() << "Epsilon Reaction Field: " << epsilonRF << endl;
			_domain->setepsilonRF(epsilonRF);
			xmlconfig.changecurrentnode("..");
		} else {
			global_log->error() << "Electrostatics section for reaction field setup missing." << endl;
			Simulation::exit(1);
		}

		if (xmlconfig.changecurrentnode("electrostatic[@type='FastMultipoleMethod']")) {
			_FMM = new bhfmm::FastMultipoleMethod();
			_FMM->readXML(xmlconfig);
			xmlconfig.changecurrentnode("..");
		}

		/* parallelisation */
		if(xmlconfig.changecurrentnode("parallelisation")) {
			string parallelisationtype("DomainDecomposition");
			xmlconfig.getNodeValue("@type", parallelisationtype);
			global_log->info() << "Parallelisation type: " << parallelisationtype << endl;
		#ifdef ENABLE_MPI
			/// @todo Dummy Decomposition now included in DecompBase - still keep this name?
			if(parallelisationtype == "DummyDecomposition") {
				global_log->error() << "DummyDecomposition not available in parallel mode." << endl;
			}
			else if(parallelisationtype == "DomainDecomposition") {
				delete _domainDecomposition;
				_domainDecomposition = new DomainDecomposition();
			}
			else if(parallelisationtype == "KDDecomposition") {
				delete _domainDecomposition;
				_domainDecomposition = new KDDecomposition(getcutoffRadius(), _domain, _ensemble->getComponents()->size());
			} else if (parallelisationtype == "GeneralDomainDecomposition") {
				double skin = 0.;
				// we need the skin here, so we extract it from the AutoPas container's xml,
				// because the ParticleContainer needs to be instantiated later. :/
				xmlconfig.changecurrentnode("..");
				if (xmlconfig.changecurrentnode("datastructure")) {
					string datastructuretype;
					xmlconfig.getNodeValue("@type", datastructuretype);
					if (datastructuretype == "AutoPas" or datastructuretype == "AutoPasContainer") {
						xmlconfig.getNodeValue("skin", skin);
						global_log->info() << "Using skin = " << skin << " for the GeneralDomainDecomposition." << std::endl;
					} else {
						global_log->error() << "Using the GeneralDomainDecomposition is only supported when using "
											   "AutoPas, but the configuration file does not use it."
											<< endl;
						Simulation::exit(2);
					}
				} else {
					global_log->error() << "Datastructure section missing" << endl;
					Simulation::exit(1);
				}
				if(not xmlconfig.changecurrentnode("../parallelisation")){
					global_log->error() << "Could not go back to parallelisation path. Aborting." << endl;
					Simulation::exit(1);
				}
				delete _domainDecomposition;
				_domainDecomposition = new GeneralDomainDecomposition(getcutoffRadius() + skin, _domain);
			} else {
				global_log->error() << "Unknown parallelisation type: " << parallelisationtype << endl;
				Simulation::exit(1);
			}
		#else /* serial */
			if(parallelisationtype != "DummyDecomposition") {
				global_log->warning()
						<< "Executable was compiled without support for parallel execution: "
						<< parallelisationtype
						<< " not available. Using serial mode." << endl;
				//Simulation::exit(1);
			}
			//_domainDecomposition = new DomainDecompBase();  // already set in initialize()
		#endif
			_domainDecomposition->readXML(xmlconfig);

			string loadTimerStr("SIMULATION_COMPUTATION");
			xmlconfig.getNodeValue("timerForLoad", loadTimerStr);
			global_log->info() << "Using timer " << loadTimerStr << " for the load calculation." << std::endl;
			_timerForLoad = timers()->getTimer(loadTimerStr);
			if (not _timerForLoad) {
				global_log->error() << "'timerForLoad' set to a timer that does not exist('" << loadTimerStr
									<< "')! Aborting!" << std::endl;
				exit(1);
			}

			xmlconfig.changecurrentnode("..");
		}
		else {
		#ifdef ENABLE_MPI
			global_log->error() << "Parallelisation section missing." << endl;
			Simulation::exit(1);
		#else /* serial */
			// set _timerForLoad, s.t. it always exists.
			_timerForLoad = timers()->getTimer("SIMULATION_COMPUTATION");
			//_domainDecomposition = new DomainDecompBase(); // already set in initialize()
		#endif
		}

		/* datastructure */
		if(xmlconfig.changecurrentnode("datastructure")) {
			string datastructuretype;
			xmlconfig.getNodeValue("@type", datastructuretype);
			global_log->info() << "Datastructure type: " << datastructuretype << endl;
			if(datastructuretype == "LinkedCells") {
#ifdef MARDYN_AUTOPAS
				global_log->fatal()
					<< "LinkedCells not compiled (use AutoPas instead, or compile with disabled autopas mode)!"
					<< std::endl;
				Simulation::exit(33);
#else
				_moleculeContainer = new LinkedCells();
				/** @todo Review if we need to know the max cutoff radius usable with any datastructure. */
				global_log->info() << "Setting cell cutoff radius for linked cell datastructure to " << _cutoffRadius << endl;
				_moleculeContainer->setCutoff(_cutoffRadius);
#endif
			} else if(datastructuretype == "AdaptiveSubCells") {
				global_log->warning() << "AdaptiveSubCells no longer supported." << std::endl;
				Simulation::exit(-1);
			} else if(datastructuretype == "AutoPas" || datastructuretype == "AutoPasContainer") {
#ifdef MARDYN_AUTOPAS
				global_log->info() << "Using AutoPas container." << std::endl;
				_moleculeContainer = new AutoPasContainer(_cutoffRadius);
				global_log->info() << "Setting cell cutoff radius for AutoPas container to " << _cutoffRadius << endl;
#else
				global_log->fatal() << "AutoPas not compiled (use LinkedCells instead, or compile with enabled autopas mode)!" << std::endl;
				Simulation::exit(33);
#endif
			}
			else {
				global_log->error() << "Unknown data structure type: " << datastructuretype << endl;
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
			global_log->error() << "Datastructure section missing" << endl;
			Simulation::exit(1);
		}

		// TODO: move parts to readXML in TemperatureControl?
		/* thermostats */
		if(xmlconfig.changecurrentnode("thermostats")) {
			long numThermostats = 0;
			XMLfile::Query query = xmlconfig.query("thermostat");
			numThermostats = query.card();
			global_log->info() << "Number of thermostats: " << numThermostats << endl;
			if(numThermostats > 1) {
				global_log->info() << "Enabling component wise thermostat" << endl;
				_velocityScalingThermostat.enableComponentwise();
			}
			string oldpath = xmlconfig.getcurrentnodepath();
			XMLfile::Query::const_iterator thermostatIter;
			for( thermostatIter = query.begin(); thermostatIter; thermostatIter++ ) {
				xmlconfig.changecurrentnode( thermostatIter );
				string thermostattype;
				xmlconfig.getNodeValue("@type", thermostattype);
				if(thermostattype == "VelocityScaling") {
					double temperature = _ensemble->T();
					xmlconfig.getNodeValue("temperature", temperature);
					string componentName("global");
					xmlconfig.getNodeValue("@componentId", componentName);
					if(componentName == "global"){
						_domain->setGlobalTemperature(temperature);
						global_log->info() << "Adding global velocity scaling thermostat, T = " << temperature << endl;
					}
					else {
						int componentId = 0;
						componentId = getEnsemble()->getComponent(componentName)->ID();
						int thermostatID = _domain->getThermostat(componentId);
						_domain->setTargetTemperature(thermostatID, temperature);
						global_log->info() << "Adding velocity scaling thermostat for component '" << componentName << "' (ID: " << componentId << "), T = " << temperature << endl;
					}
				}
				else if(thermostattype == "TemperatureControl") {
                    if (nullptr == _temperatureControl) {
                        _temperatureControl = new TemperatureControl();
                        _temperatureControl->readXML(xmlconfig);
                    } else {
                        global_log->error() << "Instance of TemperatureControl allready exist! Programm exit ..."
                                            << endl;
                        Simulation::exit(-1);
                    }
                }
				else
				{
					global_log->warning() << "Unknown thermostat " << thermostattype << endl;
					continue;
				}
			}
			xmlconfig.changecurrentnode(oldpath);
			xmlconfig.changecurrentnode("..");
		}
		else {
			global_log->warning() << "Thermostats section missing." << endl;
		}

		/* long range correction */
		if(xmlconfig.changecurrentnode("longrange") )
		{
			std::string type;
			if( !xmlconfig.getNodeValue("@type", type) )
			{
				global_log->error() << "LongRangeCorrection: Missing type specification. Program exit ..." << endl;
				Simulation::exit(-1);
			}
			if("planar" == type)
			{
				unsigned int nSlabs = 10;
				delete _longRangeCorrection;
				_longRangeCorrection = new Planar(_cutoffRadius, _LJCutoffRadius, _domain, _domainDecomposition, _moleculeContainer, nSlabs, global_simulation);
				_longRangeCorrection->readXML(xmlconfig);
			}
			else if("homogeneous" == type)
			{
				/*
				 * Needs to be initialized later for some reason, done in Simulation::prepare_start()
				 * TODO: perhaps work on this, to make it more robust
				 *
				 *_longRangeCorrection = new Homogeneous(_cutoffRadius, _LJCutoffRadius,_domain,global_simulation);
                 */
			}
			else
			{
				global_log->error() << "LongRangeCorrection: Wrong type. Expected type == homogeneous|planar. Program exit ..." << endl;
                Simulation::exit(-1);
			}
			xmlconfig.changecurrentnode("..");
		}

		xmlconfig.changecurrentnode(".."); /* algorithm section */
	}
	else {
		global_log->error() << "Algorithm section missing." << endl;
	}

	global_log -> info() << "Registering default plugins..." << endl;
    // REGISTERING/ENABLING PLUGINS
	PluginFactory<PluginBase> pluginFactory;
    pluginFactory.registerDefaultPlugins();
	global_log -> info() << "Successfully registered plugins." << endl;

    long numPlugs = 0;
	numPlugs += pluginFactory.enablePlugins(_plugins, xmlconfig, "plugin", _domain);
	numPlugs += pluginFactory.enablePlugins(_plugins, xmlconfig, "output/outputplugin", _domain);
    global_log -> info() << "Number of enabled Plugins: " << numPlugs << endl;


    string oldpath = xmlconfig.getcurrentnodepath();

	if(xmlconfig.changecurrentnode("ensemble/phasespacepoint/file")) {
		global_log->info() << "Reading phase space from file." << endl;
		string pspfiletype;
		xmlconfig.getNodeValue("@type", pspfiletype);
		global_log->info() << "Face space file type: " << pspfiletype << endl;

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
		else {
			global_log->error() << "Unknown phase space file type" << endl;
			Simulation::exit(-1);
		}
	}
	xmlconfig.changecurrentnode(oldpath);

	oldpath = xmlconfig.getcurrentnodepath();
	if(xmlconfig.changecurrentnode("ensemble/phasespacepoint/generator")) {
		string generatorName;
		xmlconfig.getNodeValue("@name", generatorName);
		global_log->info() << "Initializing phase space using generator: " << generatorName << endl;
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
		else {
			global_log->error() << "Unknown generator: " << generatorName << endl;
			Simulation::exit(1);
		}
		_inputReader->readXML(xmlconfig);
	}
	xmlconfig.changecurrentnode(oldpath);

	/** Prepare start options, affecting behavior of method prepare_start() */
	_prepare_start_opt.refreshIDs = false;

	oldpath = xmlconfig.getcurrentnodepath();
	if(xmlconfig.changecurrentnode("options")) {
		unsigned long numOptions = 0;
		XMLfile::Query query = xmlconfig.query("option");
		numOptions = query.card();
		global_log->info() << "Number of prepare start options: " << numOptions << endl;

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
					global_log->info() << "Particle IDs will be refreshed before simulation start." << endl;
				else
					global_log->info() << "Particle IDs will NOT be refreshed before simulation start." << endl;
			}
			else
			{
				global_log->warning() << "Unknown option '" << strOptionName << "'" << endl;
			}
		}
	}
	xmlconfig.changecurrentnode(oldpath);
}


void Simulation::readConfigFile(string filename) {
	string extension(getFileExtension(filename.c_str()));
	global_log->debug() << "Found config filename extension: " << extension << endl;
	if (extension == "xml") {
		initConfigXML(filename);
	}
	else {
		global_log->error() << "Unknown config file extension '" << extension << "'." << endl;
		Simulation::exit(1);
	}
}

void Simulation::initConfigXML(const string& inputfilename) {
	global_log->info() << "Initializing XML config file: " << inputfilename << endl;

	try{
		XMLfileUnits inp(inputfilename);

		global_log->debug() << "Input XML:" << endl << string(inp) << endl;

		if(inp.changecurrentnode("/mardyn") < 0) {
			global_log->error() << "Cound not find root node /mardyn in XML input file." << endl;
			global_log->fatal() << "Not a valid MarDyn XML input file." << endl;
			Simulation::exit(1);
		}

		string version("unknown");
		inp.getNodeValue("@version", version);
		global_log->info() << "MarDyn XML config file version: " << version << endl;

		if(inp.changecurrentnode("simulation")) {
			readXML(inp);
			inp.changecurrentnode("..");
		} // simulation-section
		else {
			global_log->error() << "Simulation section missing" << endl;
			Simulation::exit(1);
		}
	} catch (const std::exception& e) {
		global_log->error() << "Error in XML config. Please check your input file!" << std::endl;
		global_log->error() << "Exception: " << e.what() << std::endl;
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
	timers()->registerTimer("PHASESPACE_CREATION",  vector<string>{"SIMULATION_IO"}, new Timer());
	timers()->setOutputString("PHASESPACE_CREATION", "Phasespace creation took:");
	timers()->getTimer("PHASESPACE_CREATION")->start();

	_inputReader->readPhaseSpace(_moleculeContainer, _domain, _domainDecomposition);
	timers()->getTimer("PHASESPACE_CREATION")->stop();

	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	_domain->updateglobalNumMolecules(_moleculeContainer, _domainDecomposition);
	unsigned long globalNumMolecules = _domain->getglobalNumMolecules();
	double rho_global = globalNumMolecules/ _ensemble->V();
	global_log->info() << "Setting domain class parameters: N_global: " << globalNumMolecules << ", rho_global: " << rho_global << ", T_global: " << _ensemble->T() << endl;
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
	global_log->info() << "Initializing simulation" << endl;

	global_log->info() << "Initialising cell processor" << endl;
	if (!_legacyCellProcessor) {
#ifndef ENABLE_REDUCED_MEMORY_MODE
		global_log->info() << "Using vectorized cell processor." << endl;
		_cellProcessor = new VectorizedCellProcessor( *_domain, _cutoffRadius, _LJCutoffRadius);
#else
		global_log->info() << "Using reduced memory mode (RMM) cell processor." << endl;
		_cellProcessor = new VCP1CLJRMM( *_domain, _cutoffRadius, _LJCutoffRadius);
#endif
	} else {
		global_log->info() << "Using legacy cell processor." << endl;
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
	if(dynamic_cast<KDDecomposition*>(_domainDecomposition) != nullptr){
		static_cast<KDDecomposition*>(_domainDecomposition)->fillTimeVecs(&_cellProcessor);
	}
#endif

	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();
	global_log->info() << "Updating domain decomposition" << endl;

	if(getMemoryProfiler()) {
		getMemoryProfiler()->doOutput("without halo copies");
	}

	// temporary addition until MPI communication is parallelized with OpenMP
	// we don't actually need the mpiOMPCommunicationTimer here -> deactivate it..
	global_simulation->timers()->deactivateTimer("SIMULATION_MPI_OMP_COMMUNICATION");
	updateParticleContainerAndDecomposition(1.0);
	global_simulation->timers()->activateTimer("SIMULATION_MPI_OMP_COMMUNICATION");

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

	global_log->info() << "Performing initial force calculation" << endl;
	global_simulation->timers()->start("SIMULATION_FORCE_CALCULATION");
	_moleculeContainer->traverseCells(*_cellProcessor);

	global_simulation->timers()->stop("SIMULATION_FORCE_CALCULATION");
	global_log->info() << "Performing initial FLOP count (if necessary)" << endl;

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
		global_log->info() << "Performing initial FMM force calculation" << endl;
		_FMM->computeElectrostatics(_moleculeContainer);
	}

	/** Init TemperatureControl beta_trans, beta_rot log-files, register as observer if plugin DistControl is in use. */
	if(nullptr != _temperatureControl)
		_temperatureControl->prepare_start();  // Has to be called before plugin initialization (see below): plugin->init(...)

	// initializing plugins and starting plugin timers
	for (auto& plugin : _plugins) {
		global_log->info() << "Initializing plugin " << plugin->getPluginName() << endl;
		plugin->init(_moleculeContainer, _domainDecomposition, _domain);
		string timer_name = plugin->getPluginName();
		// TODO: real timer
		global_simulation->timers()->registerTimer(timer_name, vector<string>{"SIMULATION_PER_STEP_IO"}, new Timer());
		string timer_plugin_string = string("Plugin ") + timer_name + string(" took:");
		global_simulation->timers()->setOutputString(timer_name, timer_plugin_string);
	}

	//afterForces Plugin Call
	global_log->debug() << "[AFTER FORCES] Performing AfterForces plugin call"
						<< endl;
	for (auto plugin : _plugins) {
		global_log->debug() << "[AFTER FORCES] Plugin: "
							<< plugin->getPluginName() << endl;
		plugin->afterForces(_moleculeContainer, _domainDecomposition, _simstep);
	}

#ifndef MARDYN_AUTOPAS
	// clear halo
	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();
#endif

	if (_longRangeCorrection == nullptr){
		_longRangeCorrection = new Homogeneous(_cutoffRadius, _LJCutoffRadius,_domain,this);
	}

	_longRangeCorrection->calculateLongRange();
	// here we have to call calcFM() manually, otherwise force and moment are not
	// updated inside the molecule (actually this is done in upd_postF)
	// integrator->eventForcesCalculated should not be called, since otherwise the velocities would already be updated.
	//updateForces();

	global_log->info() << "Calculating global values" << endl;
	_domain->calculateThermostatDirectedVelocity(_moleculeContainer);

	_domain->calculateVelocitySums(_moleculeContainer);

	_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer,
			true, 1.0);
	global_log->debug() << "Calculating global values finished." << endl;

	_ensemble->prepare_start();

	_simstep = _initSimulation = (unsigned long) round(_simulationTime / _integrator->getTimestepLength() );
	global_log->info() << "Set initial time step to start from to " << _initSimulation << endl;
	global_log->info() << "System initialised with " << _domain->getglobalNumMolecules() << " molecules." << endl;

	/** refresh particle IDs */
	if(_prepare_start_opt.refreshIDs)
		this->refreshParticleIDs();
}

void Simulation::simulate() {
	global_log->info() << "Started simulation" << endl;

	_ensemble->updateGlobalVariable(_moleculeContainer, NUM_PARTICLES);
	global_log->debug() << "Number of particles in the Ensemble: " << _ensemble->N() << endl;
// 	ensemble.updateGlobalVariable(ENERGY);
// 	global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E() << endl;
// 	ensemble.updateGlobalVariable(TEMPERATURE);
// 	global_log->debug() << "Temperature of the Ensemble: " << ensemble.T() << endl;*/

	/***************************************************************************/
	/* BEGIN MAIN LOOP                                                         */
	/***************************************************************************/

	global_simulation->timers()->setOutputString("SIMULATION_LOOP", "Computation in main loop took:");
	global_simulation->timers()->setOutputString("SIMULATION_DECOMPOSITION", "Decomposition took:");
	global_simulation->timers()->setOutputString("SIMULATION_COMPUTATION", "Computation took:");
	global_simulation->timers()->setOutputString("SIMULATION_PER_STEP_IO", "IO in main loop took:");
	global_simulation->timers()->setOutputString("SIMULATION_FORCE_CALCULATION", "Force calculation took:");
	global_simulation->timers()->setOutputString("SIMULATION_MPI_OMP_COMMUNICATION", "Communication took:");
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
		global_log->debug() << "timestep: " << getSimulationStep() << endl;
		global_log->debug() << "simulation time: " << getSimulationTime() << endl;
		global_simulation->timers()->incrementTimerTimestepCounter();

		computationTimer->start();

        // beforeEventNewTimestep Plugin Call
        global_log -> debug() << "[BEFORE EVENT NEW TIMESTEP] Performing beforeEventNewTimestep plugin call" << endl;
        for (auto plugin : _plugins) {
            global_log -> debug() << "[BEFORE EVENT NEW TIMESTEP] Plugin: " << plugin->getPluginName() << endl;
            plugin->beforeEventNewTimestep(_moleculeContainer, _domainDecomposition, _simstep);
        }

        _ensemble->beforeEventNewTimestep(_moleculeContainer, _domainDecomposition, _simstep);

		_integrator->eventNewTimestep(_moleculeContainer, _domain);

        // beforeForces Plugin Call
        global_log -> debug() << "[BEFORE FORCES] Performing BeforeForces plugin call" << endl;
        for (auto plugin : _plugins) {
            global_log -> debug() << "[BEFORE FORCES] Plugin: " << plugin->getPluginName() << endl;
            plugin->beforeForces(_moleculeContainer, _domainDecomposition, _simstep);
        }

		computationTimer->stop();



#if defined(ENABLE_MPI) && defined(ENABLE_OVERLAPPING)
		bool overlapCommComp = true;
#else
		bool overlapCommComp = false;
#endif


		if (overlapCommComp) {
			double currentTime = _timerForLoad->get_etime();
			performOverlappingDecompositionAndCellTraversalStep(currentTime - previousTimeForLoad);
			previousTimeForLoad = currentTime;
		}
		else {
			decompositionTimer->start();
			// ensure that all Particles are in the right cells and exchange Particles
			global_log->debug() << "Updating container and decomposition" << endl;
			double currentTime = _timerForLoad->get_etime();
#if defined(ENABLE_MPI)
			updateParticleContainerAndDecomposition(currentTime - previousTimeForLoad); //  - *(_domainDecomposition->getProcessTimerPointer())
			// FIXME: Time from new Timer needs to be subtracted incl. Reset

#else
			updateParticleContainerAndDecomposition(currentTime - previousTimeForLoad);
#endif
			previousTimeForLoad = currentTime;

			decompositionTimer->stop();

			double startEtime = computationTimer->get_etime();
			// Force calculation and other pair interaction related computations
			global_log->debug() << "Traversing pairs" << endl;
			computationTimer->start();
			forceCalculationTimer->start();

			_moleculeContainer->traverseCells(*_cellProcessor);

			// siteWiseForces Plugin Call
			global_log -> debug() << "[SITEWISE FORCES] Performing siteWiseForces plugin call" << endl;
			for (auto plugin : _plugins) {
				global_log -> debug() << "[SITEWISE FORCES] Plugin: " << plugin->getPluginName() << endl;
				plugin->siteWiseForces(_moleculeContainer, _domainDecomposition, _simstep);
			}

			// Update forces in molecules so they can be exchanged
			updateForces();

			forceCalculationTimer->stop();
			computationTimer->stop();

			decompositionTimer->start();
			// Exchange forces if it's required by the cell container.
			if(_moleculeContainer->requiresForceExchange()){
				global_log->debug() << "Exchanging Forces" << std::endl;
				_domainDecomposition->exchangeForces(_moleculeContainer, _domain);
			}
			decompositionTimer->stop();
			_loopCompTime += computationTimer->get_etime() - startEtime;
			_loopCompTimeSteps ++;
		}
		computationTimer->start();


		if (_FMM != nullptr) {
			global_log->debug() << "Performing FMM calculation" << endl;
			_FMM->computeElectrostatics(_moleculeContainer);
		}

		//afterForces Plugin Call
		global_log -> debug() << "[AFTER FORCES] Performing AfterForces plugin call" << endl;
		for (auto plugin : _plugins) {
			global_log -> debug() << "[AFTER FORCES] Plugin: " << plugin->getPluginName() << endl;
			plugin->afterForces(_moleculeContainer, _domainDecomposition, _simstep);
		}

		_ensemble->afterForces(_moleculeContainer, _domainDecomposition, _cellProcessor, _simstep);

		// TODO: test deletions and insertions
		global_log->debug() << "Deleting outer particles / clearing halo." << endl;
#ifndef MARDYN_AUTOPAS
		_moleculeContainer->deleteOuterParticles();
#endif

		if (!(_simstep % _collectThermostatDirectedVelocity))
			_domain->calculateThermostatDirectedVelocity(_moleculeContainer);
		_longRangeCorrection->calculateLongRange();
		_longRangeCorrection->writeProfiles(_domainDecomposition, _domain, _simstep);

		_ensemble->beforeThermostat(_simstep, _initStatistics);

		global_log->debug() << "Inform the integrator (forces calculated)" << endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);

		// calculate the global macroscopic values from the local values
		global_log->debug() << "Calculate macroscopic values" << endl;
		_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer,
				(!(_simstep % _collectThermostatDirectedVelocity)), Tfactor(_simstep));

		// scale velocity and angular momentum
        // TODO: integrate into Temperature Control
		if ( !_domain->NVE() && _temperatureControl == nullptr) {
			if (_thermostatType ==VELSCALE_THERMOSTAT) {
				global_log->debug() << "Velocity scaling" << endl;
				if (_domain->severalThermostats()) {
					_velocityScalingThermostat.enableComponentwise();
					for(unsigned int cid = 0; cid < global_simulation->getEnsemble()->getComponents()->size(); cid++) {
						int thermostatId = _domain->getThermostat(cid);
						_velocityScalingThermostat.setBetaTrans(thermostatId, _domain->getGlobalBetaTrans(thermostatId));
						_velocityScalingThermostat.setBetaRot(thermostatId, _domain->getGlobalBetaRot(thermostatId));
						global_log->debug() << "Thermostat for CID: " << cid << " thermID: " << thermostatId
								<< " B_trans: " << _velocityScalingThermostat.getBetaTrans(thermostatId)
								<< " B_rot: " << _velocityScalingThermostat.getBetaRot(thermostatId) << endl;
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
		// global_log->debug() << "Number of particles in the Ensemble: " << _ensemble->N() << endl;
		/*
		ensemble.updateGlobalVariable(ENERGY);
		global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E() << endl;
		ensemble.updateGlobalVariable(TEMPERATURE);
		global_log->debug() << "Temperature of the Ensemble: " << ensemble.T() << endl;
		*/
		/* END PHYSICAL SECTION */


		computationTimer->stop();
		perStepIoTimer->start();

		// CALL ALL PLUGIN ENDSTEP METHODS
		pluginEndStepCall(_simstep);

		if( (_forced_checkpoint_time > 0) && (loopTimer->get_etime() >= _forced_checkpoint_time) ) {
			/* force checkpoint for specified time */
			string cpfile(_outputPrefix + ".timed.restart.dat");
			global_log->info() << "Writing timed, forced checkpoint to file '" << cpfile << "'" << endl;
			_domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition, _simulationTime);
			_forced_checkpoint_time = -1; /* disable for further timesteps */
		}
		perStepIoTimer->stop();
	}
	loopTimer->stop();
	/***************************************************************************/
	/* END MAIN LOOP                                                           */
	/***************************************************************************/


	global_simulation->timers()->registerTimer("SIMULATION_FINAL_IO", vector<string>{"SIMULATION_IO"}, new Timer());
	global_simulation->timers()->setOutputString("SIMULATION_FINAL_IO", "Final IO took:");
	global_simulation->timers()->getTimer("SIMULATION_FINAL_IO")->start();
    if( _finalCheckpoint ) {
        /* write final checkpoint */
        string cpfile(_outputPrefix + ".restart.dat");
        global_log->info() << "Writing final checkpoint to file '" << cpfile << "'" << endl;
        _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition, _simulationTime, false);
    }

	global_log->info() << "Finish plugins" << endl;
	for (auto plugin : _plugins) {
		plugin->finish(_moleculeContainer, _domainDecomposition, _domain);
	}
	global_simulation->timers()->getTimer("SIMULATION_FINAL_IO")->stop();

	global_log->info() << "Timing information:" << endl;
	global_simulation->timers()->printTimers();
	global_simulation->timers()->resetTimers();
	if(getMemoryProfiler()) {
		getMemoryProfiler()->doOutput();
	}
	global_log->info() << endl;

#ifdef WITH_PAPI
	global_log->info() << "PAPI counter values for loop timer:"  << endl;
	for(int i = 0; i < loopTimer->get_papi_num_counters(); i++) {
		global_log->info() << "  " << papi_event_list[i] << ": " << loopTimer->get_global_papi_counter(i) << endl;
	}
#endif /* WITH_PAPI */
}

void Simulation::pluginEndStepCall(unsigned long simstep) {

	std::list<PluginBase*>::iterator pluginIter;
	for (pluginIter = _plugins.begin(); pluginIter != _plugins.end(); pluginIter++) {
		PluginBase* plugin = (*pluginIter);
		global_log->debug() << "Plugin end of step: " << plugin->getPluginName() << endl;
		global_simulation->timers()->start(plugin->getPluginName());
		plugin->endStep(_moleculeContainer, _domainDecomposition, _domain, simstep);
		global_simulation->timers()->stop(plugin->getPluginName());
	}


	if (_domain->thermostatWarning())
		global_log->warning() << "Thermostat!" << endl;
	/* TODO: thermostat */
	global_log->info() << "Simstep = " << simstep << "\tT = "
					   << _domain->getGlobalCurrentTemperature() << "\tU_pot = "
					   << _domain->getGlobalUpot() << "\tp = "
					   << _domain->getGlobalPressure() << endl;
	using std::isnan;
	if (isnan(_domain->getGlobalCurrentTemperature()) || isnan(_domain->getGlobalUpot()) || isnan(_domain->getGlobalPressure())) {
		global_log->error() << "NaN detected, exiting." << std::endl;
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

void Simulation::updateParticleContainerAndDecomposition(double lastTraversalTime) {
	// The particles have moved, so the neighborhood relations have
	// changed and have to be adjusted
	_moleculeContainer->update();
	//_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain);
	bool forceRebalancing = false;
	global_simulation->timers()->start("SIMULATION_MPI_OMP_COMMUNICATION");
	_domainDecomposition->balanceAndExchange(lastTraversalTime, forceRebalancing, _moleculeContainer, _domain);
	global_simulation->timers()->stop("SIMULATION_MPI_OMP_COMMUNICATION");
	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateMoleculeCaches();
}

void Simulation::performOverlappingDecompositionAndCellTraversalStep(double etime) {
	bool forceRebalancing = false;

	#ifdef ENABLE_MPI
		#ifdef ENABLE_OVERLAPPING
			NonBlockingMPIHandlerBase* nonBlockingMPIHandler =
					new NonBlockingMPIMultiStepHandler(static_cast<DomainDecompMPIBase*>(_domainDecomposition), _moleculeContainer, _domain, _cellProcessor);
		#else
			NonBlockingMPIHandlerBase* nonBlockingMPIHandler =
					new NonBlockingMPIHandlerBase(static_cast<DomainDecompMPIBase*>(_domainDecomposition), _moleculeContainer, _domain, _cellProcessor);
		#endif

		nonBlockingMPIHandler->performOverlappingTasks(forceRebalancing, etime);
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
	global_log->info() << "Creating alibi domain decomposition ... " << endl;
	_domainDecomposition = new DomainDecompBase();
#else
	global_log->info() << "Creating standard domain decomposition ... " << endl;
	_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
#endif

	_outputPrefix.append(gettimestring());

	global_log->info() << "Creating domain ..." << endl;
	_domain = new Domain(ownrank);
	global_log->info() << "Creating ParticlePairs2PotForceAdapter ..." << endl;
	_particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);

	global_log->info() << "Initialization done" << endl;
}

bool Simulation::keepRunning() {

	// Simstep Criterion
	if (_simstep >= _numberOfTimesteps){
		global_log->info() << "Maximum Simstep reached: " << _simstep << std::endl;
		return false;
	}
	// WallTime Criterion, elapsed time since Simulation constructor
	else if(_wallTimeEnabled && _timeFromStart.get_etime_running() > _maxWallTime){
		global_log->info() << "Maximum Walltime reached (s): " << _maxWallTime << std::endl;
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
	uint64_t prevMaxID = 0;  // max ID of previous process
	int ownRank = _domainDecomposition->getRank();

#ifdef ENABLE_MPI
	int numProcs = _domainDecomposition->getNumProcs();
	if (ownRank != 0) {
		MPI_Recv(&prevMaxID, 1, MPI_UNSIGNED_LONG, (ownRank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifndef NDEBUG
		cout << "Process " << ownRank << " received maxID=" << prevMaxID << " from process " << (ownRank-1) << "." << endl;
#endif
	}
#endif

	uint64_t tmpID = prevMaxID;

#ifdef ENABLE_MPI
	if(ownRank < (numProcs-1) ) {
		prevMaxID += _moleculeContainer->getNumberOfParticles();
		MPI_Send(&prevMaxID, 1, MPI_UNSIGNED_LONG, (ownRank+1), 0, MPI_COMM_WORLD);
	}
#endif

#ifndef NDEBUG
	cout << "["<<ownRank<<"]tmpID=" << tmpID << endl;
#endif
	for (auto pit = _moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		pit->setid(++tmpID);
	}
}
