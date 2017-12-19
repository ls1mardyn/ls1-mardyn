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
#include "parallel/DomainDecompBase.h"
#include "parallel/NonBlockingMPIHandlerBase.h"
#include "parallel/NonBlockingMPIMultiStepHandler.h"
#include "molecules/Molecule.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/adapter/VCP1CLJRMM.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "integrators/LeapfrogRMM.h"
#include "molecules/Wall.h"
#include "molecules/Mirror.h"

#include "io/OutputBase.h"
#include "io/OutputPluginFactory.h"

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

#include "ensemble/GrandCanonical.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/PressureGradient.h"
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

#include "particleContainer/adapter/VectorizationTuner.h"

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
	_thermostatType(VELSCALE_THERMOSTAT),
	_nuAndersen(0.0),
	_numberOfTimesteps(1),
	_simstep(0),
	_initSimulation(0),
	_initCanonical(_initSimulation + 1),  // simulate canonical (with Tfactor == 1.) from begin on!
	_initGrandCanonical(10000000),
	_initStatistics(20000),
	_ensemble(nullptr),
	_pressureGradient(nullptr),
	_moleculeContainer(nullptr),
	_particlePairsHandler(nullptr),
	_cellProcessor(nullptr),
	_domainDecomposition(nullptr),
	_integrator(nullptr),
	_domain(nullptr),
	_inputReader(nullptr),
	_outputPrefix("mardyn"),
	_doAlignCentre(false),
	_componentSpecificAlignment(false),
	_alignmentInterval(25),
	_alignmentCorrection(0.0),
	_applyWallFun_LJ_9_3(false),
	_applyWallFun_LJ_10_4(false),
	_applyMirror(false),
	_wall(nullptr),
	_mirror(nullptr),
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
	_loopCompTimeSteps(0.0),
	_nFmaxOpt(CFMAXOPT_NO_CHECK),
	_nFmaxID(0),
	_dFmaxInit(0.0),
	_dFmaxThreshold(0.0),
	_virialRequired(false)
{
	_ensemble = new CanonicalEnsemble();
	initialize();
}

Simulation::~Simulation() {
	delete _ensemble;
	_ensemble = nullptr;
	delete _pressureGradient;
	_pressureGradient = nullptr;
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
	delete _wall;
	_wall = nullptr;
	delete _mirror;
	_mirror = nullptr;
	delete _longRangeCorrection;
	_longRangeCorrection = nullptr;
	delete _temperatureControl;
	_temperatureControl = nullptr;
	delete _FMM;
	_FMM = nullptr;
	/* destruct output plugins and remove them from the output plugin list */
	_outputPlugins.remove_if([](OutputBase *pluginPtr) {delete pluginPtr; return true;} );
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
			if(_ensemble != NULL) delete _ensemble;
			_ensemble = new CanonicalEnsemble();
		} else if (ensembletype == "muVT") {
			global_log->error() << "muVT ensemble not completely implemented via XML input." << endl;
			Simulation::exit(1);
			// _ensemble = new GrandCanonicalEnsemble();
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

		/// @todo Move the phasespace check into an output plugin executed once at the start!
		/* check phasespacepoint */
		bool bInputOk = true;
		_nFmaxOpt = CFMAXOPT_NO_CHECK;
		std::string strType = "unknown";
		std::string strOption = "unknown";
		_dFmaxInit = 0.;
		_dFmaxThreshold = 0.;
		bInputOk = bInputOk && xmlconfig.getNodeValue("phasespacepoint/check/@type", strType);
		bInputOk = bInputOk && xmlconfig.getNodeValue("phasespacepoint/check/@option", strOption);
		if(true == bInputOk)
		{
			if("show-only" == strOption)
				_nFmaxOpt = CFMAXOPT_SHOW_ONLY;
			else if("check-greater" == strOption)
			{
				_nFmaxOpt = CFMAXOPT_CHECK_GREATER;
				bInputOk = bInputOk && xmlconfig.getNodeValue("phasespacepoint/check/Fmax", _dFmaxThreshold);
				global_log->info() << "Checking if initial max. force (Fmax) is less than:"
						" " << _dFmaxThreshold << endl;
			}
			else
				global_log->error() << "XML-I/O: Wrong option in attribute: ensemble/phasespacepoint/check/@option" << endl;

			if("Fmax" == strType && CFMAXOPT_NO_CHECK != _nFmaxOpt)
				global_log->info() << "Checking max. initial force (Fmax)." << endl;
		}
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
			}
			else {
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
			xmlconfig.changecurrentnode("..");
		}
		else {
		#ifdef ENABLE_MPI
			global_log->error() << "Parallelisation section missing." << endl;
			Simulation::exit(1);
		#else /* serial */
			//_domainDecomposition = new DomainDecompBase(); // already set in initialize()
		#endif
		}

		/* datastructure */
		if(xmlconfig.changecurrentnode("datastructure")) {
			string datastructuretype;
			xmlconfig.getNodeValue("@type", datastructuretype);
			global_log->info() << "Datastructure type: " << datastructuretype << endl;
			if(datastructuretype == "LinkedCells") {
							_moleculeContainer = new LinkedCells();
							/** @todo Review if we need to know the max cutoff radius usable with any datastructure. */
							global_log->info() << "Setting cell cutoff radius for linked cell datastructure to " << _cutoffRadius << endl;
							LinkedCells *lc = static_cast<LinkedCells*>(_moleculeContainer);
							lc->setCutoff(_cutoffRadius);
			}else if(datastructuretype == "AdaptiveSubCells") {
				global_log->warning() << "AdaptiveSubCells no longer supported." << std::endl;
				Simulation::exit(-1);
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
				else if(thermostattype == "TemperatureControl")
				{
					if(NULL == _temperatureControl)
					{
						_temperatureControl = new TemperatureControl();
						_temperatureControl->readXML(xmlconfig);
					}
					else
					{
						global_log->error() << "Instance of TemperatureControl allready exist! Programm exit ..." << endl;
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
				if(NULL != _longRangeCorrection)
					delete _longRangeCorrection;
				_longRangeCorrection = new Planar(_cutoffRadius, _LJCutoffRadius, _domain, _domainDecomposition, _moleculeContainer, nSlabs, global_simulation);
				_longRangeCorrection->readXML(xmlconfig);
			}
			else if("homogeneous" == type)
			{
				/*
				 * Needs to be initialized later for some reason, done in Simulation::prepare_start()
				 * TODO: perhabs work on this, to make it more robust 
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

		/* features */
		if(xmlconfig.changecurrentnode("features")) {
			int numFeatures = 0;
			XMLfile::Query query = xmlconfig.query("feature");
			numFeatures = query.card();
			global_log->info() << "Number of features: " << numFeatures << endl;
			if(numFeatures < 1) {
				global_log->warning() << "No feature specified." << endl;
			}
			string oldpath = xmlconfig.getcurrentnodepath();
			XMLfile::Query::const_iterator featureIter;
			for( featureIter = query.begin(); featureIter; featureIter++ ) {
				xmlconfig.changecurrentnode( featureIter );
				string featureName;
				xmlconfig.getNodeValue("@name", featureName);
				if("Wallfun" == featureName) {
					_wall = new Wall();
					_wall->readXML(xmlconfig);
					_applyWallFun_LJ_9_3 = true;
				}
				else
				{
					global_log->warning() << "Unknown feature " << featureName << endl;
					continue;
				}
			}
			xmlconfig.changecurrentnode(oldpath);
			xmlconfig.changecurrentnode("..");
		}

		xmlconfig.changecurrentnode(".."); /* algorithm section */
	}
	else {
		global_log->error() << "Algorithm section missing." << endl;
	}

	/* output */
	long numOutputPlugins = 0;
	XMLfile::Query query = xmlconfig.query("output/outputplugin");
	numOutputPlugins = query.card();
	global_log->info() << "Number of output plugins: " << numOutputPlugins << endl;
	if(numOutputPlugins < 1) {
		global_log->warning() << "No output plugins specified." << endl;
	}

	string oldpath = xmlconfig.getcurrentnodepath();
	OutputPluginFactory outputPluginFactory;
	for( auto outputPluginIter = query.begin(); outputPluginIter; ++outputPluginIter ) {
		xmlconfig.changecurrentnode( outputPluginIter );
		string pluginname("");
		xmlconfig.getNodeValue("@name", pluginname);
		bool enabled = true;
		xmlconfig.getNodeValue("@enabled", enabled);
		if(not enabled) {
			global_log->debug() << "Skipping disabled output plugin: " << pluginname << endl;
			continue;
		}
		global_log->info() << "Enabling output plugin: " << pluginname << endl;


		OutputBase *outputPlugin = outputPluginFactory.create(pluginname);
		if(outputPlugin == nullptr) {
			global_log->warning() << "Could not create output plugin using factory: " << pluginname << endl;
		}


		if(pluginname == "MmpldWriter") {
			/** @todo this should be handled in the MMPLD Writer readXML() */
			std::string sphere_representation = "simple";
			xmlconfig.getNodeValue("@type", sphere_representation);
			delete outputPlugin;
			if("simple" == sphere_representation) {
				outputPlugin = new MmpldWriterSimpleSphere();
			} else if("multi" == sphere_representation) {
				outputPlugin = new MmpldWriterMultiSphere ();
			} else {
				global_log->error() << "[MMPLD Writer] Unknown sphere representation type: " << sphere_representation << endl;
				Simulation::exit(-1);
			}
		}
		/* temporary */
		else if(pluginname == "VectorizationTuner") {
			outputPlugin = new VectorizationTuner(_cutoffRadius, _LJCutoffRadius, &_cellProcessor);
		}
		else if(pluginname == "DomainProfiles") {
			outputPlugin = outputPluginFactory.create("DensityProfileWriter");
			/** @todo This is ugly. Needed as the vectorized code does not support the virial ...*/
			if(xmlconfig.getNodeValue_int("/mardyn/simulation/output/outputplugin[@name='DomainProfiles']/options/option[@keyword='profileVirial']", 0) > 0) {
				_virialRequired = true;
			}
			_domain->readXML(xmlconfig);
		}

		if(nullptr != outputPlugin) {
			outputPlugin->readXML(xmlconfig);
			_outputPlugins.push_back(outputPlugin);
		} else {
			global_log->warning() << "Unknown plugin " << pluginname << endl;
		}
	}
	xmlconfig.changecurrentnode(oldpath);

	oldpath = xmlconfig.getcurrentnodepath();
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
		temp->initCommunicationPartners(_cutoffRadius, _domain);
	}
#endif

	// read particle data (or generate particles, if a generator is chosen)
	timers()->registerTimer("PHASESPACE_CREATION",  vector<string>{"SIMULATION_IO"}, new Timer());
	timers()->setOutputString("PHASESPACE_CREATION", "Phasespace creation took:");
	timers()->getTimer("PHASESPACE_CREATION")->start();
	unsigned long globalNumMolecules = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain, _domainDecomposition);
	timers()->getTimer("PHASESPACE_CREATION")->stop();

	double rho_global = globalNumMolecules/ _ensemble->V();
	global_log->info() << "Setting domain class parameters: N_global: " << globalNumMolecules << ", rho_global: " << rho_global << ", T_global: " << _ensemble->T() << endl;
	_domain->setglobalNumMolecules(globalNumMolecules);
	_domain->setGlobalTemperature(_ensemble->T());
	_domain->setglobalRho(rho_global);

	_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
	//domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif
	unsigned idi = _lmu.size();
	unsigned j = 0;
	std::list<ChemicalPotential>::iterator cpit;
	for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
		cpit->setIncrement(idi);
		double tmp_molecularMass = global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->m();
		cpit->setSystem(_domain->getGlobalLength(0),
				_domain->getGlobalLength(1), _domain->getGlobalLength(2),
				tmp_molecularMass);
		cpit->setGlobalN(global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->getNumMolecules());
		cpit->setNextID(j + (int) (1.001 * (256 + globalNumMolecules)));

		cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0),
				_moleculeContainer->getBoundingBoxMax(0),
				_moleculeContainer->getBoundingBoxMin(1),
				_moleculeContainer->getBoundingBoxMax(1),
				_moleculeContainer->getBoundingBoxMin(2),
				_moleculeContainer->getBoundingBoxMax(2));
		/* TODO: thermostat */
		double Tcur = _domain->getCurrentTemperature(0);
		/* FIXME: target temperature from thermostat ID 0 or 1?  */
		double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;
		cpit->submitTemperature(Tcur);
		if (h != 0.0)
			cpit->setPlanckConstant(h);
		j++;
	}
}

void Simulation::updateForces() {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const ParticleIterator begin = _moleculeContainer->iteratorBegin();
		const ParticleIterator end = _moleculeContainer->iteratorEnd();

		if(_simstep==0 && (CFMAXOPT_SHOW_ONLY == _nFmaxOpt || CFMAXOPT_CHECK_GREATER == _nFmaxOpt)) {

			uint64_t id;
			uint32_t cid;
			double r[3];
			double F[3];
			double Fabs=0.;
			double FabsSquared=0.;
			double FmaxInitSquared=_dFmaxInit*_dFmaxInit;
			double FmaxThresholdSquared = _dFmaxThreshold*_dFmaxThreshold;

			for (ParticleIterator i = begin; i != end; ++i){
				i->calcFM();

				id=i->id();
				cid=i->componentid()+1;
				for(uint8_t d=0; d<3; ++d)
				{
					r[d]=i->r(d);
					F[d]=i->F(d);
				}
				FabsSquared=(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
				if(FabsSquared > FmaxInitSquared)
				{
					Fabs=sqrt(FabsSquared);
					_dFmaxInit=Fabs;
					FmaxInitSquared=_dFmaxInit*_dFmaxInit;
					_nFmaxID=id;
				}
				if(CFMAXOPT_CHECK_GREATER == _nFmaxOpt && FabsSquared > FmaxThresholdSquared) {
					Fabs=sqrt(FabsSquared);
					global_log->warning()<<"Fabs="<<Fabs<<" > "<<_dFmaxThreshold<<" (threshold) for molecule: id="<<id<<", cid="<<cid<<", "
						"x,y,z="<<r[0]<<", "<<r[1]<<", "<<r[2]<<endl;
				}
			}
			global_log->info()<<"Max. initial force is found for molecule: id="<<_nFmaxID<<", Fmax="<<_dFmaxInit<<endl;
		}
		else {
			for (ParticleIterator i = begin; i != end; ++i){
				i->calcFM();
			}
		}
	} // end pragma omp parallel
}

void Simulation::prepare_start() {
	global_log->info() << "Initializing simulation" << endl;

	global_log->info() << "Initialising cell processor" << endl;
#if ENABLE_VECTORIZED_CODE
#ifndef ENABLE_REDUCED_MEMORY_MODE
	if(_virialRequired) {
		global_log->warning() << "Using legacy cell processor. (The vectorized code does not support the virial tensor and the localized virial profile.)" << endl;
		_cellProcessor = new LegacyCellProcessor(_cutoffRadius, _LJCutoffRadius, _particlePairsHandler);
	} else if (nullptr != getOutputPlugin("RDF")) {
		global_log->warning() << "Using legacy cell processor. (The vectorized code does not support rdf sampling.)" << endl;
		_cellProcessor = new LegacyCellProcessor(_cutoffRadius, _LJCutoffRadius, _particlePairsHandler);
	} else {
		global_log->info() << "Using vectorized cell processor." << endl;
		_cellProcessor = new VectorizedCellProcessor( *_domain, _cutoffRadius, _LJCutoffRadius);
	}
#else
	global_log->info() << "Using reduced memory mode (RMM) cell processor." << endl;
	_cellProcessor = new VCP1CLJRMM( *_domain, _cutoffRadius, _LJCutoffRadius);
#endif // ENABLE_REDUCED_MEMORY_MODE
#else
	global_log->info() << "Using legacy cell processor." << endl;
	_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _particlePairsHandler);
#endif // ENABLE_VECTORIZED_CODE

	if (_FMM != NULL) {

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
		_FMM->init(globalLength, bBoxMin, bBoxMax,
				dynamic_cast<LinkedCells*>(_moleculeContainer)->cellLength(), _moleculeContainer);

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
	measureFLOPRate(_moleculeContainer, 0);
	
	/*
	// Update forces in molecules so they can be exchanged - future
	updateForces(); 
	
	// Exchange forces if it's required by the cell container.
	if(_moleculeContainer->requiresForceExchange()){
		_domainDecomposition->exchangeForces(_moleculeContainer, _domain);
	}
	*/


#ifdef ENABLE_REDUCED_MEMORY_MODE
	// now set vcp1clj_wr_cellProcessor::_dtInvm back.
	vcp1clj_wr_cellProcessor->setDtInvm(dt_inv_m);
#endif /* ENABLE_REDUCED_MEMORY_MODE */

	_loopCompTime = global_simulation->timers()->getTime("SIMULATION_FORCE_CALCULATION");
	global_simulation->timers()->reset("SIMULATION_FORCE_CALCULATION");
	++_loopCompTimeSteps;


	if (_FMM != NULL) {
		global_log->info() << "Performing initial FMM force calculation" << endl;
		_FMM->computeElectrostatics(_moleculeContainer);
	}

    // clear halo
    global_log->info() << "Clearing halos" << endl;
    _moleculeContainer->deleteOuterParticles();

    if (_longRangeCorrection == NULL){
        _longRangeCorrection = new Homogeneous(_cutoffRadius, _LJCutoffRadius,_domain,this);
    }


    _longRangeCorrection->calculateLongRange();
	// here we have to call calcFM() manually, otherwise force and moment are not
	// updated inside the molecule (actually this is done in upd_postF)
	// integrator->eventForcesCalculated should not be called, since otherwise the velocities would already be updated.
	updateForces();

	if (_pressureGradient->isAcceleratingUniformly()) {
		global_log->info() << "Initialising uniform acceleration." << endl;
		unsigned long uCAT = _pressureGradient->getUCAT();
		global_log->info() << "uCAT: " << uCAT << " steps." << endl;
		_pressureGradient->determineAdditionalAcceleration(
				_domainDecomposition, _moleculeContainer, uCAT
						* _integrator->getTimestepLength());
		global_log->info() << "Uniform acceleration initialised." << endl;
	}

	global_log->info() << "Calculating global values" << endl;
	_domain->calculateThermostatDirectedVelocity(_moleculeContainer);

	_domain->calculateVelocitySums(_moleculeContainer);

	_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer,
			true, 1.0);
	global_log->debug() << "Calculating global values finished." << endl;

	if (_lmu.size() + _mcav.size() > 0) {
		/* TODO: thermostat */
		double Tcur = _domain->getGlobalCurrentTemperature();
		/* FIXME: target temperature from thermostat ID 0 or 1? */
		double
				Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1)
								: _domain->getTargetTemperature(0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;

		list<ChemicalPotential>::iterator cpit;
		if (h == 0.0)
			h = sqrt(6.2831853 * Ttar);
		for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
			cpit->submitTemperature(Tcur);
			cpit->setPlanckConstant(h);
		}
                map<unsigned, CavityEnsemble>::iterator ceit;
		for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {
		   ceit->second.submitTemperature(Tcur);
		}
	}

	_simstep = _initSimulation = (unsigned long) round(_simulationTime / _integrator->getTimestepLength() );
	global_log->info() << "Set initial time step to start from to " << _initSimulation << endl;

	global_log->info() << "Initializing output plugins and corresponding output timers" << endl;
	for (auto& outputPlugin : _outputPlugins) {
		global_log->info() << "Initializing output plugin " << outputPlugin->getPluginName() << endl;
		outputPlugin->initOutput(_moleculeContainer, _domainDecomposition, _domain);
		string timer_name = outputPlugin->getPluginName();
		global_simulation->timers()->registerTimer(timer_name,  vector<string>{"SIMULATION_PER_STEP_IO"}, new Timer());
		string timer_output_string = string("Output Plugin ") + timer_name + string(" took:");
		global_simulation->timers()->setOutputString(timer_name, timer_output_string);
	}

	global_log->info() << "System initialised with " << _domain->getglobalNumMolecules() << " molecules." << endl;

	/** Init TemperatureControl beta_trans, beta_rot log-files*/
	if(NULL != _temperatureControl)
		_temperatureControl->InitBetaLogfiles();
}

void Simulation::simulate() {
	global_log->info() << "Started simulation" << endl;

	// (universal) constant acceleration (number of) timesteps
	unsigned uCAT = _pressureGradient->getUCAT();

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
	output(_initSimulation);

	Timer perStepTimer;
	perStepTimer.reset();
	for (_simstep = _initSimulation + 1; _simstep <= _numberOfTimesteps; _simstep++) {
		global_log->debug() << "timestep: " << getSimulationStep() << endl;
		global_log->debug() << "simulation time: " << getSimulationTime() << endl;
		global_simulation->timers()->incrementTimerTimestepCounter();

		computationTimer->start();
		perStepTimer.start();
		/** @todo What is this good for? Where come the numbers from? Needs documentation */
		if (_simstep >= _initGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
				if (!((_simstep + 2 * j + 3) % cpit->getInterval())) {
					cpit->prepareTimestep(_moleculeContainer, _domainDecomposition);
				}
				j++;
			}
		}
		if (_simstep >= _initStatistics) {
		   map<unsigned, CavityEnsemble>::iterator ceit;
		   for(ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++) {
			  if (!((_simstep + 2 * ceit->first + 3) % ceit->second.getInterval())) {
				 ceit->second.preprocessStep();
			  }
		   }
		}

		_integrator->eventNewTimestep(_moleculeContainer, _domain);

		/** @todo FIXME: here two things are done:
		 * 1. the actual kernel to collect rdf data in the handler is activated which would have to be done only once.
		 * 2. the number of sampling steps is incremented, which has to be done each time step.
		 */
		RDF* rdf;
		if ( (_simstep >= _initStatistics) && (nullptr != (rdf = static_cast<RDF*>(getOutputPlugin("RDF")))) ) {
			global_log->debug() << "Activating the RDF sampling" << endl;
			rdf->tickRDF();
			_particlePairsHandler->setRDF(rdf);
			rdf->accumulateNumberOfMolecules(*(getEnsemble()->getComponents()));
		}

		/*! by Stefan Becker <stefan.becker@mv.uni-kl.de> 
		 *realignment tools borrowed from Martin Horsch, for the determination of the centre of mass 
		 *the halo MUST NOT be present*/
#ifndef NDEBUG 
#ifndef ENABLE_MPI
		particleNoTest = _moleculeContainer->getNumberOfParticles();
		global_log->debug()<<"particles before determine shift-methods, halo not present:" << particleNoTest<< "\n";
#endif
#endif
        if(_doAlignCentre && !(_simstep % _alignmentInterval)) {
			if(_componentSpecificAlignment) {
				//! !!! the sequence of calling the two methods MUST be: FIRST determineXZShift() THEN determineYShift() !!!
				_domain->determineXZShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
				_domain->determineYShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
			// edited by Michaela Heier --> realign can be used when LJ93-Potential will be used. Only the shift in the xz-plane will be used. 
			else if(_doAlignCentre && _applyWallFun_LJ_9_3){
				global_log->info() << "realign in the xz-plane without a shift in y-direction\n";
				_domain->determineXZShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
				_domain->noYShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
			else if(_doAlignCentre && _applyWallFun_LJ_10_4){
				global_log->info() << "realign in the xz-plane without a shift in y-direction\n";
				_domain->determineXZShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
				_domain->noYShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
			else{
				_domain->determineShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
#ifndef NDEBUG 
#ifndef ENABLE_MPI			
			particleNoTest = _moleculeContainer->getNumberOfParticles();
			global_log->info()<<"particles after determine shift-methods, halo not present:" << particleNoTest<< "\n";
#endif
#endif
		}
		perStepTimer.stop();
		computationTimer->stop();



#if defined(ENABLE_MPI) && defined(ENABLE_OVERLAPPING)
		bool overlapCommComp = true;
#else
		bool overlapCommComp = false;
#endif

		if (overlapCommComp) {
			performOverlappingDecompositionAndCellTraversalStep(perStepTimer.get_etime());
		}
		else {
			decompositionTimer->start();
			// ensure that all Particles are in the right cells and exchange Particles
			global_log->debug() << "Updating container and decomposition" << endl;
			updateParticleContainerAndDecomposition(perStepTimer.get_etime());
			decompositionTimer->stop();
			perStepTimer.reset();
			double startEtime = computationTimer->get_etime();
			// Force calculation and other pair interaction related computations
			global_log->debug() << "Traversing pairs" << endl;
			computationTimer->start();
			perStepTimer.start();
			forceCalculationTimer->start();
			_moleculeContainer->traverseCells(*_cellProcessor);

			// Update forces in molecules so they can be exchanged
			// TODO: we need some forces to be updated to handle force exchange... but 
			updateForces();
			forceCalculationTimer->stop();
			perStepTimer.stop();
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
		perStepTimer.start();

		measureFLOPRate(_moleculeContainer, _simstep);

		if (_FMM != NULL) {
			global_log->debug() << "Performing FMM calculation" << endl;
			_FMM->computeElectrostatics(_moleculeContainer);
		}



		if(_wall && _applyWallFun_LJ_9_3){
		  _wall->calcTSLJ_9_3(_moleculeContainer, _domain);
		}

		if(_wall && _applyWallFun_LJ_10_4){
		  _wall->calcTSLJ_10_4(_moleculeContainer, _domain);
		}

		if(_mirror && _applyMirror){
		  _mirror->VelocityChange(_moleculeContainer, _domain);
		}

		/** @todo For grand canonical ensemble? Should go into appropriate ensemble class. Needs documentation. */
		// test deletions and insertions
		if (_simstep >= _initGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
				if (!((_simstep + 2 * j + 3) % cpit->getInterval())) {
					global_log->debug() << "Grand canonical ensemble(" << j << "): test deletions and insertions"
							<< endl;
					this->_domain->setLambda(cpit->getLambda());
					this->_domain->setDensityCoefficient(cpit->getDensityCoefficient());
					double localUpotBackup = _domain->getLocalUpot();
					double localVirialBackup = _domain->getLocalVirial();
					cpit->grandcanonicalStep(_moleculeContainer, _domain->getGlobalCurrentTemperature(), this->_domain,
							_cellProcessor);
					_domain->setLocalUpot(localUpotBackup);
					_domain->setLocalVirial(localVirialBackup);
#ifndef NDEBUG
					/* check if random numbers inserted are the same for all processes... */
					cpit->assertSynchronization(_domainDecomposition);
#endif

					int localBalance = cpit->getLocalGrandcanonicalBalance();
					int balance = cpit->grandcanonicalBalance(_domainDecomposition);
					global_log->debug() << "   b[" << ((balance > 0) ? "+" : "") << balance << "("
							<< ((localBalance > 0) ? "+" : "") << localBalance << ")" << " / c = "
							<< cpit->getComponentID() << "]   " << endl;
					_domain->Nadd(cpit->getComponentID(), balance, localBalance);
				}

				j++;
			}
		}

		if(_simstep >= _initStatistics) {
			map<unsigned, CavityEnsemble>::iterator ceit;
			for(ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++) {

				unsigned cavityComponentID = ceit->first;
				CavityEnsemble & cavEns = ceit->second;

				if (!((_simstep + 2 * cavityComponentID + 3) % cavEns.getInterval())) {
					global_log->debug() << "Cavity ensemble for component " << cavityComponentID << ".\n";

					cavEns.cavityStep(this->_moleculeContainer);
				}

				if( (!((_simstep + 2 * cavityComponentID + 7) % cavEns.getInterval())) ||
					(!((_simstep + 2 * cavityComponentID + 3) % cavEns.getInterval())) ||
					(!((_simstep + 2 * cavityComponentID - 1) % cavEns.getInterval())) ) {

					// warning, return value is ignored!
					/*unsigned long ret = */ cavEns.communicateNumCavities(this->_domainDecomposition);
				}
			}
		}

		global_log->debug() << "Deleting outer particles / clearing halo." << endl;
		_moleculeContainer->deleteOuterParticles();

		/** @todo For grand canonical ensemble? Sould go into appropriate ensemble class. Needs documentation. */
		if (_simstep >= _initGrandCanonical) {
			_domain->evaluateRho(_moleculeContainer->getNumberOfParticles(), _domainDecomposition);
		}

		if (!(_simstep % _collectThermostatDirectedVelocity))
			_domain->calculateThermostatDirectedVelocity(_moleculeContainer);
		if (_pressureGradient->isAcceleratingUniformly()) {
			if (!(_simstep % uCAT)) {
				global_log->debug() << "Determine the additional acceleration" << endl;
				_pressureGradient->determineAdditionalAcceleration(
						_domainDecomposition, _moleculeContainer, uCAT
						* _integrator->getTimestepLength());
			}
			global_log->debug() << "Process the uniform acceleration" << endl;
			_integrator->accelerateUniformly(_moleculeContainer, _domain);
			_pressureGradient->adjustTau(this->_integrator->getTimestepLength());
		}
		_longRangeCorrection->calculateLongRange();
		_longRangeCorrection->writeProfiles(_domainDecomposition, _domain, _simstep);

		/* radial distribution function */
		if (_simstep >= _initStatistics) {
			if (this->_lmu.size() == 0) {
				this->_domain->record_cv();
			}
		}

		global_log->debug() << "Inform the integrator (forces calculated)" << endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);

		// calculate the global macroscopic values from the local values
		global_log->debug() << "Calculate macroscopic values" << endl;
		_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, 
				(!(_simstep % _collectThermostatDirectedVelocity)), Tfactor(_simstep));

		// scale velocity and angular momentum
		if ( !_domain->NVE() && _temperatureControl == NULL) {
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
			else if(_thermostatType == ANDERSEN_THERMOSTAT) { //! the Andersen Thermostat
				//global_log->info() << "Andersen Thermostat" << endl;
				double nuDt = _nuAndersen * _integrator->getTimestepLength();
				//global_log->info() << "Timestep length = " << _integrator->getTimestepLength() << " nuDt = " << nuDt << "\n";
				unsigned numPartThermo = 0; // for testing reasons
				double tTarget;
				double stdDevTrans, stdDevRot;
				if(_domain->severalThermostats()) {
					for (ParticleIterator tM = _moleculeContainer->iteratorBegin(); tM != _moleculeContainer->iteratorEnd(); ++tM) {
						if (_rand.rnd() < nuDt) {
							numPartThermo++;
							int thermostat = _domain->getThermostat(tM->componentid());
							tTarget = _domain->getTargetTemperature(thermostat);
							stdDevTrans = sqrt(tTarget/tM->mass());
							for(unsigned short d = 0; d < 3; d++) {
								stdDevRot = sqrt(tTarget*tM->getI(d));
								tM->setv(d,_rand.gaussDeviate(stdDevTrans));
								tM->setD(d,_rand.gaussDeviate(stdDevRot));
							}
						}
					}
				}
				else{
					tTarget = _domain->getTargetTemperature(0);
					for (ParticleIterator tM = _moleculeContainer->iteratorBegin(); tM != _moleculeContainer->iteratorEnd(); ++tM) {
						if (_rand.rnd() < nuDt) {
							numPartThermo++;
							// action of the anderson thermostat: mimic a collision by assigning a maxwell distributed velocity
							stdDevTrans = sqrt(tTarget/tM->mass());
							for(unsigned short d = 0; d < 3; d++) {
								stdDevRot = sqrt(tTarget*tM->getI(d));
								tM->setv(d,_rand.gaussDeviate(stdDevTrans));
								tM->setD(d,_rand.gaussDeviate(stdDevRot));
							}
						}
					}
				}
				//global_log->info() << "Andersen Thermostat: n = " << numPartThermo ++ << " particles thermostated\n";
			}
		}

		// mheinen 2015-07-27 --> TEMPERATURE_CONTROL
        else if ( _temperatureControl != NULL) {
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


		perStepTimer.stop();
		computationTimer->stop();
		perStepIoTimer->start();

		output(_simstep);

		//! TODO: this should be moved! it is definitely not I/O
		/*! by Stefan Becker <stefan.becker@mv.uni-kl.de> 
		  * realignment tools borrowed from Martin Horsch
		  * For the actual shift the halo MUST NOT be present!
		  */
		if(_doAlignCentre && !(_simstep % _alignmentInterval)) {
			_domain->realign(_moleculeContainer);
#ifndef NDEBUG 
#ifndef ENABLE_MPI
			particleNoTest = 0;
			particleNoTest = _moleculeContainer->getNumberOfParticles();
			cout <<"particles after realign(), halo absent: " << particleNoTest<< "\n";
#endif
#endif
		}

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
	global_log->info() << "Finish output from output plugins" << endl;
	for (auto outputPlugin : _outputPlugins) {
		outputPlugin->finishOutput(_moleculeContainer, _domainDecomposition, _domain);
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

void Simulation::output(unsigned long simstep) {

	int mpi_rank = _domainDecomposition->getRank();

	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		OutputBase* output = (*outputIter);
		global_log->debug() << "Output from " << output->getPluginName() << endl;
		global_simulation->timers()->start(output->getPluginName());
		output->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep, &(_lmu), &(_mcav));
		global_simulation->timers()->stop(output->getPluginName());
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
	if (_FMM != NULL) {
		_FMM->printTimers();
		bhfmm::VectorizedLJP2PCellProcessor * temp = dynamic_cast<bhfmm::VectorizedLJP2PCellProcessor*>(_cellProcessor);
		temp->printTimers();
	}

#ifdef TASKTIMINGPROFILE
	std::string outputFileName = "taskTimings_"
								 + std::to_string(std::time(nullptr))
								 + ".csv";
	_taskTimingProfiler->dump(outputFileName);
#endif /* TASKTIMINGPROFILE */

	if (_domainDecomposition != NULL) {
		delete _domainDecomposition;
		_domainDecomposition = NULL;
	}
	global_simulation = NULL;
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

	global_log->info() << "Creating PressureGradient ... " << endl;
	_pressureGradient = new PressureGradient(ownrank);

	global_log->info() << "Creating domain ..." << endl;
	_domain = new Domain(ownrank, this->_pressureGradient);
	global_log->info() << "Creating ParticlePairs2PotForceAdapter ..." << endl;
	_particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);

	this->_mcav = map<unsigned, CavityEnsemble>();

	global_log->info() << "Initialization done" << endl;
}

OutputBase* Simulation::getOutputPlugin(const std::string& name)  {
	for(auto& outputPlugin : _outputPlugins) {
		if(name.compare(outputPlugin->getPluginName()) == 0) {
			return outputPlugin;
		}
	}
	return nullptr;
}

void Simulation::measureFLOPRate(ParticleContainer* cont, unsigned long simstep) {
	OutputBase * flopRateBase = getOutputPlugin("FlopRateWriter");
	if (flopRateBase == nullptr) {
		return;
	}

	FlopRateWriter * flopRateWriter = dynamic_cast<FlopRateWriter * >(flopRateBase);
	mardyn_assert(flopRateWriter != nullptr);

	flopRateWriter->measureFLOPS(cont, simstep);
}

unsigned long Simulation::getNumberOfTimesteps() const {
	return _numberOfTimesteps;
}
