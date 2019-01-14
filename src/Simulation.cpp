#define SIMULATION_SRC
#include "Simulation.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "Common.h"
#include "Domain.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/AdaptiveSubCells.h"
#include "parallel/DomainDecompBase.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#else
#include "parallel/DomainDecompDummy.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/adapter/LJFlopCounter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"

#include "io/io.h"
#include "io/GeneratorFactory.h"
#include "io/RDF.h"
#include "io/TcTS.h"
#include "io/Mkesfera.h"

#include "ensemble/GrandCanonical.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/PressureGradient.h"
#include "ensemble/CavityEnsemble.h"

#include "thermostats/VelocityScalingThermostat.h"

#include "utils/OptionParser.h"
#include "utils/Timer.h"
#include "utils/Logger.h"

#include "molecules/Molecule.h"

using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;

Simulation* global_simulation;

Simulation::Simulation()
	: _simulationTime(0),
	_initStatistics(0),
	_rdf(NULL),
	_ljFlopCounter(NULL),
	_domainDecomposition(NULL),
	_forced_checkpoint_time(0) {
	_ensemble = new CanonicalEnsemble();
	initialize();
}

Simulation::~Simulation() {
	delete _domainDecomposition;
	delete _pressureGradient;
	delete _domain;
	delete _particlePairsHandler;
	delete _cellProcessor;
	delete _moleculeContainer;
	delete _integrator;
	delete _inputReader;
	delete _ljFlopCounter;
	delete _thismol;
}

void Simulation::exit(int exitcode) {
#ifdef ENABLE_MPI
	// terminate all mpi processes and return exitcode
	MPI_Abort(MPI_COMM_WORLD, exitcode);
#else
	// call global exit
	::exit(exitcode);
#endif
}



void Simulation::readXML(XMLfileUnits& xmlconfig) {
	/* integrator */
	string integratorType;
	xmlconfig.getNodeValue("integrator@type", integratorType);
	global_log->info() << "Integrator type: " << integratorType << endl;
	if(integratorType == "Leapfrog") {
		_integrator = new Leapfrog();
	} else {
		global_log-> error() << "Unknown integrator " << integratorType << endl;
		this->exit(1);
	}
	if(xmlconfig.changecurrentnode("integrator")) {
		_integrator->readXML(xmlconfig);
		_integrator->init();
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->error() << "Integrator section missing." << endl;
	}

	/* steps */
	xmlconfig.getNodeValue("run/production/steps", _numberOfTimesteps);
	global_log->info() << "Number of timesteps: " << _numberOfTimesteps << endl;
	xmlconfig.getNodeValue("run/equilibration/steps", _initStatistics);
	global_log->info() << "Number of equilibration steps: " << _initStatistics << endl;
	xmlconfig.getNodeValueReduced("run/currenttime", _simulationTime);
	global_log->info() << "Simulation start time: " << _simulationTime << endl;

	/* enseble */
	string ensembletype;
	xmlconfig.getNodeValue("ensemble@type", ensembletype);
	global_log->info() << "Ensemble: " << ensembletype<< endl;
	if (ensembletype == "NVT") {
		_ensemble = new CanonicalEnsemble();
	} else if (ensembletype == "muVT") {
		global_log->error() << "muVT ensemble not completely implemented via XML input." << endl;
		this->exit(1);
// 		_ensemble = new GrandCanonicalEnsemble();
	} else {
		global_log->error() << "Unknown ensemble type: " << ensembletype << endl;
		this->exit(1);
	}
	if(xmlconfig.changecurrentnode("ensemble")) {
		_ensemble->readXML(xmlconfig);
		/* store data in the _domain member as long as we do not use the ensemble everywhere */
		for (int d = 0; d < 3; d++) {
			_domain->setGlobalLength(d, _ensemble->domain()->length(d));
		}
		_domain->setGlobalTemperature(_ensemble->T());
		xmlconfig.changecurrentnode("..");
	}
	else {
		global_log->error() << "Ensemble section missing." << endl;
	}

	/* algorithm */
	if(xmlconfig.changecurrentnode("algorithm")) {
		/* cutoffs */
		if (xmlconfig.getNodeValueReduced("cutoffs/radiusLJ", _LJCutoffRadius)) {
			_cutoffRadius = _LJCutoffRadius;
			global_log->info() << "dimensionless LJ cutoff radius:\t"
					<< _LJCutoffRadius << endl;
			global_log->info() << "dimensionless cutoff radius:\t"
					<< _cutoffRadius << endl;
		} else {
			global_log->error() << "Cutoff section missing." << endl;
			this->exit(1);
		}

		double epsilonRF = 0;
		xmlconfig.getNodeValueReduced("electrostatic[@type='ReactionField']/epsilon", epsilonRF);
		global_log->info() << "Epsilon Reaction Field: " << epsilonRF << endl;
		_domain->setepsilonRF(epsilonRF);

		/* parallelization */
		string parallelisationtype("DomainDecomposition");
		xmlconfig.getNodeValue("parallelisation@type", parallelisationtype);
		global_log->info() << "Parallelisation type: " << parallelisationtype << endl;
		if(parallelisationtype == "DomainDecomposition") {
	#ifdef ENABLE_MPI
			_domainDecomposition = new DomainDecomposition();
	#else
		global_log->error() << "DomainDecomposition not available in sequential mode." << endl;
	#endif
		}
		else if(parallelisationtype == "KDDecomposition") {
	#ifdef ENABLE_MPI
			_domainDecomposition = new KDDecomposition(getcutoffRadius(), _domain);
	#else
		global_log->error() << "KDDecomposition not available in sequential mode." << endl;
	#endif
		}
		else if(parallelisationtype == "DummyDecomposition") {
	#ifdef ENABLE_MPI
		global_log->error() << "DummyDecomposition not available in parallel mode." << endl;
	#else
			_domainDecomposition = new DomainDecompDummy();
	#endif
		}
		else {
			global_log->error() << "Unknown parallelisation type: " << parallelisationtype << endl;
			this->exit(1);
		}

		/* datastructure */
		string datastructuretype;
		xmlconfig.getNodeValue("datastructure@type", datastructuretype);
		global_log->info() << "Datastructure type: " << datastructuretype << endl;
		if(datastructuretype == "LinkedCells") {
			_moleculeContainer = new LinkedCells();
			_particleContainerType = LINKED_CELL; /* TODO: Necessary? */
			global_log->info() << "Setting cell cutoff radius for linked cell datastructure to " << _cutoffRadius << endl;
			LinkedCells *lc = static_cast<LinkedCells*>(_moleculeContainer);
			lc->setCutoff(_cutoffRadius);
		}
		else if(datastructuretype == "AdaptiveSubCells") {
			_moleculeContainer = new AdaptiveSubCells();
			_particleContainerType = ADAPTIVE_LINKED_CELL; /* TODO: Necessary? */
		}
		else {
			global_log->error() << "Unknown data structure type: " << datastructuretype << endl;
			this->exit(1);
		}
		if(xmlconfig.changecurrentnode("datastructure")) {
			_moleculeContainer->readXML(xmlconfig);
			double bBoxMin[3];
			double bBoxMax[3];
			/* TODO: replace Domain with DomainBase. */
			_domainDecomposition->getBoundingBoxMinMax(_domain, bBoxMin, bBoxMax);
			_moleculeContainer->rebuild(bBoxMin, bBoxMax);
			xmlconfig.changecurrentnode("..");
		} else {
			global_log->error() << "Datastructure section missing" << endl;
			this->exit(1);
		}

	#ifndef ENABLE_MPI
		if (parallelisationtype != "DummyDecomposition") {
			global_log->warning()
				<< "Input demands parallelization, but the current compilation doesn't support parallel execution."
				<< endl;
		}
	#endif

		if(xmlconfig.changecurrentnode("parallelisation")) {
			_domainDecomposition->readXML(xmlconfig);
			xmlconfig.changecurrentnode("..");
		}
		else {
			global_log->warning() << "Parallelisation section missing." << endl;
		}

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
						componentId = getEnsemble()->component(componentName)->ID();
						int thermostatID = _domain->getThermostat(componentId);
						_domain->setTargetTemperature(thermostatID, temperature);
						global_log->info() << "Adding velocity scaling thermostat for component '" << componentName << "' (ID: " << componentId << "), T = " << temperature << endl;
					}
				}
				else {
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
	XMLfile::Query::const_iterator outputPluginIter;
	for( outputPluginIter = query.begin(); outputPluginIter; outputPluginIter++ ) {
		xmlconfig.changecurrentnode( outputPluginIter );
		OutputBase *outputPlugin;
		string pluginname("");
		xmlconfig.getNodeValue("@name", pluginname);
		global_log->info() << "Enabling output plugin: " << pluginname << endl;
		if(pluginname == "CheckpointWriter") {
			outputPlugin = new CheckpointWriter();
		}
		else if(pluginname == "DecompWriter") {
			outputPlugin = new DecompWriter();
		}
		else if(pluginname == "LJFLOPCounter") {
			/** @todo  Make the LJ Flop counter a real output plugin */
			_ljFlopCounter = new LJFlopCounter(_LJCutoffRadius);
			continue;
		}
		else if(pluginname == "MmspdWriter") {
			outputPlugin = new MmspdWriter();
		}
		else if(pluginname == "PovWriter") {
			outputPlugin = new PovWriter();
		}
		else if(pluginname == "RDF") {
			_rdf = new RDF();
			outputPlugin = _rdf;
		}
		else if(pluginname == "Resultwriter") {
			outputPlugin = new ResultWriter();
		}
		else if(pluginname == "SysMonOutput") {
			outputPlugin = new SysMonOutput();
		}
		else if(pluginname == "VISWriter") {
			outputPlugin = new VISWriter();
		}
#ifdef VTK
		else if(pluginname == "VTKMoleculeWriter") {
			outputPlugin = new VTKMoleculeWriter();
		}
		else if(pluginname == "VTKGridWriter") {
			outputPlugin = new VTKGridWriter();
		}
#endif /* VTK */
		else if(pluginname == "XyzWriter") {
			outputPlugin = new XyzWriter();
		}
		/* temporary */
		else if(pluginname == "MPICheckpointWriter") {
			outputPlugin = new MPICheckpointWriter();
		}
		else {
			global_log->warning() << "Unknown plugin " << pluginname << endl;
			continue;
		}
		outputPlugin->readXML(xmlconfig);
		_outputPlugins.push_back(outputPlugin);
	}
	xmlconfig.changecurrentnode(oldpath);
}

void Simulation::readConfigFile(string filename) {
	if (filename.rfind(".xml") == filename.size() - 4) {
		global_log->info() << "command line config file type is XML (*.xml)"
				<< endl;
		initConfigXML(filename);
	} else if (filename.rfind(".cfg") == filename.size() - 4) {
		global_log->info()
				<< "command line config file type is oldstyle (*.cfg)" << endl;
		initConfigOldstyle(filename);
	} else {
		global_log->info()
				<< "command line config file type is unknown: trying oldstyle"
				<< endl;
		initConfigOldstyle(filename);
	}
}

void Simulation::initConfigXML(const string& inputfilename) {
	global_log->info() << "init XML config file: " << inputfilename << endl;
	XMLfileUnits inp(inputfilename);

	global_log->debug() << "Input XML:" << endl << string(inp) << endl;

	inp.changecurrentnode("/mardyn");
	string version("unknown");
	inp.getNodeValue("@version", version);
	global_log->info() << "MarDyn XML config file version: " << version << endl;

	if (inp.changecurrentnode("simulation")) {
		string siminpfile;
		int numsimpfiles = inp.getNodeValue("input", siminpfile);
		if (numsimpfiles == 1) {
			string siminptype;
			global_log->info() << "Reading input file: " << siminpfile << endl;
			inp.getNodeValue("input@type", siminptype);
			global_log->info() << "Input file type: " << siminptype << endl;
			if (siminptype == "oldstyle") {
				initConfigOldstyle(siminpfile);
				/* Skip the rest of the xml config for old cfg files. */
				return;
			} else {
				global_log->error() << "Unknown input file type: " << siminptype << endl;
				this->exit(1);;
			}
		} else if (numsimpfiles > 1) {
			global_log->error() << "Multiple input file sections are not supported." << endl;
			this->exit(1);
		}

		readXML(inp);

		string pspfile;
		if (inp.getNodeValue("ensemble/phasespacepoint/file", pspfile)) {
			pspfile.insert(0, inp.getDir());
			global_log->info() << "phasespacepoint description file:\t"
					<< pspfile << endl;

			string pspfiletype("ASCII");
			inp.getNodeValue("ensemble/phasespacepoint/file@type", pspfiletype);
			global_log->info() << "       phasespacepoint file type:\t"
					<< pspfiletype << endl;
			if (pspfiletype == "ASCII") {
				_inputReader = (InputBase*) new InputOldstyle();
				_inputReader->setPhaseSpaceFile(pspfile);
			}
		}
		string oldpath = inp.getcurrentnodepath();
		if(inp.changecurrentnode("ensemble/phasespacepoint/generator")) {
			string generatorName;
			inp.getNodeValue("@name", generatorName);
			global_log->info() << "Generator: " << generatorName << endl;
			if(generatorName == "GridGenerator") {
				_inputReader = new GridGenerator();
			}
			else if(generatorName == "mkesfera") {
				_inputReader = new MkesferaGenerator();
			}
			else if(generatorName == "mkTcTS") {
				_inputReader = new MkTcTSGenerator();
			}
			else {
				global_log->error() << "Unknown generator: " << generatorName << endl;
				exit(1);
			}
			_inputReader->readXML(inp);
		}
		inp.changecurrentnode(oldpath);


		inp.changecurrentnode("..");
	} // simulation-section
	else {
		global_log->error() << "Simulation section missing" << endl;
		exit(1);
	}


	// read particle data
	unsigned long maxid = _inputReader->readPhaseSpace(_moleculeContainer,
			&_lmu, _domain, _domainDecomposition);


	_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
	_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

	// test new Decomposition
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
		double tmp_molecularMass = global_simulation->getEnsemble()->component(cpit->getComponentID())->m();
		cpit->setSystem(_domain->getGlobalLength(0),
				_domain->getGlobalLength(1), _domain->getGlobalLength(2),
				tmp_molecularMass);
		cpit->setGlobalN(global_simulation->getEnsemble()->component(cpit->getComponentID())->getNumMolecules());
		cpit->setNextID(j + (int) (1.001 * (256 + maxid)));

		cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0),
				_moleculeContainer->getBoundingBoxMax(0),
				_moleculeContainer->getBoundingBoxMin(1),
				_moleculeContainer->getBoundingBoxMax(1),
				_moleculeContainer->getBoundingBoxMin(2),
				_moleculeContainer->getBoundingBoxMax(2));
		/* TODO: thermostat */
		double Tcur = _domain->getCurrentTemperature(0);
		/* FIXME: target temperature from thermostat ID 0 or 1?  */
		double
				Ttar =
						_domain->severalThermostats() ? _domain->getTargetTemperature(
								1)
								: _domain->getTargetTemperature(0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;
		cpit->submitTemperature(Tcur);
		if (h != 0.0)
			cpit->setPlanckConstant(h);

		j++;
	}
}

void Simulation::initConfigOldstyle(const string& inputfilename) {
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif

	global_log->info() << "init oldstyle config file: " << inputfilename
			<< endl;
			
	// initialize result folder
	mkdir("./Results", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);  

	// open filestream to the input file
	ifstream inputfilestream(inputfilename.c_str());
	if (!inputfilestream.is_open()) {
		global_log->error() << "Could not open file " << inputfilename << endl;
		exit(1);
	}

	//  std::string inputPath;
	//  unsigned int lastIndex = inputfilename.find_last_of('/',inputfilename.size()-1);
	//  if (lastIndex == string::npos)
	//    inputPath="";
	//  else
	//    inputPath = inputfilename.substr(0, lastIndex+1);


	// used to store one token of the inputfilestream
	string token;
	string token2;

	double timestepLength;
	unsigned cosetid = 0;
        bool widom = false;
        bool ultra = false;
        map<unsigned, double*> cavity_grid;

	// The first line of the config file has to contain the token "MDProjectConfig"
	inputfilestream >> token;
	if ((token != "mardynconfig") && (token != "MDProjectConfig")) {
		global_log->error() << "Not a mardynconfig file! First token: "
				<< token << endl;
		exit(1);
	}

	while (inputfilestream) {
		token.clear();
		inputfilestream >> token;
		global_log->debug() << " [[" << token << "]]" << endl;

		if (token.substr(0, 1) == "#") {
			inputfilestream.ignore(std::numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		if (token == "phaseSpaceFile") {
			string phaseSpaceFileFormat;
			inputfilestream >> phaseSpaceFileFormat;

			if (timestepLength == 0.0) {
				global_log->error() << "timestep missing." << endl;
				exit(1);
			}
			if (phaseSpaceFileFormat == "OldStyle") {
				string phaseSpaceFileName;
				inputfilestream >> phaseSpaceFileName;
				_inputReader = (InputBase*) new InputOldstyle();
				_inputReader->setPhaseSpaceFile(phaseSpaceFileName);
				_inputReader->setPhaseSpaceHeaderFile(phaseSpaceFileName);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else if (phaseSpaceFileFormat == "Generator") {
				global_log->info() << "phaseSpaceFileFormat is Generator!"
						<< endl;
				string generatorName;  // name of the library to load
				string inputFile;  // name of the input file for the generator

				string line;
				getline(inputfilestream, line);
				stringstream lineStream(line);
				lineStream >> generatorName >> inputFile;
				_inputReader = GeneratorFactory::loadGenerator(generatorName,
						inputFile);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else {
				global_log->error() << "Don't recognize phasespaceFile reader "
						<< phaseSpaceFileFormat << endl;
				exit(1);
			}
			if (_LJCutoffRadius == 0.0)
				_LJCutoffRadius = _cutoffRadius;
			_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
		} else if (token == "timestepLength") {
			inputfilestream >> timestepLength;
			setTimeStepLength(timestepLength);
		} else if (token == "cutoffRadius") {
			inputfilestream >> _cutoffRadius;
		} else if (token == "LJCutoffRadius") {
			inputfilestream >> _LJCutoffRadius;
		} else if ((token == "parallelization") || (token == "parallelisation")) {
#ifndef ENABLE_MPI
			global_log->warning()
					<< "Input file demands parallelization, but the current compilation doesn't\n\tsupport parallel execution.\n"
					<< endl;
			inputfilestream >> token;
#else
			inputfilestream >> token;
			if (token == "DomainDecomposition") {
				// default DomainDecomposition is already set in initialize();
				//_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
			}
			else if(token == "KDDecomposition") {
				delete _domainDecomposition;
				int updateFrequency = 100;
				int fullSearchThreshold = 3;
				string line;
				getline(inputfilestream, line);
				stringstream lineStream(line);
				lineStream >> updateFrequency >> fullSearchThreshold;
				_domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, updateFrequency, fullSearchThreshold);
			}
#endif
		} else if (token == "datastructure") {

			if (_domainDecomposition == NULL) {
				global_log->error()
						<< "_domainDecomposition is NULL! Probably you compiled for MPI, but didn't specify line \"parallelization\" before line \"datastructure\"!"
						<< endl;
				exit(1);
			}

			inputfilestream >> token;
			if (token == "LinkedCells") {
				_particleContainerType = LINKED_CELL;  /* TODO: Necessary? */
				int cellsInCutoffRadius;
				inputfilestream >> cellsInCutoffRadius;
				double bBoxMin[3];
				double bBoxMax[3];
				for (int i = 0; i < 3; i++) {
					bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i,
							_domain);
					bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i,
							_domain);
				}
				if (this->_LJCutoffRadius == 0.0)
					_LJCutoffRadius = this->_cutoffRadius;
				_moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius,
				        cellsInCutoffRadius);
			} else if (token == "AdaptiveSubCells") {
				_particleContainerType = ADAPTIVE_LINKED_CELL; /* TODO: Necessary? */
				double bBoxMin[3];
				double bBoxMax[3];
				for (int i = 0; i < 3; i++) {
					bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i,_domain);
					bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i,_domain);
				}
				// creates a new Adaptive SubCells datastructure
				if (_LJCutoffRadius == 0.0)
					_LJCutoffRadius = _cutoffRadius;
					_moleculeContainer = new AdaptiveSubCells(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius);
			} else {
				global_log->error() << "UNKOWN DATASTRUCTURE: " << token
						<< endl;
				exit(1);
			}
		} else if (token == "output") {
			inputfilestream >> token;
			// default output file format setting
			string all ("all");
			_domain->specifyOutputFormat(all);
			if (token == "ResultWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new ResultWriter(writeFrequency,
						outputPathAndPrefix));
				global_log->debug() << "ResultWriter '" << outputPathAndPrefix
						<< "'.\n";
			} else if (token == "XyzWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new XyzWriter(writeFrequency,
						outputPathAndPrefix, true));
				global_log->debug() << "XyzWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "justCavityFiles") { // suboption of XyzWriter: suppresses the output of xyz-profiles --> data reduction
				_noXYZ = true;
			} else if (token == "CheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new CheckpointWriter(writeFrequency,
						outputPathAndPrefix, true));
				global_log->debug() << "CheckpointWriter " << writeFrequency
						<< " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "PovWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new PovWriter(writeFrequency,
						outputPathAndPrefix, true));
				global_log->debug() << "POVWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "DecompWriter") {
				unsigned long writeFrequency;
				string mode;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> mode
						>> outputPathAndPrefix;
				_outputPlugins.push_back(new DecompWriter(writeFrequency, mode,
						outputPathAndPrefix, true));
				global_log->debug() << "DecompWriter " << writeFrequency
						<< " '" << outputPathAndPrefix << "'.\n";
			} else if ((token == "VisittWriter") || (token == "VISWriter")) {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new VISWriter(writeFrequency,
						outputPathAndPrefix, true));
				global_log->debug() << "VISWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "VTKWriter") {
#ifdef VTK

				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new VTKMoleculeWriter(writeFrequency,
						outputPathAndPrefix));
				global_log->debug() << "VTKWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
#else
				Log::global_log->error() << std::endl << "VTK-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			} else if (token == "VTKGridWriter") {
#ifdef VTK
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;

				if (_particleContainerType == LINKED_CELL) {
					_outputPlugins.push_back(new VTKGridWriter(writeFrequency,
							outputPathAndPrefix));
					global_log->debug() << "VTKGridWriter " << writeFrequency
							<< " '" << outputPathAndPrefix << "'.\n";
				} else {
					global_log->warning()
							<< "VTKGridWriter only supported with LinkedCells!"
							<< std::endl;
					global_log->warning()
							<< "Generating no VTK output for the grid!"
							<< std::endl;
				}
#else
				Log::global_log->error() << std::endl << "VTK-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			}
			// by Stefan Becker <stefan.becker@mv.uni-kl.de>
			// output for the MegaMol Simple Particle Data File Format (*.mmspd)
			else if (token == "MmspdWriter"){
			      unsigned long writeFrequency = 0;
			      string outputPathAndPrefix;
			      inputfilestream >> writeFrequency >> outputPathAndPrefix;
			      _outputPlugins.push_back(new MmspdWriter(writeFrequency, outputPathAndPrefix));
			      global_log->debug() << "MmspdWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			}
			// temporary
			else if (token == "MPICheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new MPICheckpointWriter(writeFrequency,
						outputPathAndPrefix, true));
				global_log->debug() << "MPICheckpointWriter " << writeFrequency
						<< " '" << outputPathAndPrefix << "'.\n";
			}
			else {
				global_log->warning() << "Unknown output plugin " << token << endl;
			}
		} else if (token == "accelerate") {
			cosetid++;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "comp") {
				global_log->error() << "Expected 'comp' instead of '" << token
						<< "'.\n";
				exit(1);
			}
			int cid = 0;
			while (cid >= 0) {
				inputfilestream >> cid;
				if (cid > 0)
					global_log->info() << "acc. for component " << cid << endl;
				cid--;
				_pressureGradient->assignCoset((unsigned) cid, cosetid);
			}
			double v;
			inputfilestream >> v;
			global_log->debug() << "velocity " << v << endl;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "towards") {
				global_log->error() << "Expected 'towards' instead of '"
						<< token << "'.\n";
				exit(1);
			}
			double dir[3];
			double dirnorm = 0;
			for (unsigned d = 0; d < 3; d++) {
				inputfilestream >> dir[d];
				dirnorm += dir[d] * dir[d];
			}
			dirnorm = 1.0 / sqrt(dirnorm);
			for (unsigned d = 0; d < 3; d++)
				dir[d] *= dirnorm;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "within") {
				global_log->error() << "Expected 'within' instead of '"
						<< token << "'.\n";
				exit(1);
			}
			double tau;
			inputfilestream >> tau;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "from") {
				global_log->error() << "Expected 'from' instead of '" << token
						<< "'.\n";
				exit(1);
			}
			double ainit[3];
			for (unsigned d = 0; d < 3; d++)
				inputfilestream >> ainit[3];
			if (timestepLength == 0.0) {
				global_log->error() << "timestep missing." << endl;
				exit(1);
			}
			_pressureGradient->specifyComponentSet(cosetid, dir, tau, ainit,
					timestepLength);
		} else if (token == "constantAccelerationTimesteps") {
			unsigned uCAT;
			inputfilestream >> uCAT;
			_pressureGradient->setUCAT(uCAT);
		} else if (token == "zetaFlow") {
			double zeta;
			inputfilestream >> zeta;
			_pressureGradient->setZetaFlow(zeta);
		} else if (token == "tauPrimeFlow") {
			double tauPrime;
			inputfilestream >> tauPrime;
			if (timestepLength == 0.0) {
				global_log->error() << "timestep missing." << endl;
				exit(1);
			}
			_pressureGradient->specifyTauPrime(tauPrime, timestepLength);
		} else if (token == "RDF") {
			double interval;
			unsigned bins;
			inputfilestream >> interval >> bins;
			if (global_simulation->getEnsemble()->components()->size() <= 0) {
				global_log->error()
						<< "PhaseSpaceFile-Specifiation has to occur befor RDF-Token!"
						<< endl;
				exit(-1);
			}
			_rdf = new RDF(interval, bins, global_simulation->getEnsemble()->components());
			_outputPlugins.push_back(_rdf);
		} else if (token == "RDFOutputTimesteps") { /* TODO: suboption of RDF */
			unsigned int RDFOutputTimesteps;
			inputfilestream >> RDFOutputTimesteps;
			_rdf->setOutputTimestep(RDFOutputTimesteps);
		} else if (token == "RDFOutputPrefix") { /* TODO: suboption of RDF */
			std::string RDFOutputPrefix;
			inputfilestream >> RDFOutputPrefix;
			_rdf->setOutputPrefix(RDFOutputPrefix);
		} else if(token == "OutputFormat") { // determination of the output format
			string outputFormat;
			inputfilestream >> outputFormat;
			_domain->specifyOutputFormat(outputFormat);
			string matlab ("matlab");
			string vtk ("vtk");
			string all ("all");
			if(matlab.compare(outputFormat) == 0 || vtk.compare(outputFormat) == 0 || all.compare(outputFormat) == 0)
			  global_log->info() << "File format for the output is " << outputFormat << endl;
			else{
			  global_log->info() << "Invalid type of output file format \'" << outputFormat << "\' found. [Allowed: matlab, vtk, all]" << endl;
			  global_log->info() << "Continuing with output file format 'all'." << endl;
			  global_log->error() << "Invalid type of output file format \'" << outputFormat << "\' found. [Allowed: matlab, vtk, all]" << endl;
			}
		} else if (token == "EnergyOutputPerComponent") {
			_doEnergyOutputPerComponent = true;	
		} else if (token == "profile") {   // 
			unsigned xun, yun, zun;
			inputfilestream >> xun >> yun >> zun;
			_domain->setupProfile(xun, yun, zun);
			_doRecordProfile = true;
		} else if (token == "profileRecordingTimesteps") { /* TODO: suboption of profile */
			inputfilestream >> _profileRecordingTimesteps;
		} else if (token == "profileOutputTimesteps") { /* TODO: suboption of profile */
			inputfilestream >> _profileOutputTimesteps;	
		} else if (token == "profiledComponent") { /* TODO: suboption of profile, check if required to enable output in general */
			unsigned cid;
			inputfilestream >> cid;
			cid--;
			_domain->considerComponentInProfile(cid);
		} else if (token == "profileOutputPrefix") { /* TODO: suboption of profile */
			inputfilestream >> _profileOutputPrefix;
		} else if (token == "confinementDensity"){	/* variation of profile: records discrete information about temperature, density and velocity averaged over in-plane-direction in polar coordinates around the lower asperity; needs: _domain->setupProfile(xun, yun, zun), _profileRecordingTimesteps, _profileOutputTimesteps */
			double radius1, radius2, xCentre, yCentre;
			inputfilestream >> radius1 >> radius2 >> xCentre >> yCentre;
			this->_domain->confinementDensity(radius1, radius2, xCentre, yCentre);				  
		} else if (token == "directedVelocityTimeSpan") {	// if this option is given directed velocities are by default calculated as a moving average over a certain time span; if an additional suboption is chosen for the directed velocities, here, the time span for the averaging is determined
			unsigned dirVelTime;
			inputfilestream >> dirVelTime;
			_directedVelocityTime = dirVelTime;
			_boolDirectedVel = true;
			global_log->info() << "Directed velocity is calculated each " << _directedVelocityTime << " timesteps." << endl;
		} else if (token == "movingAverageSigma2D") {	// if this option is given directed velocities are calculated as a moving average in size 1*1*z_max over a certain time span
			_doMovingAverageSigma2D = true;
			global_log->info() << "Directed velocity is calculated as moving average 'Sigma 2D'." << endl;	
		} else if (token == "simpleAverage") {	// if this option is given directed velocities are calculated as a simple average over a certain time span; good for large parallel simulations
			_doSimpleAverage = true;
			global_log->info() << "Directed velocity is calculated as simple average." << endl;
		} else if (token == "simpleAverageSigma2D") {	// if this option is given directed velocities are calculated as a simple average in size 1*1*z_max over a certain time span; good for large parallel simulations
			_doSimpleAverageSigma2D = true;
			global_log->info() << "Directed velocity is calculated as simple average 'Sigma 2D'." << endl;
		} else if (token == "simpleAverageSigma3D") {	// if this option is given directed velocities are calculated as a simple average in size 1*1*1 over a certain time span; good for large parallel simulations
			_doSimpleAverageSigma3D = true;
			global_log->info() << "Directed velocity is calculated as simple average 'Sigma 3D'." << endl;	
		} else if (token == "neighbourListAverage") {	// if this option is given directed velocities are calculated as a moving averaged over a certain time span and over the neighbour list; good for large parallel simulations
			_doNeighbourAverage = true;
			global_log->info() << "Directed velocity is calculated as moving average over the neighbour list." << endl;
		} else if (token == "slabProfile") {	// records discrete information about temperature, density and velocity averaged over in-plane-direction
			unsigned xun, yun, zun;
			inputfilestream >> xun >> yun >> zun;
			_domain->setupSlabProfile(xun, yun, zun);
			_doRecordSlabProfile = true;	
		} else if (token == "slabProfileRecordingTimesteps") { /* TODO: suboption of slabProfile */
			inputfilestream >> _slabProfileRecordingTimesteps;
		} else if (token == "slabProfileOutputTimesteps") { /* TODO: suboption of slabProfile */
			inputfilestream >> _slabProfileOutputTimesteps;	
		} else if (token == "slabProfiledComponent") { /* TODO: suboption of slabProfile, check if required to enable output in general */
			unsigned cid;
			inputfilestream >> cid;
			cid--;
			_domain->considerComponentInProfileSlab(cid);
		} else if (token == "slabProfileOutputPrefix") { /* TODO: suboption of slabProfile */
			inputfilestream >> _slabProfileOutputPrefix;
		} else if (token == "slabNoVelocity") { /* TODO: suboption of slabProfile */
			_reduceDataSlab = true;
		} else if (token == "stress") {	   // calculates virial stresses in the desired component
			unsigned xun, yun, zun;
			string stress, weightingFunc;
			bool properties[6];
			for (int prop_n = 0; prop_n < 6; prop_n++)
				properties[prop_n] = true;
			
			inputfilestream >> xun >> yun >> zun >> stress;
			_doRecordStressProfile = true;
			if(stress == "Hardy" && zun == 1){ // calculates Hardy Stresses; just possible in 2D(x,y)
			  _HardyStress = true;
			  inputfilestream >> token;
			  if (token != "Linear" && token != "Pyramide" ) {
				global_log->error() << "Expected Hardy stress option 'Linear' or 'Pyramide' instead of '"
						<< token << "' for the weighting function.\n";
				exit(1);
			  }
			  weightingFunc = token;
			  _weightingStress = weightingFunc;
			}else if(stress == "Hardy" && zun != 1){
				global_log->info() << "Hardy stresses can just be calculated in 2D so far. Therefore, the calculation method is switched to 'Virial' automatically!\n";
				stress == "Virial";
			}
			global_log->info() << "Stress profile: " << stress << endl;
			_domain->setupStressProfile(xun, yun, zun, properties);
		} else if (token == "stressProperties"){ // calculate properties: p = hydrodynamic pressure, normalStress = normal stresses, shearStress = shear stresses, misesStress = Mises stress, q = heatflux, D = diffusion
			bool properties[6];
			int numberOfProperties = 0;
			int count = 0;
			for (int prop_n = 0; prop_n < 6; prop_n++)
				properties[prop_n] = false;
			  inputfilestream >> numberOfProperties;
			  for (count = 0; count < numberOfProperties; count++){
				inputfilestream >> token2;
				if (token2 == "p") {
					properties[0] = true;
					global_log->info() << "Stress property: " << token2 << endl;
				}else if (token2 == "normalStress") {
					properties[1] = true;
					global_log->info() << "Stress property: " << token2 << endl;
				}else if (token2 == "shearStress") {
					properties[2] = true;  
					global_log->info() << "Stress property: " << token2 << endl;
				}else if (token2 == "misesStress") {
					properties[3] = true;   
					global_log->info() << "Stress property: " << token2 << endl;
				}else if (token2 == "q") {
					properties[4] = true;
					global_log->info() << "Stress property: " << token2 << endl;
				}else if (token2 == "D") {
					properties[5] = true;  
					global_log->info() << "Stress property: " << token2 << endl;
				}else{
					global_log->error() << "Number of calculated properties doesn't match the expectation or you chose the wrong indentifier!\n";
					global_log->error() << "Accepted indentifier: p normalStress shearStress misesStress q D\n";
					exit(1);
				}
			  }
			_domain->setupStressProperties(properties);
		} else if (token == "stressProfileRecordingTimesteps") { /* TODO: suboption of stressProfile */
			inputfilestream >> _stressProfileRecordingTimesteps;
			_domain->setupStressRecTime(_stressProfileRecordingTimesteps);
		} else if (token == "stressProfileOutputTimesteps") { /* TODO: suboption of stressProfile */
			inputfilestream >> _stressProfileOutputTimesteps;	
			_domain->setupStressOutTime(_stressProfileOutputTimesteps);
		} else if (token == "stressProfiledComponent") { /* TODO: suboption of stressProfile, check if required to enable output in general */
			unsigned cid;
			inputfilestream >> cid;
			cid--;
			_domain->enableStressCalculation(cid);
		} else if (token == "stressProfileOutputPrefix") { /* TODO: suboption of stressProfile */
			inputfilestream >> _stressProfileOutputPrefix;	
		} else if (token == "bulkPressure") {	// calculates the pressure, temperature and density in the bulk fluid (due to desired component and control volume)
			double xmin, xmax, ymin, ymax;
			unsigned cid;
			inputfilestream >> xmin >> xmax >> ymin >> ymax >> cid;
			cid--;
			_doRecordBulkPressure = true;
			_domain->setupBulkPressure(xmin, xmax, ymin, ymax, cid);
		} else if (token == "bulkPressureOutputTimesteps") { /* TODO: suboption of stressProfile */
			inputfilestream >> _bulkPressureOutputTimesteps;				
		} else if (token == "confinementProperties") {	// 
			double wallThickness, horDist, vertDist, radius2, xmax, ymax, zmax;
			int cid;
			unsigned long upperID, lowerID;
			inputfilestream >> wallThickness >> horDist >> vertDist >> radius2 >> cid >> xmax >> ymax >> zmax >> upperID >> lowerID;
			cid--;
			_doRecordConfinement = true;
			_domain->setupConfinementProperties(wallThickness, horDist, vertDist, radius2, cid, xmax, ymax, zmax, upperID, lowerID);
		} else if (token == "confinementPropertiesFixedArea") {	// 
			double xmin, xmax, ymin, ymax, zmax;
			int cid;
			inputfilestream >> xmin >> xmax >> ymin >> ymax >> zmax >> cid;
			cid--;
			_doRecordConfinement = true;
			_domain->setupConfinementPropertiesFixed(xmin, xmax, ymin, ymax, zmax, cid);
		} else if (token == "confinementProfile") { /* TODO: suboption of confinementProperties */
			unsigned xun, yun;
			double correlationLength = 1.0;
			string stress, weightingFunc;
			inputfilestream >> xun >> yun >> stress;
			_domain->specifyStressCalcMethodConfinement(stress);
			if(stress == "Hardy") // calculates Hardy Stresses
			  _HardyConfinement = true;
			
			if(stress == "Hardy"){ // calculates Hardy Stresses
			  _HardyConfinement = true;
			  inputfilestream >> token;
			  if (token != "Linear" && token != "Pyramide" ) {
				global_log->error() << "Expected Hardy stress option 'Linear' or 'Pyramide' instead of '"
						<< token << "' for the weighting function.\n";
				exit(1);
			  }
			  weightingFunc = token;
			  _weightingConfinement = weightingFunc;
			  inputfilestream >> correlationLength;
			  if (correlationLength <= 0 || correlationLength >= 100) {
				global_log->error() << "Expected Hardy stress option for the correlation length in the range [0;100] instead of '"
						<< correlationLength << "'.\n";
				exit(1);
			  }
			}
			_domain->setupConfinementProfile(xun, yun, correlationLength); 
			global_log->info() << "Confinement stress profile: " << stress << endl;
		} else if (token == "confinementRecordingTimesteps") { /* TODO: suboption of confinementProfile */
			inputfilestream >> _confinementRecordingTimesteps; 	
			_domain->setupConfinementRecTime(_confinementRecordingTimesteps);
		} else if (token == "confinementOutputTimesteps") { /* TODO: suboption of confinementProfile */
			inputfilestream >> _confinementOutputTimesteps;
			_domain->setupConfinementOutTime(_confinementOutputTimesteps);
		} else if (token == "confinementProfilePrefix") { /* TODO: suboption of confinementProfile */
			inputfilestream >> _confinementProfileOutputPrefix;
		} else if (token == "confinementSplitStress") { /* TODO: suboption of confinementProfile */
			_doConfinementSplitStress = true;	
		} else if (token == "ThWallLayer") {
			double layerThickness;
			inputfilestream >> layerThickness;
			_domain->setThermostatWallLayer(layerThickness);
			// specify a thermostat for each layer in the components 'fixed' and 'moved'
			if( !_domain->severalThermostats() )
				_domain->enableLayerwiseWallThermostat();
		} else if (token == "shearRate") {	// applies shear by keeping the two mostouter layers of atoms fixed and setting a target velocity in the middle layer
			double xmin, xmax, ymin, ymax, shearRate, shearWidth;
			unsigned cid;
			inputfilestream >> xmin >> xmax >> ymin >> ymax >> cid >> shearRate >> shearWidth;
			cid--;
			_doShearRate = true;
			if(_doShearForce == true){
				global_log->error() << "There is already a shear force applied. Please, just select one of the two options!\n";
				exit(1);
			}
			_pressureGradient->setupShearRate(xmin, xmax, ymin, ymax, cid, shearRate, shearWidth, false);
		} else if (token == "shearForce") {	// applies shear by keeping the two mostouter layers of atoms fixed and setting a target force profile over the rest of the system
			double xmin, xmax, ymin, ymax, shearRate, shearWidth;
			unsigned cid;
			inputfilestream >> xmin >> xmax >> ymin >> ymax >> cid >> shearRate >> shearWidth;
			cid--;
			_doShearForce = true;
			if(_doShearRate == true){
				global_log->error() << "There is already a shear rate applied. Please, just select one of the two options!\n";
				exit(1);
			}
			_pressureGradient->setupShearRate(xmin, xmax, ymin, ymax, cid, shearRate, shearWidth, true);		
		} else if (token == "collectThermostatDirectedVelocity") { /* subotion of the thermostate replace with directe thermostate */
			inputfilestream >> _collectThermostatDirectedVelocity;
		} else if (token == "zOscillator") {
			_zoscillation = true;
			inputfilestream >> _zoscillator;
		}
		// chemicalPotential <mu> component <cid> [control <x0> <y0> <z0>
		// to <x1> <y1> <z1>] conduct <ntest> tests every <nstep> steps
		else if (token == "chemicalPotential") {
			double imu;
			inputfilestream >> imu;
			inputfilestream >> token;
			if (token != "component") {
				global_log->error() << "Expected 'component' instead of '"
						<< token << "'.\n";
				global_log->debug()
						<< "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
						<< "conduct <ntest> tests every <nstep> steps" << endl;
				exit(1);
			}
			unsigned icid;
			inputfilestream >> icid;
			icid--;
			inputfilestream >> token;
			double x0, y0, z0, x1, y1, z1;
			bool controlVolume = false;
			if (token == "control") {
				controlVolume = true;
				inputfilestream >> x0 >> y0 >> z0;
				inputfilestream >> token;
				if (token != "to") {
					global_log->error() << "Expected 'to' instead of '"
							<< token << "'.\n";
					global_log->debug()
							<< "Syntax: chemicalPotential <mu> component <cid> "
							<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
							<< "conduct <ntest> tests every <nstep> steps\n";
					exit(1);
				}
				inputfilestream >> x1 >> y1 >> z1;
				inputfilestream >> token;
			}
			if (token != "conduct") {
				global_log->error() << "Expected 'conduct' instead of '"
						<< token << "'.\n";
				global_log->debug()
						<< "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
						<< "conduct <ntest> tests every <nstep> steps\n";
				exit(1);
			}
			unsigned intest;
			inputfilestream >> intest;
			inputfilestream >> token;
			if (token != "tests") {
				global_log->error() << "Expected 'tests' instead of '" << token
						<< "'.\n";
				global_log->debug()
						<< "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
						<< "conduct <ntest> tests every <nstep> steps" << endl;
				exit(1);
			}
			inputfilestream >> token;
			if (token != "every") {
				global_log->error() << "Expected 'every' instead of '" << token
						<< "'.\n";
				global_log->debug()
						<< "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
						<< "conduct <ntest> tests every <nstep> steps" << endl;
				exit(1);
			}
			unsigned instep;
			inputfilestream >> instep;
			inputfilestream >> token;
			if (token != "steps") {
				global_log->error() << "Expected 'steps' instead of '" << token
						<< "'.\n";
				global_log->debug()
						<< "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
						<< "conduct <ntest> tests every <nstep> steps" << endl;
				exit(1);
			}
			ChemicalPotential tmu = ChemicalPotential();
			tmu.setMu(icid, imu);
			tmu.setInterval(instep);
			tmu.setOriginalInterval(instep);
			tmu.setInstances(intest);
			tmu.setOriginalInstances(intest);
			if (controlVolume)
				tmu.setControlVolume(x0, y0, z0, x1, y1, z1);
			global_log->info() << setprecision(6) << "chemical Potential "
					<< imu << " component " << icid + 1 << " (internally "
					<< icid << ") conduct " << intest << " tests every "
					<< instep << " steps: ";
			
			inputfilestream >> token;
			if (token == "barostat") {
			  double targetPress;
			  tmu.setGCMD_barostat(true);
			  inputfilestream >> targetPress;
			  tmu.setTargetPressure(targetPress);
			  tmu.setCountBarostating(0);
			  _domain->setupBarostat(tmu.getControl_bottom(0), tmu.getControl_top(0), tmu.getControl_bottom(1), tmu.getControl_top(1), tmu.getControl_bottom(2), tmu.getControl_top(2), tmu.getComponentID());
			}
			global_log->info() << flush;
			_lmu.push_back(tmu);
			global_log->info() << " pushed back." << endl;
		} else if (token == "planckConstant") {
			inputfilestream >> h;
		} else if(token == "Widom") {
                        widom = true;
		} else if(token == "cavity") {
                        unsigned cavity_cid;
                        inputfilestream >> cavity_cid;
                        cavity_cid--;
                        inputfilestream >> token;
                        if (token != "radius") {
                                global_log->error() << "Expected 'radius' instead of '"
                                                << token << "'.\n";
                                global_log->debug()
                                                << "Syntax: cavity <cid> radius <R> <coord_max> grid <N_x> <N_y> <N_z> every <n_step> steps\n";
                                exit(1);
                        }
                        double cavity_radius;
                        inputfilestream >> cavity_radius;
                        int neighbours;
                        inputfilestream >> neighbours;
                        inputfilestream >> token;
                        if (token != "grid") {
                                global_log->error() << "Expected 'grid' instead of '"
                                                << token << "'.\n";
                                global_log->debug()
                                                << "Syntax: cavity <cid> radius <R> <coord_max> grid <N_x> <N_y> <N_z> every <n_step> steps\n";
                                exit(1);
                        }
                        unsigned gridx, gridy, gridz;
                        inputfilestream >> gridx >> gridy >> gridz;
                        inputfilestream >> token;
                        if (token != "every") {
                                global_log->error() << "Expected 'every' instead of '"
                                                << token << "'.\n";
                                global_log->debug()
                                                << "Syntax: cavity <cid> radius <R> <coord_max> grid <N_x> <N_y> <N_z> every <n_step> steps\n";
                                exit(1);
                        }
                        unsigned cavity_steps;
                        inputfilestream >> cavity_steps;
                        inputfilestream >> token;
                        if (token != "steps") {
                                global_log->error() << "Expected 'steps' instead of '"
                                                << token << "'.\n";
                                global_log->debug()
                                                << "Syntax: cavity <cid> radius <R> <coord_max> grid <N_x> <N_y> <N_z> every <n_step> steps\n";
                                exit(1);
                        }
                        this->_mcav[cavity_cid] = CavityEnsemble();
                        this->_mcav[cavity_cid].setInterval(cavity_steps);
                        this->_mcav[cavity_cid].setMaxNeighbours(neighbours, cavity_radius*cavity_radius);
                        cavity_grid[cavity_cid] = new double[3];
                        (cavity_grid[cavity_cid])[0] = gridx;
                        (cavity_grid[cavity_cid])[1] = gridy;
                        (cavity_grid[cavity_cid])[2] = gridz;
		} else if (token == "NVE") {
			/* TODO: Documentation, what it does (no "Enerstat" at the moment) */
			_domain->thermostatOff();
			global_log->error() << "Not implemented" << endl;
			this->exit(1);
		} else if (token == "initCanonical") {
			inputfilestream >> _initCanonical;
		} else if (token == "initGrandCanonical") { /* suboption of chemical potential */
			inputfilestream >> _initGrandCanonical;
			inputfilestream >> token;
			if (token != "end") {
				global_log->error() << "Expected 'end' instead of '" << token
						<< "'.\n";
				exit(1);
			}
			inputfilestream >> _endGrandCanonical;
			global_log->info() << "Grand canonical ensemble in time slot: " << _initGrandCanonical << " until " << _endGrandCanonical << endl;
		} else if (token == "initStatistics") {
			inputfilestream >> _initStatistics;
		} else if (token == "cutoffRadius") {
			double rc;
			inputfilestream >> rc;
			this->setcutoffRadius(rc);
		} else if (token == "LJCutoffRadius") {
			double rc;
			inputfilestream >> rc;
			this->setLJCutoff(rc);
		} else if (token == "tersoffCutoffRadius") {
			double rc;
			inputfilestream >> rc;
			this->setTersoffCutoff(rc);
		} else if (token == "ultra") {
                        ultra = true;
		} else if (token == "cancelMomentum") { // cancels momentum of the total system
                        _cancelMomentum = true;
		} else {
			if (token != "")
				global_log->warning() << "Did not process unknown token "
						<< token << endl;
		}
	}
	
	// read particle data
	maxid = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain,
			_domainDecomposition);

	if (this->_LJCutoffRadius == 0.0)
		_LJCutoffRadius = this->_cutoffRadius;
	_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

	// @todo comment
	_integrator = new Leapfrog(timestepLength);

	// test new Decomposition
	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	unsigned idi = _lmu.size();
	unsigned j = 0;
	std::list<ChemicalPotential>::iterator cpit;
	for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
                if(widom) cpit->enableWidom();
		cpit->setIncrement(idi);
		double tmp_molecularMass = global_simulation->getEnsemble()->component(cpit->getComponentID())->m();
		cpit->setSystem(_domain->getGlobalLength(0),
				_domain->getGlobalLength(1), _domain->getGlobalLength(2),
				tmp_molecularMass);
		cpit->setGlobalN(global_simulation->getEnsemble()->component(cpit->getComponentID())->getNumMolecules());
		cpit->setNextID(j + (int) (1.001 * (256 + maxid)));

		cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0),
				_moleculeContainer->getBoundingBoxMax(0),
				_moleculeContainer->getBoundingBoxMin(1),
				_moleculeContainer->getBoundingBoxMax(1),
				_moleculeContainer->getBoundingBoxMin(2),
				_moleculeContainer->getBoundingBoxMax(2));
		/* TODO: thermostat */
		double Tcur = _domain->getCurrentTemperature(0);
		/* FIXME: target temperature from thermostat ID 0 or 1?  */
		double
				Ttar =_domain->severalThermostats() ? _domain->getTargetTemperature(1)
								: _domain->getTargetTemperature(0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;
		cpit->submitTemperature(Tcur);
		if (h != 0.0)
			cpit->setPlanckConstant(h);

		j++;
	}
	
 	unsigned long Nbasis = 3 * (this->_domain->getglobalNumMolecules() + 3 * this->_lmu.size());
	map<unsigned, CavityEnsemble>::iterator ceit;
        for(ceit = _mcav.begin(); ceit != _mcav.end(); ceit++)
        {
           if(ultra) ceit->second.enableUltra();
           ceit->second.setSystem(_domain->getGlobalLength(0), _domain->getGlobalLength(1), _domain->getGlobalLength(2));
           ceit->second.setIdOffset((1 + ceit->first) * Nbasis);

           ceit->second.setSubdomain( ownrank, _moleculeContainer->getBoundingBoxMin(0),
                                               _moleculeContainer->getBoundingBoxMax(0),
                                               _moleculeContainer->getBoundingBoxMin(1),
                                               _moleculeContainer->getBoundingBoxMax(1),
                                               _moleculeContainer->getBoundingBoxMin(2),
                                               _moleculeContainer->getBoundingBoxMax(2)  );
           double Tcur = _domain->getCurrentTemperature(0);
           double Ttar =_domain->severalThermostats()? _domain->getTargetTemperature(1)
                                                     : _domain->getTargetTemperature(0);
           if ((Tcur < 0.667 * Ttar) || (Tcur > 1.5 * Ttar)) Tcur = Ttar;
           ceit->second.submitTemperature(Tcur);
           
           ceit->second.init(
              global_simulation->getEnsemble()->component(ceit->first),
              (cavity_grid[ceit->first])[0], (cavity_grid[ceit->first])[1], (cavity_grid[ceit->first])[2]
           );
           ceit->second.communicateNumCavities(this->_domainDecomposition);
        }
}

void Simulation::prepare_start() {
	global_log->info() << "Initializing simulation" << endl;

	global_log->info() << "Initialising cell processor" << endl;
#if ENABLE_VECTORIZED_CODE
	global_log->debug() << "Checking if vectorized cell processor can be used" << endl;
	bool lj_present = false;
	bool charge_present = false;
	bool dipole_present = false;
	bool quadrupole_present = false;
	bool tersoff_present = false;

	const vector<Component> components = *(global_simulation->getEnsemble()->components());
	for (size_t i = 0; i < components.size(); i++) {
		lj_present |= components[i].numLJcenters() != 0;
		charge_present |= components[i].numCharges() != 0;
		dipole_present |= components[i].numDipoles() != 0;
		quadrupole_present |= components[i].numQuadrupoles() != 0;
		tersoff_present |= components[i].numTersoff() != 0;
	}

	if (charge_present || dipole_present || quadrupole_present || tersoff_present) {
		global_log->warning() << "Using legacy cell processor. (Vectorized code not yet available for charges, dipoles, quadrupoles and tersoff interactions.)" << endl;
		global_log->debug() << "xx lj present: " << lj_present << endl;
		global_log->debug() << "xx charge present: " << charge_present << endl;
		global_log->debug() << "xx dipole present: " << dipole_present << endl;
		global_log->debug() << "xx quadrupole present: " << quadrupole_present << endl;
		global_log->debug() << "xx tersoff present: " << tersoff_present << endl;
		_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius, _particlePairsHandler);
	} else {
		global_log->info() << "Using vectorized cell processor." << endl;
		_cellProcessor = new VectorizedCellProcessor( *_domain,_LJCutoffRadius);
	}
#else
	global_log->info() << "Using legacy cell processor." << endl;
	_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius, _particlePairsHandler);
#endif

	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();
	global_log->info() << "Updating domain decomposition" << endl;
	updateParticleContainerAndDecomposition();
	global_log->info() << "Performing inital force calculation" << endl;
	
	_moleculeContainer->traverseCells(*_cellProcessor);

	/* If enabled count FLOP rate of LS1. */
	if( NULL != _ljFlopCounter ) {
		_moleculeContainer->traverseCells(*_ljFlopCounter);
	}

	// clear halo
	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();

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
	
	// once at the beginning of the calculation the target velocity of accelerateInstantaneously is transferred (from to pressureGradient to Domain)
	if (_domain->getPG()->isAcceleratingInstantaneously(_domain->getNumberOfComponents()))
		_domain->setTargetVelocityAcceleration();
		
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

	if (_zoscillation) {
		global_log->debug() << "Initializing z-oscillators" << endl;
		_integrator->init1D(_zoscillator, _moleculeContainer);
	}

	// initialize output
	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter
			!= _outputPlugins.end(); outputIter++) {
		(*outputIter)->initOutput(_moleculeContainer, _domainDecomposition,
				_domain);
	}
	
	// Initialize spring force  
	for(_thismol = _moleculeContainer->begin(); _thismol != _moleculeContainer->end(); _thismol = _moleculeContainer->next())
	    for(int d = 0; d < 3; d++)
		_thismol->setF_Spring(d, 0.0);


	global_log->info() << "System initialised\n" << endl;
	global_log->info() << "System contains "
			<< _domain->getglobalNumMolecules() << " molecules." << endl;
}

void Simulation::simulate() {

	global_log->info() << "Started simulation" << endl;
	
	// (universal) constant acceleration (number of) timesteps
	unsigned uCAT = _pressureGradient->getUCAT();
	double Force[3][3], SpringForce[3][3];
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j< 3; j++){
	    Force[i][j] = 0.0;
	    SpringForce[i][j] = 0.0;
	  }
	  
	// reset of heat flux average
	for (int tht_id = 0; tht_id <= _domain->maxThermostat(); tht_id++)
		_domain->setThT_heatFlux(tht_id, 0.0);
	  
	// initialize result files
	    std::ofstream ForceData;
	    std::ofstream BulkPressure;
	    std::ofstream Confinement;
	    std::ofstream HeatFlux;
	    std::ofstream Energy;
	    
	// writing header for result files    
	if(_domainDecomposition->getRank() == 0){
		if(ForceData.is_open()) ForceData.close();
		if(BulkPressure.is_open()) BulkPressure.close();
		if(Confinement.is_open()) Confinement.close();
		if(HeatFlux.is_open()) HeatFlux.close();
		if(Energy.is_open()) Energy.close();
	  
		ForceData.open("./Results/ForceData.dat");
		if (!ForceData.is_open()) {
		global_log->error() << "Could not open file " << "./Results/ForceData.dat" << endl;
		exit(1);
		}
		ForceData << "#step\t\t\tForce_x Force_y Force_z\t\tSpringForce_x SpringForce_y SpringForce_z\t\t";
		ForceData << "FrictionForce(x)\t\t";
		ForceData << "Before Acc: VelAverage_x VelAverage_y VelAverage_z\t\t";
		ForceData << "After Acc: VelAverage_x VelAverage_y VelAverage_z\t\tAfter ThT: VelAverage_x VelAverage_y VelAverage_z\n";
		
		BulkPressure.open("./Results/BulkPressure.dat");
		if (!BulkPressure.is_open()) {
		global_log->error() << "Could not open file " << "./Results/BulkPressure.dat" << endl;
		exit(1);
		}
		BulkPressure << "#step\t\t\tBulk Pressure\t\t\tDensity\t\t\tTemperature\n";
		
		Confinement.open("./Results/Confinement.dat");
		if (!Confinement.is_open()) {
		global_log->error() << "Could not open file " << "./Results/Confinement.dat" << endl;
		exit(1);
		}
		Confinement << "#step\t\t\tNumberOfMolecules\tFluidVolume\tSolidSurface\tdMidpoints\tdAverage\t";
		Confinement << "Density\t\tPressure\tTemperature1\tTemperature2\t\t\tForce_x\tForce_y\tEffectiveViscosity\n";
		
		HeatFlux.open("./Results/HeatFlux.dat");
		if (!HeatFlux.is_open()) {
		global_log->error() << "Could not open file " << "./Results/HeatFlux.dat" << endl;
		exit(1);
		}
		HeatFlux << "#A_Comp_2 (upper plate) " << 2*_domain->getGlobalLength(0)*_domain->getGlobalLength(2) << "\n" << "#A_Comp_3 (lower plate)" << 2*((_domain->getConfinementEdge(1)-_domain->getConfinementEdge(0)) + _domain->get_confinementMidPoint(3))*_domain->getGlobalLength(2) << "\n\n";
		if(_domain->isThermostatWallLayer()){
			for(int tht = 1; tht <= _domain->maxThermostat(); tht++) 
				HeatFlux << "#Layer " << _domain->getThermostat(tht) << "\n";
			HeatFlux << "\n";
		}else{
			for(unsigned cid = 0; cid < _domain->getNumberOfComponents(); cid++) 
				HeatFlux << "#Comp " << cid << " " << _domain->getThermostat(cid) << "\n";
			HeatFlux << "\n";
		}
		HeatFlux << "#step\t\t\tQ_ThT_ges\tQ_ThT_1\tQ_ThT_2\tQ_ThT_3";
		if(_doShearRate)
			HeatFlux << "\tdEnergy_noDirVel\tdEnergy_Default\tshearVelDelta";
		HeatFlux << "\n\n";
		
		Energy.open("./Results/Energy.dat");
		if (!Energy.is_open()) {
		global_log->error() << "Could not open file " << "./Results/Energy.dat" << endl;
		exit(1);
		}
		
		Energy << "#step\t\t\t";
		for(unsigned cid = 0; cid < _domain->getNumberOfComponents(); cid++) 
			Energy << "Upot[cid]\t\tUkin[cid]\t\t";
		Energy << "\n";	
		//Check discjunction of thermostat layers
		if(_domain->isThermostatLayer() == true)
		  _domain->checkThermostatLayerDisjunct();
		cout << "LOST" << endl;
	}
	
		
	_initSimulation = (unsigned long) (getSimulationTime() / _integrator->getTimestepLength());

	/* demonstration for the usage of the new ensemble class */
	CanonicalEnsemble ensemble(_moleculeContainer, global_simulation->getEnsemble()->components());
	ensemble.updateGlobalVariable(NUM_PARTICLES);
	global_log->debug() << "Number of particles in the Ensemble: "
			<< ensemble.N() << endl;
	ensemble.updateGlobalVariable(ENERGY);
	global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E()
		<< endl;
	
	ensemble.updateGlobalVariable(TEMPERATURE);
	global_log->debug() << "Temperature of the Ensemble: " << ensemble.T()
		<< endl;

	/***************************************************************************/
	/* BEGIN MAIN LOOP                                                         */
	/***************************************************************************/
	// all timers except the ioTimer measure inside the main loop
	Timer loopTimer;
	Timer decompositionTimer;
	Timer perStepIoTimer;
	Timer ioTimer;

	loopTimer.start();
	for (_simstep = _initSimulation; _simstep <= _numberOfTimesteps; _simstep++) {
		if (_simstep >= _initGrandCanonical && _simstep <= _endGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
				if ((!((_simstep + 2 * j + 3) % cpit->getInterval()) && _domain->getDifferentBarostatInterval() == false) || (_simstep == cpit->getInterval() && _domain->getDifferentBarostatInterval() == true)) {
					cpit->prepareTimestep(_moleculeContainer,
							_domainDecomposition);
				}
				j++;
			}
		}
		if (_simstep >= _initStatistics) {
                   map<unsigned, CavityEnsemble>::iterator ceit;
                   for(ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++)
                   {
                      if (!((_simstep + 2 * ceit->first + 3) % ceit->second.getInterval()))
                      {
                         ceit->second.preprocessStep();
                      }
                   }
                }
		
		// initialize Hardy stresses
		  if(_HardyStress && _simstep >= this->_initStatistics)
		    for(_thismol = _moleculeContainer->begin(); _thismol != _moleculeContainer->end(); _thismol = _moleculeContainer->next()){
		      _thismol->setHardyStress(true);
		      _thismol->setHardyStressWeighting(_weightingStress);
		    }
		  if(_HardyConfinement && _simstep >= this->_initStatistics)
		    for(_thismol = _moleculeContainer->begin(); _thismol != _moleculeContainer->end(); _thismol = _moleculeContainer->next()){
		      _thismol->setHardyConfinement(true);
		      _thismol->setHardyConfinementWeighting(_weightingConfinement);
		    }
		  
		
		global_log->debug() << "timestep: " << getSimulationStep() << endl;
		global_log->debug() << "simulation time: " << getSimulationTime() << endl;

		_integrator->eventNewTimestep(_moleculeContainer, _domain);
		// activate RDF sampling
		if ((_simstep >= this->_initStatistics) && this->_rdf != NULL) {
			this->_rdf->tickRDF();
			this->_particlePairsHandler->setRDF(_rdf);
			this->_rdf->accumulateNumberOfMolecules(*(global_simulation->getEnsemble()->components()));
		}
		// ensure that all Particles are in the right cells and exchange Particles
		global_log->debug() << "Updating container and decomposition" << endl;
		loopTimer.stop();
		decompositionTimer.start();
                if (_simstep >= _initStatistics) {
                   map<unsigned, CavityEnsemble>::iterator ceit;
                   for(ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++)
                   {
                      if( (!((_simstep + 2 * ceit->first + 4) % ceit->second.getInterval())) ||
                          (!((_simstep + 2 * ceit->first + 2) % ceit->second.getInterval())) ||
                          (!((_simstep + 2 * ceit->first) % ceit->second.getInterval())) )
                      {
                         _domainDecomposition->exchangeCavities(&ceit->second, _domain);
                      }
                   }
                }
                
		updateParticleContainerAndDecomposition();
		decompositionTimer.stop();
		loopTimer.start();

		// Force calculation
		global_log->debug() << "Traversing pairs" << endl;
		_moleculeContainer->traverseCells(*_cellProcessor);
		
		//Force Record of spring force and total force
		unsigned long timeStepForce = 1000;
		if(_simstep%timeStepForce == 0 /*&& _simstep >= _initStatistics*/){
		    string moved ("moved");
		    unsigned cid = _domain->getCidMovement(moved, _domain->getNumberOfComponents());
		    cid--;
		    if (cid >= 0 && cid < _domain->getNumberOfComponents()){
		      if(cid >= 0 && cid < (unsigned)_domain->getNumberOfComponents())
		       _pressureGradient->collectForcesOnComponent(_domainDecomposition, cid);
		      if(_simstep-_initSimulation == timeStepForce && _domain->isSpringDamped())
			_pressureGradient->resetSpringForcesOnComponent(cid);
		      
		      if(_simstep > _initStatistics){
			for(int d = 0; d < 3; d++)
			      Force[d][cid] = _pressureGradient->getGlobalForceSum(d, cid)/timeStepForce;
			if(_domain->isSpringDamped()){
			      _pressureGradient->collectSpringForcesOnComponent(_domainDecomposition, cid);

			    if(_domainDecomposition->getRank() == 0)   
			      for(int d = 0; d < 3; d++)
			     	      SpringForce[d][cid] = _pressureGradient->getGlobalSpringForceSum(d, cid)/timeStepForce;
			}
			
			if(_domainDecomposition->getRank()==0){  
			  ForceData << _simstep << "\t\t\t" << Force[0][cid] << " " << Force[1][cid] << " " << Force[2][cid];
			  if(_domain->isSpringDamped())
			    ForceData << "\t\t\t" << SpringForce[0][cid] << " " << SpringForce[1][cid] << " " << SpringForce[2][cid];
			}
			if(cid >= 0 && cid < (unsigned)_domain->getNumberOfComponents())
			  _pressureGradient->resetForcesOnComponent(cid);
			if(_domain->isSpringDamped())
			  _pressureGradient->resetSpringForcesOnComponent(cid);
		      }
		    }
		}
		
		// test deletions and insertions for the grand canonical ensemble and the barostat
		if (_simstep >= _initGrandCanonical && _simstep <= _endGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
				if ((!((_simstep + 2 * j + 3) % cpit->getInterval()) && _domain->getDifferentBarostatInterval() == false) || (_simstep == cpit->getInterval() && _domain->getDifferentBarostatInterval() == true)) {
					double currentTemp = _domain->getGlobalCurrentTemperature();
					double currentPressure = 1.0;
					double targetPressure = 0.0;
					double dampFac = 1.0;
					_domain->setDifferentBarostatInterval(true);
					if(cpit->isGCMD_barostat()){
					    double imu;
					    unsigned icid = cpit->getComponentID();
					    long double p_Vir = _domain->getPressureVirial_barostat();
					    long double p_Kin = _domain->getPressureKin_barostat();
					    double Volume = cpit->getVolume_Barostat();
					    if (_doRecordConfinement){
					      Volume = _domain->get_confinedVolume();
					    }
					    unsigned long N = _domain->getN_barostat();
					    targetPressure = cpit->getTargetPressure();
					    currentPressure = (double)((p_Vir + p_Kin)/(3 * Volume));
					    currentTemp = (double)(p_Kin/(3 * N));
					    // artificial prefactor (~damping function) applied to the number number of test insertions and time interval between test
					    dampFac = ((long double)(_endGrandCanonical - _simstep))/((long double)(_endGrandCanonical - _initGrandCanonical)); 
					    dampFac = 0.1 + dampFac * dampFac;
					    
					    if(_domainDecomposition->getRank()==0)
					      cout << " t " << _simstep << " targetPressure " << targetPressure << " currentPressure " << currentPressure << " currentTemp " << currentTemp 
					      << " Volume " << Volume << " N " << N << " damping " << dampFac;
					    
					    double deltaPressure = abs(currentPressure-targetPressure);
					    if(deltaPressure == 0)
					      deltaPressure = 0.01;
					    unsigned newInstances = (unsigned)floor(cpit->getOriginalInstances()*10*deltaPressure*dampFac);
					    // amount of instances are reglemented to be not bigger than 3 times the amount of original instances
					    if(newInstances > 3 * cpit->getOriginalInstances())
					      newInstances = 3 * cpit->getOriginalInstances();
					    
					    if(_domainDecomposition->getRank()==0)
					      cout << " newInst " << newInstances << endl;
					  if(0.99*targetPressure > currentPressure){
					    imu = 50;
					    cpit->setMu(icid, imu);
					    cpit->setInstances(newInstances);
					  }else if(1.01*targetPressure < currentPressure){
					    imu = -5;
					    cpit->setMu(icid, imu);
					    cpit->setInstances(newInstances);
					  }else{
					    imu = -0.1;
					    cpit->setMu(icid, imu);
					    cpit->setInstances(5);
					  }
					  // reset barostat specific state variables
					  _domain->resetBarostat();
					}
					global_log->debug() << "Grand canonical ensemble(" << j
							<< "): test deletions and insertions" << endl;
                                        this->_domain->setLambda(cpit->getLambda());
                                        this->_domain->setDensityCoefficient(cpit->getDensityCoefficient());
					_moleculeContainer->grandcanonicalStep(&(*cpit),
							currentTemp, this->_domain, *_cellProcessor);
#ifndef NDEBUG
					/* check if random numbers inserted are the same for all processes... */
					cpit->assertSynchronization(_domainDecomposition);
#endif
					int localBalance =
							_moleculeContainer->localGrandcanonicalBalance();
					int balance = _moleculeContainer->grandcanonicalBalance(
							_domainDecomposition);
					global_log->debug() << "   b["
							<< ((balance > 0) ? "+" : "") << balance << "("
							<< ((localBalance > 0) ? "+" : "") << localBalance
							<< ")" << " / c = " << cpit->getComponentID()
							<< "]   " << endl;
					_domain->Nadd(cpit->getComponentID(), balance, localBalance);
					
					// adaptable timestep for barostat: the bigger the difference between targetPressure and currentPressure, the shorter the timestep
					if(cpit->isGCMD_barostat()){
					  cpit->addCountBarostating(1);
					  double deltaPressure = abs(currentPressure-targetPressure);
					  
					  if(deltaPressure == 0)
					    deltaPressure = 0.01;
					  unsigned newInterval = (unsigned)floor(cpit->getOriginalInterval()/deltaPressure/dampFac);
					  if(_doRecordConfinement && _confinementRecordingTimesteps <= 10)
					    newInterval = newInterval * 0.5 * _confinementRecordingTimesteps;
					  else if(_doRecordConfinement && _confinementRecordingTimesteps > 10)
					    newInterval = newInterval * 0.5 * 10;
					  
					  // after 7 times barostating a break of 50000 timesteps is forced for equilibration
					  if(cpit->getCountBarostating()%7 == 0)
					    newInterval = 50000;
					  
					  if(newInterval > 50000)
					    newInterval = 50000;
					  if(newInterval < 5 * _confinementRecordingTimesteps)
					    newInterval = 5 * _confinementRecordingTimesteps;
					  newInterval += _simstep;
					  
					  cpit->setInterval(newInterval);
					  if(_domainDecomposition->getRank()==0)
					    cout << " newI " << newInterval << endl;
					}
				}
				j++;
			}
		}
		
		if(_simstep >= _initStatistics)
                {
                   map<unsigned, CavityEnsemble>::iterator ceit;
                   for(ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++)
                   {
                      if (!((_simstep + 2 * ceit->first + 3) % ceit->second.getInterval()))
                      {
                         global_log->debug() << "Cavity ensemble for component " << ceit->first << ".\n";
                         
                         this->_moleculeContainer->cavityStep(
                            &ceit->second, _domain->getGlobalCurrentTemperature(), this->_domain, *_cellProcessor
                         );
                      }
                           
                      if( (!((_simstep + 2 * ceit->first + 7) % ceit->second.getInterval())) ||
                          (!((_simstep + 2 * ceit->first + 3) % ceit->second.getInterval())) ||
                          (!((_simstep + 2 * ceit->first - 1) % ceit->second.getInterval())) )
                      {                                   
                         this->_moleculeContainer->numCavities(&ceit->second, this->_domainDecomposition);
                      }
                   }
                }
                
		global_log->debug() << "Delete outer particles / clearing halo" << endl;
		_moleculeContainer->deleteOuterParticles();
		
		// Collect forces in the confinement for the calculation of the viscosity
		if ((_simstep > _initStatistics) && _doRecordConfinement)
		    _domain->collectForcesOnComponentConfinement(_moleculeContainer);
		
		if ((_simstep > _initStatistics)) {
		  // collects ID of component that moves due to an external force field
		   string moved ("moved");
		   unsigned cid = _domain->getCidMovement(moved, _domain->getNumberOfComponents());
		   cid--;
		   if(cid >= 0 && cid < _domain->getNumberOfComponents())
		       _pressureGradient->calculateForcesOnComponent(_moleculeContainer, cid);
		   if(_domain->isSpringDamped())
		       _pressureGradient->calculateSpringForcesOnComponent(_moleculeContainer, cid);
		}
		
		if (_simstep >= _initGrandCanonical && _simstep <= _endGrandCanonical) {
			_domain->evaluateRho(_moleculeContainer->getNumberOfParticles(),
					_domainDecomposition);
		}
		if (!(_simstep % _collectThermostatDirectedVelocity))
			_domain->calculateThermostatDirectedVelocity(_moleculeContainer);
		/*
		 * radial distribution function
		 */
		if (_simstep >= _initStatistics) {
			if (this->_lmu.size() == 0) {
				this->_domain->record_cv();
			}
		}
		
		/*
                 * cavity (cluster) analysis
                 */
                if (_simstep >= _initStatistics) {
                   map<unsigned, CavityEnsemble>::iterator ceit;
                   for(ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++)
                   {
                      if( (!((_simstep + 2 * ceit->first + 5) % ceit->second.getInterval())) ||
                          (!((_simstep + 2 * ceit->first + 2) % ceit->second.getInterval())) ||
                          (!((_simstep + 2 * ceit->first - 1) % ceit->second.getInterval())) )
                      {
                         ceit->second.detectClusters();
                      }
                   }
                }

		if (_zoscillation) {
			global_log->debug() << "alert z-oscillators" << endl;
			_integrator->zOscillation(_zoscillator, _moleculeContainer);
		}
		// Inform the integrator about the calculated forces
		global_log->debug() << "Inform the integrator" << endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);
		
		if(_domainDecomposition->getRank() == 0 && (_simstep % 10) == 0) 
			for (int t = 0; t <= _domain->maxThermostat(); t++){
				cout << " T" << t << ": " << _domain->getCurrentTemperature(_domain->getThermostat(t));
				if (t == _domain->maxThermostat())
					cout << endl;
			}
		// Acceleration of the component type "moved"; independent of the shear rate calculation
		if (_pressureGradient->isAcceleratingUniformly()) {
			if (!(_simstep % uCAT)) {
				global_log->debug() << "Determine the additional acceleration"
						<< endl;
				_pressureGradient->determineAdditionalAcceleration(
						_domainDecomposition, _moleculeContainer, uCAT
						* _integrator->getTimestepLength());
			}
			global_log->debug() << "Process the uniform acceleration" << endl;
			_integrator->accelerateUniformly(_moleculeContainer, _domain);
			_pressureGradient->adjustTau(this->_integrator->getTimestepLength());
		}else if (_pressureGradient->isAcceleratingInstantaneously(_domain->getNumberOfComponents())){
		  if (_simstep >= _initStatistics){
		    	global_log->debug() << "Determine the instantaneous acceleration" << endl;
			_integrator->accelerateInstantaneously(_domainDecomposition, _moleculeContainer, _domain);
			// collects ID of component that is moving due to external force field
			string moved ("moved");
			unsigned cid = _domain->getCidMovement(moved, _domain->getNumberOfComponents());
			cid--;
			// initialization for time averaging of velocity sums of moved component
			if(_simstep == _initStatistics){
			  for(int d = 0; d < 3; d++){
			      _pressureGradient->setGlobalVelSumBeforeAcc(d, cid, 0);
			      _pressureGradient->setGlobalVelSumAfterAcc(d, cid, 0);
			  }
			}
			if(cid >= 0 && cid < _domain->getNumberOfComponents() && _domainDecomposition->getRank()==0 && _simstep%timeStepForce == 0){
			    double tmp_molecularMass = global_simulation->getEnsemble()->component(cid)->m();
			    ForceData << "\t\t\t\t" << (_pressureGradient->getGlobalVelSumAfterAcc(0, cid)-_pressureGradient->getGlobalVelSumBeforeAcc(0, cid))/(timeStepForce*getTimeStepLength())*_pressureGradient->getGlobalN(cid)*tmp_molecularMass;
			    ForceData << "\t\t\t\t" << _pressureGradient->getGlobalVelSumBeforeAcc(0, cid)/timeStepForce << " " << _pressureGradient->getGlobalVelSumBeforeAcc(1, cid)/timeStepForce << " " << _pressureGradient->getGlobalVelSumBeforeAcc(2, cid)/timeStepForce;
			    ForceData << "\t\t\t" << _pressureGradient->getGlobalVelSumAfterAcc(0, cid)/timeStepForce << " " << _pressureGradient->getGlobalVelSumAfterAcc(1, cid)/timeStepForce << " " << _pressureGradient->getGlobalVelSumAfterAcc(2, cid)/timeStepForce;
			    for(int d = 0; d < 3; d++){
			      _pressureGradient->setGlobalVelSumBeforeAcc(d, cid, 0);
			      _pressureGradient->setGlobalVelSumAfterAcc(d, cid, 0);
			    }
			}
		  }
		}
				
		// shear rate accelerates the total system with a velocity gradient; for _doShearForce just the margin layers are fixed
		if((_doShearRate||_doShearForce) && _simstep >= _initStatistics)
			_integrator->shearRate(_domainDecomposition, _moleculeContainer, _domain);
		// calculate the global macroscopic values from the local values
		global_log->debug() << "Calculate macroscopic values" << endl;
		_domain->calculateGlobalValues(_domainDecomposition,
				_moleculeContainer, (!(_simstep % _collectThermostatDirectedVelocity)), Tfactor(
								_simstep));
		if(_domainDecomposition->getRank()==0 && _doEnergyOutputPerComponent == true){
			for(unsigned comp = 0; comp < _domain->getNumberOfComponents(); comp++)
			Energy << " Upot[" << comp << "] = " << _domain->getGlobalUpot(comp) << " Ukin[" << comp << "] = " << _domain->getGlobalUkin(comp);
			Energy << endl;
		}
		// scale velocity and angular momentum
		if (!_domain->NVE()) {
			global_log->debug() << "Velocity scaling" << endl;
			if (_domain->severalThermostats()) {
				_velocityScalingThermostat.enableComponentwise();
				for(int thermostatId = 1; thermostatId <= _domain->maxThermostat(); thermostatId++) {
					unsigned int cid = _domain->getCidToThermostat(thermostatId);
					if(_domain->isScaling1Dim(thermostatId) && _domain->getAlphaTransCorrection(thermostatId) == false){
					      _velocityScalingThermostat.setAlphaTrans(thermostatId, _domain->getGlobalAlphaTrans(thermostatId));
					      _velocityScalingThermostat.setBetaTrans(thermostatId, _domain->getGlobalBetaTrans(thermostatId));
					      _velocityScalingThermostat.setBetaRot(thermostatId, _domain->getGlobalBetaRot(thermostatId));
					}
					else{
					    _velocityScalingThermostat.setBetaTrans(thermostatId, _domain->getGlobalBetaTrans(thermostatId));
					    _velocityScalingThermostat.setBetaRot(thermostatId, _domain->getGlobalBetaRot(thermostatId));
					}
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
			_velocityScalingThermostat.apply(_moleculeContainer, _domainDecomposition);
			// Control of velocity after thermostat scaling
			if(_simstep >= _initStatistics){
			    // collects ID of the component that is allowed to move
			    string moved ("moved");
			    unsigned cid = _domain->getCidMovement(moved, _domain->getNumberOfComponents());
			    cid--;
			    if(_simstep == _initStatistics){
			      for(int d = 0; d < 3; d++)
				  _pressureGradient->setGlobalVelSumAfterThT(d, cid, 0);
			    }
			    if(cid >= 0 && cid < _domain->getNumberOfComponents()){  
			      _pressureGradient->prepare_getMissingVelocity(_domainDecomposition, _moleculeContainer, cid, _domain->getNumberOfComponents(), getDirectedVelocityTime());
			      for(int d = 0; d < 3; d++)
				_pressureGradient->addGlobalVelSumAfterThT(d, cid, _pressureGradient->getGlobalVelSum(d,cid) / _pressureGradient->getGlobalN(cid));
			
			      if(_domainDecomposition->getRank()==0 && _simstep >= _initStatistics && _simstep%timeStepForce == 0){
				ForceData << "\t\t\t\t" << _pressureGradient->getGlobalVelSumAfterThT(0, cid)/timeStepForce << " " << _pressureGradient->getGlobalVelSumAfterThT(1, cid)/timeStepForce << " " << _pressureGradient->getGlobalVelSumAfterThT(2, cid)/timeStepForce << endl;
				for(int d = 0; d < 3; d++)
				    _pressureGradient->setGlobalVelSumAfterThT(d, cid, 0);
			      }
			      
			    }
			}
		}
		// Calculate directed velocities
		if(_boolDirectedVel == true){
			if(_doSimpleAverage == true)
				_velocityScalingThermostat.calculateDirectedVelocitiesSimple(_moleculeContainer, _domainDecomposition);
			else if(_doSimpleAverageSigma2D == true)
				_velocityScalingThermostat.calculateDirectedVelocitiesSimple(_moleculeContainer, _domainDecomposition);
			else if(_doSimpleAverageSigma3D == true)
				_velocityScalingThermostat.calculateDirectedVelocitiesSimpleSigma3D(_moleculeContainer, _domainDecomposition);
			else if(_doNeighbourAverage == true)
				_velocityScalingThermostat.calculateDirectedVelocitiesNeighbourList(_moleculeContainer, _domainDecomposition);
			else if(_doMovingAverageSigma2D == true)
				_velocityScalingThermostat.calculateDirectedVelocities(_moleculeContainer, _domainDecomposition);
			else
				_velocityScalingThermostat.calculateDirectedVelocities(_moleculeContainer, _domainDecomposition);
		}

		advanceSimulationTime(_integrator->getTimestepLength());

		/* BEGIN PHYSICAL SECTION:
		 * the system is in a consistent state so we can extract global variables
		 */
		ensemble.updateGlobalVariable(NUM_PARTICLES);
		global_log->debug() << "Number of particles in the Ensemble: "
				<< ensemble.N() << endl;
		ensemble.updateGlobalVariable(ENERGY);
		global_log->debug() << "Kinetic energy in the Ensemble: "
				<< ensemble.E() << endl;
		ensemble.updateGlobalVariable(TEMPERATURE);
		global_log->debug() << "Temperature of the Ensemble: " << ensemble.T()
			<< endl;
		/* END PHYSICAL SECTION */

		// measure per timestep IO
		loopTimer.stop();
		perStepIoTimer.start();
		output(_simstep);

		// Writing confinement Properties as output
		if((_simstep > _initStatistics) && _doRecordConfinement && !(_simstep % _confinementOutputTimesteps)){
		  _domain->collectConfinementProperties(_domainDecomposition);
		  if(_domainDecomposition->getRank()==0){
		    unsigned long N = _domain->get_universalNConfinement();
		    long double p_Vir = _domain->get_universalPressureVirial_Confinement();
		    long double p_Kin = _domain->get_universalPressureKin_Confinement();
		    long double dof = _domain->get_universalDOFProfile_Confinement();
		    long double eKin = _domain->get_universalKineticProfile_Confinement();
		    std::map<int,long double> force;
		    for (int d = 0; d < 3; d++)
		      force[d] = _domain->get_globalForce_Confinement(d);
		    double dMax = (double)_domain->get_dMax();
		    double dAverage = _domain->get_dAveraged();
		    double Volume = _domain->get_confinedVolume();
		    double Surface = _domain->get_confinedSolidSurface();
		    double pressure = (double)((p_Vir + p_Kin)/(3 * Volume));
		    double density = (double)(N/Volume);
		    double temperature1 = (double)(p_Kin/(3 * N));
		    double temperature2 = (double)(eKin/dof);
		    double velocity_x_solid =  _pressureGradient->getGlobalTargetVelocity(0,1);
		    double viscosity = (double)((force[0]/Surface)/(velocity_x_solid/dAverage));

		    Confinement << _simstep << "\t\t\t" << N << "\t\t\t" << Volume << "\t\t" << Surface << "\t\t" << dMax << "\t\t" << dAverage << "\t\t";
		    Confinement << density << "\t\t" << pressure << "\t" << temperature1 << "\t" << temperature2 << "\t\t\t";
		    Confinement << force[0] << "\t" << force[1] << "\t" << viscosity << "\n";
		    
		    ostringstream osstrm;
		    osstrm << _confinementProfileOutputPrefix << ".";
		    osstrm.fill('0');
		    osstrm.width(9);
		    osstrm << right << _simstep;
		    _domain->outputConfinementProperties(osstrm.str().c_str(), _pressureGradient);
		    osstrm.str("");
		    osstrm.clear();
		    
		    // Heat flux consumed by the thermostats
		    HeatFlux << _simstep << "\t";
		    for (int tht_id = 0; tht_id <= _domain->maxThermostat(); tht_id++)
			HeatFlux << _domain->getThT_heatFlux(tht_id) << "\t";
		    if(_doShearRate){
			    HeatFlux << "\t" << 0.5 * (_domain->getPostShearEnergyConf() - _domain->getPreShearEnergyConf()) << "\t" << 0.5 * (_domain->getPostShearEnergyDefault() - _domain->getPreShearEnergyDefault()) << "\t" << _domain->getShearVelDelta()/((double)_domain->getglobalNumMolecules()*(double)_confinementOutputTimesteps);
		    }
		    HeatFlux << "\n";
		  }		
		    _domain->resetConfinementProperties(); 
		    // reset heat flux average
		    for (int tht_id = 0; tht_id <= _domain->maxThermostat(); tht_id++)
			_domain->setThT_heatFlux(tht_id, 0.0);
		    if(_doShearRate){
			    _domain->setPostShearEnergyConf(0.0);
			    _domain->setPreShearEnergyConf(0.0);
			    _domain->setPostShearEnergyDefault(0.0);
			    _domain->setPreShearEnergyDefault(0.0);
			    _domain->setShearVelDelta(0.0);
		    }
		}
		// test deletions and insertions for barostat
		if (_simstep >= _initGrandCanonical && _simstep <= _endGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
			  if(cpit->isGCMD_barostat()){
				if(!(_doRecordConfinement))
				  _domain->recordBarostat(_moleculeContainer);
				if ((!((_simstep + 2 * j + 4) % cpit->getInterval()) && _domain->getDifferentBarostatInterval() == false) || ((_simstep+1) == cpit->getInterval() && _domain->getDifferentBarostatInterval() == true))
				  _domain->collectBarostat(_domainDecomposition);
			  }
			  j++;
			}
		}
		
		if(_domainDecomposition->getRank()==0 && /*(_simstep > _initStatistics) &&*/ _doRecordBulkPressure && !(_simstep % _bulkPressureOutputTimesteps)){
		    BulkPressure << _simstep << "\t\t\t" << _domain->getGlobalBulkPressure() << "\t\t\t" << _domain->getGlobalBulkDensity() << "\t\t\t" << _domain->getGlobalBulkTemperature() << endl;
		}
		
		// Cancelling momentum for each component seperately 
		if(_cancelMomentum)
		  _domain->cancelMomentum(_domainDecomposition, _moleculeContainer);
		
		
		if(_forced_checkpoint_time >= 0 && (loopTimer.get_etime() + ioTimer.get_etime() + perStepIoTimer.get_etime()) >= _forced_checkpoint_time) {
			/* force checkpoint for specified time */
			string cpfile(_outputPrefix + ".timed.restart.xdr");
			global_log->info() << "Writing timed, forced checkpoint to file '" << cpfile << "'" << endl;
			_domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);
			_forced_checkpoint_time = -1; /* disable for further timesteps */
		}
		perStepIoTimer.stop();
		loopTimer.start();
	}
	loopTimer.stop();
	/***************************************************************************/
	/* END MAIN LOOP                                                           */
	/*****************************//**********************************************/

    ForceData.close();
    BulkPressure.close();
    Confinement.close();
    HeatFlux.close();
    Energy.close();
    ioTimer.start();
    if( _finalCheckpoint ) {
        /* write final checkpoint */
        string cpfile(_outputPrefix + ".restart.xdr");
        global_log->info() << "Writing final checkpoint to file '" << cpfile << "'" << endl;
        _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);
    }

	// dellocate stress tensor
	if(_doRecordStressProfile == true) {
	  size_t rows = 6, cols = _domain->getUniversalNProfileUnits_Stress(0)*_domain->getUniversalNProfileUnits_Stress(1)*_domain->getUniversalNProfileUnits_Stress(2);
	  _domain->dellocStressMatrix(_domain->_localStress, rows, cols);
	  _domain->dellocStressMatrix(_domain->_universalStress, rows, cols);
	}
	if(_doRecordConfinement == true) {
	  size_t rows = 6, cols = _domain->getUniversalNProfileUnitsStressConfinement(0)*_domain->getUniversalNProfileUnitsStressConfinement(1);
	  _domain->dellocStressMatrix(_domain->_localStressConfinement, rows, cols);
	  _domain->dellocStressMatrix(_domain->_globalStressConfinement, rows, cols);
	}
    

	// finish output
	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		(*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain);
		delete (*outputIter);
	}
	ioTimer.stop();

	global_log->info() << "Computation in main loop took: "
			<< loopTimer.get_etime() << " sec" << endl;
	global_log->info() << "Decomposition took. "
			<< decompositionTimer.get_etime() << " sec" << endl;
	global_log->info() << "IO in main loop  took:         "
			<< perStepIoTimer.get_etime() << " sec" << endl;
	global_log->info() << "Final IO took:                 "
			<< ioTimer.get_etime() << " sec" << endl;
			
	unsigned long numTimeSteps = _numberOfTimesteps - _initSimulation + 1; // +1 because of <= in loop
	double elapsed_time = loopTimer.get_etime() + decompositionTimer.get_etime();
	if(NULL != _ljFlopCounter) {
		double flop_rate = _ljFlopCounter->getTotalFlopCount() * numTimeSteps / elapsed_time / (1024*1024);
		global_log->info() << "LJ-FLOP-Count per Iteration: " << _ljFlopCounter->getTotalFlopCount() << " FLOPs" <<endl;
		global_log->info() << "FLOP-rate: " << flop_rate << " MFLOPS" << endl;
	}
}

void Simulation::output(unsigned long simstep) {

	int mpi_rank = _domainDecomposition->getRank();
	std::list<OutputBase*>::iterator outputIter;
	if (simstep >= _initStatistics){
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		OutputBase* output = (*outputIter);
		global_log->debug() << "Ouptut from " << output->getPluginName() << endl;
		output->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep, &this->_lmu, &this->_mcav);
	}
	}
	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileRecordingTimesteps)) {
		_domain->recordProfile(_moleculeContainer);
	}
	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileOutputTimesteps)) {
		_domain->collectProfile(_domainDecomposition);
		if (mpi_rank == 0) {
			ostringstream osstrm;
			osstrm << _profileOutputPrefix << ".";
			osstrm.fill('0');
			osstrm.width(9);
			osstrm << right << simstep;
			if(this->_domain->isCylindrical()){
				this->_domain->outputCylindricalProfile(osstrm.str().c_str());
			}
			else{
				_domain->outputProfile(osstrm.str().c_str());
			}
			osstrm.str("");
			osstrm.clear();
		}
		_domain->resetProfile();
	}

	if ((simstep >= _initStatistics) && _doRecordSlabProfile && !(simstep % _slabProfileRecordingTimesteps)) {
		_domain->recordSlabProfile(_moleculeContainer);
	}
	if ((simstep >= _initStatistics) && _doRecordSlabProfile && !(simstep % _slabProfileOutputTimesteps)) {
		_domain->collectSlabProfile(_domainDecomposition);
		if (mpi_rank == 0) {
			ostringstream osstrm;
			osstrm << _slabProfileOutputPrefix << ".";
			osstrm.fill('0');
			osstrm.width(9);
			osstrm << right << simstep;
			_domain->outputSlabProfile(osstrm.str().c_str());
			osstrm.str("");
			osstrm.clear();
		}
		_domain->resetSlabProfile();
	}
	if ((simstep >= _initStatistics) && _doRecordStressProfile && !(simstep % _stressProfileRecordingTimesteps)) {
		_domain->recordStressProfile(_moleculeContainer, simstep, _initStatistics);
		if (simstep == _initStatistics)
		      _domain->resetStressProfile();
	}
	if ((simstep > _initStatistics) && _doRecordStressProfile && !(simstep % _stressProfileOutputTimesteps)) {
		_domain->collectStressProfile(_domainDecomposition);
		if (mpi_rank == 0) {
			ostringstream osstrm;
			osstrm << _stressProfileOutputPrefix << ".";
			osstrm.fill('0');
			osstrm.width(9);
			osstrm << right << simstep;
			_domain->outputStressProfile(osstrm.str().c_str());
			osstrm.str("");
			osstrm.clear();
		}
		_domain->resetStressProfile();
	}

	if (/*(simstep >= _initStatistics) &&*/ _doRecordBulkPressure) {
		_domain->recordBulkPressure(_moleculeContainer);
		if (simstep == 1/**_initStatistics*/)
		      _domain->resetBulkPressure();
	}

	if (/*(simstep > _initStatistics) &&*/ _doRecordBulkPressure && !(simstep % _bulkPressureOutputTimesteps)) {
		_domain->collectBulkPressure(_domainDecomposition);
		if (mpi_rank == 0) {
		    double bulkPressure = (_domain->getPressureVirial() + _domain->getPressureKin())/(3 * _domain->getBulkVolume() * _domain->getAccumulatedDatasetsBulkPressure());
		    _domain->setGlobalBulkPressure(bulkPressure);
		    _domain->setGlobalBulkDensity(_domain->getBulkN()/(_domain->getBulkVolume() * _domain->getAccumulatedDatasetsBulkPressure()));
		    _domain->setGlobalBulkTemperature(_domain->getPressureKin()/(3*_domain->getBulkN()));
		    //cout << "BP: " << bulkPressure << " PV: " << _domain->getPressureVirial()/(3*_domain->getBulkVolume() * _domain->getAccumulatedDatasetsBulkPressure()) << " PK: " << _domain->getPressureKin()/(3*_domain->getBulkVolume() * _domain->getAccumulatedDatasetsBulkPressure()) << " rho: " << _domain->getBulkN()/(_domain->getBulkVolume() * _domain->getAccumulatedDatasetsBulkPressure()) << " T: " << _domain->getPressureKin()/(3*_domain->getBulkN()) << " Data: " << _domain->getAccumulatedDatasetsBulkPressure() << endl;
		}
		_domain->resetBulkPressure();
	}
	if ((simstep >= _initStatistics) && _doRecordConfinement  && !(simstep % _confinementRecordingTimesteps)) {
		_domain->recordConfinementProperties(_domainDecomposition, _moleculeContainer, simstep, _initStatistics);
	}
	else if ((simstep == _initStatistics) && _doRecordConfinement){
	        _domain->resetConfinementProperties();
	}
	if (_domain->thermostatWarning())
		global_log->warning() << "Thermostat!" << endl;
	/* TODO: thermostat */
	global_log->info() << "Simstep = " << simstep << "\tT = "
			<< _domain->getGlobalCurrentTemperature() << "\tU_pot = "
			<< _domain->getAverageGlobalUpot() << "\tp = "
			<< _domain->getGlobalPressure() << endl;
}

void Simulation::finalize() {
	if (_domainDecomposition != NULL) {
		delete _domainDecomposition;
		_domainDecomposition = NULL;
	}
	global_simulation = NULL;
}

void Simulation::updateParticleContainerAndDecomposition() {
	// The particles have moved, so the neighbourhood relations have
	// changed and have to be adjusted
	_moleculeContainer->update();
	//_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain);
	_domainDecomposition->balanceAndExchange(true, _moleculeContainer, _domain);
	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateMoleculeCaches();
}

/* FIXME: we shoud provide a more general way of doing this */
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

	_domainDecomposition = NULL;
	_domain = NULL;
	_particlePairsHandler = NULL;
	_cellProcessor = NULL;
	_moleculeContainer = NULL;
	_integrator = NULL;
	_inputReader = NULL;
    _finalCheckpoint = true;

#ifndef ENABLE_MPI
	global_log->info() << "Initializing the alibi domain decomposition ... " << endl;
	_domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
#else
	global_log->info() << "Initializing the standard domain decomposition ... " << endl;
	_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
#endif
	global_log->info() << "Initialization done" << endl;

	/*
	 * default parameters
	 */
	_cutoffRadius = 0.0;
	_LJCutoffRadius = 0.0;
	_tersoffCutoffRadius = 3.0;
	_numberOfTimesteps = 1;
	_outputPrefix = string("mardyn");
	_outputPrefix.append(gettimestring());
	_doRecordProfile = false;
	_doRecordSlabProfile = false;
	_doRecordBulkPressure = false;
	_doRecordConfinement = false;
	_doRecordStressProfile = false;
	_doEnergyOutputPerComponent = false;
	_doShearRate = false;
	_doShearForce = false;
	_doConfinementSplitStress = false;
	_cancelMomentum = false;
	_reduceDataSlab = false;
	_noXYZ = false;
	_boolDirectedVel = false;
	_directedVelocityTime = 1000;
	_doMovingAverageSigma2D = false;
	_doNeighbourAverage = false;
	_doSimpleAverage = false;
	_doSimpleAverageSigma2D = false;
	_doSimpleAverageSigma3D = false;
	_HardyStress = false;
	_HardyConfinement = false;
	_weightingStress = string("Linear");
	_weightingConfinement = string("Linear");
	_profileRecordingTimesteps = 7;
	_profileOutputTimesteps = 12500;
	_profileOutputPrefix = "out";
	_collectThermostatDirectedVelocity = 100;
	_zoscillation = false;
	_zoscillator = 512;
	_initCanonical = 5000;
	_initGrandCanonical = 10000000;
	_endGrandCanonical = 100000000;
	_initStatistics = 20000;
	h = 0.0;
	_pressureGradient = new PressureGradient(ownrank);
	global_log->info() << "Constructing domain ..." << endl;
	_domain = new Domain(ownrank, this->_pressureGradient);
	global_log->info() << "Domain construction done." << endl;
	_particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);
        
        this->_mcav = map<unsigned, CavityEnsemble>();
}

bool Simulation::isAcceleratingInstantaneously(){ return _domain->getPG()->isAcceleratingInstantaneously(_domain->getNumberOfComponents()); }

double Simulation::getMovedVel(int d, unsigned cid){ return _domain->getPG()->getGlobalVelSum(d,cid) / _domain->getPG()->getGlobalN(cid); }