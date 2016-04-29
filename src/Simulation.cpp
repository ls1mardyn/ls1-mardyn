#define SIMULATION_SRC
#include "Simulation.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>

#include "Common.h"
#include "Domain.h"
#include "particleContainer/LinkedCells.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/NonBlockingMPIHandlerBase.h"
#include "molecules/Molecule.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/adapter/FlopCounter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "molecules/Wall.h"
#include "molecules/Mirror.h"

#include "io/io.h"
#include "io/GeneratorFactory.h"
#include "io/RDF.h"
#include "io/TcTS.h"
#include "io/Mkesfera.h"

#include "ensemble/GrandCanonical.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/PressureGradient.h"

#include "thermostats/VelocityScalingThermostat.h"
#include "thermostats/TemperatureControl.h"

#include "utils/FileUtils.h"
#include "utils/OptionParser.h"
#include "utils/Timer.h"
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
	: _simulationTime(0),
	_initStatistics(0),
	_rdf(NULL),
	_moleculeContainer(NULL),
	_particlePairsHandler(NULL),
	_cellProcessor(NULL),
	_flopCounter(NULL),
	_domainDecomposition(NULL),
	_integrator(NULL),
	_domain(NULL),
	_inputReader(NULL),
	_longRangeCorrection(NULL),
	_temperatureControl(NULL),
	_FMM(NULL),
	_forced_checkpoint_time(0)
{
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
	delete _flopCounter;
	delete _FMM;
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
	if(xmlconfig.changecurrentnode("integrator")) {
		string integratorType;
		xmlconfig.getNodeValue("@type", integratorType);
		global_log->info() << "Integrator type: " << integratorType << endl;
		if(integratorType == "Leapfrog") {
			_integrator = new Leapfrog();
		} else {
			global_log-> error() << "Unknown integrator " << integratorType << endl;
			this->exit(1);
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
			_ensemble = new CanonicalEnsemble();
		} else if (ensembletype == "muVT") {
			global_log->error() << "muVT ensemble not completely implemented via XML input." << endl;
			this->exit(1);
			// _ensemble = new GrandCanonicalEnsemble();
		} else {
			global_log->error() << "Unknown ensemble type: " << ensembletype << endl;
			this->exit(1);
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
		this->exit(1);
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
				this->exit(1);
			}
			global_log->info() << "dimensionless cutoff radius:\t" << _cutoffRadius << endl;
			xmlconfig.changecurrentnode("..");
		} else {
			global_log->error() << "Cutoff section missing." << endl;
			this->exit(1);
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
			this->exit(1);
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
			if(parallelisationtype == "DummyDecomposition") {
				global_log->error() << "DummyDecomposition not available in parallel mode." << endl;
				//_domainDecomposition = new DomainDecompDummy();
			}
			else if(parallelisationtype == "DomainDecomposition") {
				DomainDecomposition * temp = 0;
							temp = dynamic_cast<DomainDecomposition *>(_domainDecomposition);
							if (temp != 0) {
								temp->initCommunicationPartners(getcutoffRadius(), _domain);
							}
							_domainDecomposition = new DomainDecomposition();
			}
			else if(parallelisationtype == "KDDecomposition") {
				_domainDecomposition = new KDDecomposition(getcutoffRadius(), _domain);
			}
			else {
				global_log->error() << "Unknown parallelisation type: " << parallelisationtype << endl;
				this->exit(1);
			}
		#else /* serial */
			if(parallelisationtype != "DummyDecomposition") {
				global_log->warning()
						<< "Executable was compiled without support for parallel execution: "
						<< parallelisationtype
						<< " not available. Using serial mode." << endl;
				//this->exit(1);
			}
			//_domainDecomposition = new DomainDecompBase();  // already set in initialize()
		#endif
			_domainDecomposition->readXML(xmlconfig);
			xmlconfig.changecurrentnode("..");
		}
		else {
		#ifdef ENABLE_MPI
			global_log->error() << "Parallelisation section missing." << endl;
			this->exit(1);
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
			}
			else if(datastructuretype == "AdaptiveSubCells") {
				global_log->warning() << "AdaptiveSubCells no longer supported." << std::endl;
				global_simulation->exit(-1);
			}
			else {
				global_log->error() << "Unknown data structure type: " << datastructuretype << endl;
				this->exit(1);
			}
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
		
		
	/*	if(xmlconfig.changecurrentnode("planarLRC")) {
			XMLfile::Query query = xmlconfig.query("slabs");
			unsigned slabs = query.card();
			_longRangeCorrection = new Planar(_cutoffRadius, _LJCutoffRadius, _domain, _domainDecomposition, _moleculeContainer, slabs, global_simulation);
			_domainDecomposition->readXML(xmlconfig);
			xmlconfig.changecurrentnode("..");
		}
		else {
			_longRangeCorrection = new Homogeneous(_cutoffRadius, _LJCutoffRadius,_domain,global_simulation);	
		}*/

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
		else if(pluginname == "FLOPCounter") {
			/** @todo  Make the Flop counter a real output plugin */
			_flopCounter = new FlopCounter(_cutoffRadius, _LJCutoffRadius);
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
		else if(pluginname == "MmspdBinWriter") {
			outputPlugin = new MmspdBinWriter();
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
		else if(pluginname == "VectorizationTuner") {
			outputPlugin = new VectorizationTuner(_cutoffRadius, _LJCutoffRadius, &_cellProcessor);
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
	string extension(getFileExtension(filename.c_str()));
	global_log->debug() << "Found config filename extension: " << extension << endl;
	if (extension == "xml") {
		initConfigXML(filename);
	}
	else if (extension == "cfg") {
		initConfigOldstyle(filename);
	}
	else {
		global_log->error() << "Unknown config file extension '" << extension << "'." << endl;
		this->exit(1);;
	}
}

void Simulation::initConfigXML(const string& inputfilename) {
	global_log->info() << "Initializing XML config file: " << inputfilename << endl;
	XMLfileUnits inp(inputfilename);

	global_log->debug() << "Input XML:" << endl << string(inp) << endl;

	if(inp.changecurrentnode("/mardyn") < 0) {
		global_log->error() << "Cound not find root node /mardyn." << endl;
		global_log->error() << "Not a valid MarDyn XML input file." << endl;
		this->exit(1);
	}

	string version("unknown");
	inp.getNodeValue("@version", version);
	global_log->info() << "MarDyn XML config file version: " << version << endl;

	if (inp.changecurrentnode("simulation")) {
		/** @todo this is all for old input files. Remove! */
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

#ifdef ENABLE_MPI
	// if we are using the DomainDecomposition, please complete its initialization:
	{
		DomainDecomposition * temp = 0;
		temp = dynamic_cast<DomainDecomposition *>(_domainDecomposition);
		if (temp != 0) {
			temp->initCommunicationPartners(_cutoffRadius, _domain);
		}
	}
#endif

	// read particle data (or generate particles, if a generator is chosen)
	unsigned long maxid = _inputReader->readPhaseSpace(_moleculeContainer,
			&_lmu, _domain, _domainDecomposition);


	_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
//	_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

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

void Simulation::prepare_start() {
	global_log->info() << "Initializing simulation" << endl;

	global_log->info() << "Initialising cell processor" << endl;
#if ENABLE_VECTORIZED_CODE
	global_log->debug() << "Checking if vectorized cell processor can be used" << endl;
	bool lj_present = false;
	bool charge_present = false;
	bool dipole_present = false;
	bool quadrupole_present = false;

	const vector<Component> components = *(global_simulation->getEnsemble()->components());
	for (size_t i = 0; i < components.size(); i++) {
		lj_present |= (components[i].numLJcenters() != 0);
		charge_present |= (components[i].numCharges() != 0);
		dipole_present |= (components[i].numDipoles() != 0);
		quadrupole_present |= (components[i].numQuadrupoles() != 0);
	}
	global_log->debug() << "xx lj present: " << lj_present << endl;
	global_log->debug() << "xx charge present: " << charge_present << endl;
	global_log->debug() << "xx dipole present: " << dipole_present << endl;
	global_log->debug() << "xx quadrupole present: " << quadrupole_present << endl;

	if(this->_lmu.size() > 0) {
		global_log->warning() << "Using legacy cell processor. (The vectorized code does not support grand canonical simulations.)" << endl;
		_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _particlePairsHandler);
	}
	else if(this->_doRecordVirialProfile) {
		global_log->warning() << "Using legacy cell processor. (The vectorized code does not support the virial tensor and the localized virial profile.)" << endl;
		_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _particlePairsHandler);
	}
	else {
		global_log->info() << "Using vectorized cell processor." << endl;
		_cellProcessor = new VectorizedCellProcessor( *_domain, _cutoffRadius, _LJCutoffRadius);
	}
#else
	global_log->info() << "Using legacy cell processor." << endl;
	_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _particlePairsHandler);
#endif

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
				dynamic_cast<LinkedCells*>(_moleculeContainer)->cellLength());

		delete _cellProcessor;
		_cellProcessor = new bhfmm::VectorizedLJP2PCellProcessor(*_domain, _LJCutoffRadius, _cutoffRadius);
	}

	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();
	global_log->info() << "Updating domain decomposition" << endl;
	updateParticleContainerAndDecomposition();
	global_log->info() << "Performing inital force calculation" << endl;
	_moleculeContainer->traverseCells(*_cellProcessor);
	if (_FMM != NULL) {
		global_log->info() << "Performing inital FMM force calculation" << endl;
		_FMM->computeElectrostatics(_moleculeContainer);
	}

	/* If enabled count FLOP rate of LS1. */
	if( NULL != _flopCounter ) {
		_moleculeContainer->traverseCells(*_flopCounter);
	}

    // clear halo
    global_log->info() << "Clearing halos" << endl;
    _moleculeContainer->deleteOuterParticles();

    if (_longRangeCorrection == NULL){
        _longRangeCorrection = new Homogeneous(_cutoffRadius, _LJCutoffRadius,_domain,global_simulation);
    }


    _longRangeCorrection->calculateLongRange();

	// here we have to call calcFM() manually, otherwise force and moment are not
	// updated inside the molecule (actually this is done in upd_postF)
	// or should we better call the integrator->eventForcesCalculated?
	for (Molecule* tM = _moleculeContainer->begin(); tM
			!= _moleculeContainer->end(); tM = _moleculeContainer->next()) {
		tM->calcFM();
	}

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

	if (_lmu.size() > 0) {
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
	}

	// initialize output
	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter
			!= _outputPlugins.end(); outputIter++) {
		(*outputIter)->initOutput(_moleculeContainer, _domainDecomposition,
				_domain);
	}

	global_log->info() << "System initialised\n" << endl;
	global_log->info() << "System contains "
			<< _domain->getglobalNumMolecules() << " molecules." << endl;


}

void Simulation::simulate() {

	// by Stefan Becker
	Molecule* tM;

	global_log->info() << "Started simulation" << endl;

	// (universal) constant acceleration (number of) timesteps
	unsigned uCAT = _pressureGradient->getUCAT();
// 	_initSimulation = (unsigned long) (_domain->getCurrentTime()
// 			/ _integrator->getTimestepLength());
    _initSimulation = (unsigned long) (this->_simulationTime / _integrator->getTimestepLength());
	// _initSimulation = 1;
	/* demonstration for the usage of the new ensemble class */
	/*CanonicalEnsemble ensemble(_moleculeContainer, global_simulation->getEnsemble()->components());
	ensemble.updateGlobalVariable(NUM_PARTICLES);
	global_log->debug() << "Number of particles in the Ensemble: "
			<< ensemble.N() << endl;
	ensemble.updateGlobalVariable(ENERGY);
	global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E()
		<< endl;
	ensemble.updateGlobalVariable(TEMPERATURE);
	global_log->debug() << "Temperature of the Ensemble: " << ensemble.T()
		<< endl;*/

	/***************************************************************************/
	/* BEGIN MAIN LOOP                                                         */
	/***************************************************************************/

	// all timers except the ioTimer messure inside the main loop
	Timer loopTimer; /* timer for the entire simulation loop (synced) */
	Timer decompositionTimer; /* timer for decomposition */
	Timer computationTimer; /* timer for computation */
	Timer perStepIoTimer; /* timer for io in simulation loop */
	Timer ioTimer; /* timer for final io */

	loopTimer.set_sync(true);
#if WITH_PAPI
	const char *papi_event_list[] = {
		"PAPI_TOT_CYC",
		"PAPI_TOT_INS"
	//	"PAPI_VEC_DP"
	// 	"PAPI_L2_DCM"
	// 	"PAPI_L2_ICM"
	// 	"PAPI_L1_ICM"
	//	"PAPI_DP_OPS"
	// 	"PAPI_VEC_INS"
	};
	int num_papi_events = sizeof(papi_event_list) / sizeof(papi_event_list[0]);
	loopTimer.add_papi_counters(num_papi_events, (char**) papi_event_list);
#endif
	loopTimer.start();

	for (_simstep = _initSimulation; _simstep <= _numberOfTimesteps; _simstep++) {
		global_log->debug() << "timestep: " << getSimulationStep() << endl;
		global_log->debug() << "simulation time: " << getSimulationTime() << endl;

		computationTimer.start();

		/** @todo What is this good for? Where come the numbers from? Needs documentation */
		if (_simstep >= _initGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
				if (!((_simstep + 2 * j + 3) % cpit->getInterval())) {
					cpit->prepareTimestep(_moleculeContainer,
							_domainDecomposition);
				}
				j++;
			}
		}

		_integrator->eventNewTimestep(_moleculeContainer, _domain);

		// activate RDF sampling
		if ((_simstep >= _initStatistics) && _rdf != NULL) {
			global_log->info() << "Activating the RDF sampling" << endl;
			this->_rdf->tickRDF();
			this->_particlePairsHandler->setRDF(_rdf);
			this->_rdf->accumulateNumberOfMolecules(*(global_simulation->getEnsemble()->components()));
		}

		/*! by Stefan Becker <stefan.becker@mv.uni-kl.de> 
		 *realignment tools borrowed from Martin Horsch, for the determination of the centre of mass 
		 *the halo MUST NOT be present*/
#ifndef NDEBUG 
#ifndef ENABLE_MPI
			unsigned particleNoTest = 0;
                        for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next())
                        particleNoTest++;
                        global_log->info()<<"particles before determine shift-methods, halo not present:" << particleNoTest<< "\n";
#endif
#endif
		 if(_doAlignCentre && !(_simstep % _alignmentInterval))
		{
			if(_componentSpecificAlignment){
				//! !!! the sequence of calling the two methods MUST be: FIRST determineXZShift() THEN determineYShift() !!!
				_domain->determineXZShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
				_domain->determineYShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
			// edited by Michaela Heier --> realign can be used when LJ93-Potential will be used. Only the shift in the xz-plane will be used. 
			else if(_doAlignCentre && _applyWallFun){
				global_log->info() << "realign in the xz-plane without a shift in y-direction\n";
				_domain->determineXZShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
				_domain->noYShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
			else{
				_domain->determineShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
#ifndef NDEBUG 
#ifndef ENABLE_MPI			
			particleNoTest = 0;
			for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) 
			particleNoTest++;
			global_log->info()<<"particles after determine shift-methods, halo not present:" << particleNoTest<< "\n";
#endif
#endif
		}
		computationTimer.stop();

		bool overlapCommComp = true;

		if (overlapCommComp) {
			performOverlappingDecompositionAndCellTraversalStep(decompositionTimer, computationTimer);
		} else {
			decompositionTimer.start();
			// ensure that all Particles are in the right cells and exchange Particles
			global_log->debug() << "Updating container and decomposition"
					<< endl;
			updateParticleContainerAndDecomposition();
			decompositionTimer.stop();

			computationTimer.start();
			// Force calculation and other pair interaction related computations
			global_log->debug() << "Traversing pairs" << endl;
			_moleculeContainer->traverseCells(*_cellProcessor);
			computationTimer.stop();
		}
		computationTimer.start();
		if (_FMM != NULL) {
			global_log->debug() << "Performing FMM calculation" << endl;
			_FMM->computeElectrostatics(_moleculeContainer);
		}

		




		if(_wall && _applyWallFun){
		  _wall->calcTSLJ_9_3(_moleculeContainer, _domain);
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
					global_log->debug() << "Grand canonical ensemble(" << j
							<< "): test deletions and insertions" << endl;
                                        this->_domain->setLambda(cpit->getLambda());
                                        this->_domain->setDensityCoefficient(cpit->getDensityCoefficient());
                                        double localUpotBackup = _domain->getLocalUpot();
                                        double localVirialBackup = _domain->getLocalVirial();
					_moleculeContainer->grandcanonicalStep(&(*cpit),
							_domain->getGlobalCurrentTemperature(), this->_domain, *_cellProcessor);
                                        _domain->setLocalUpot(localUpotBackup);
                                        _domain->setLocalVirial(localVirialBackup);
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
				}

				j++;
			}
		}

		/*! by Stefan Becker <stefan.becker@mv.uni-kl.de> 
		  * realignment tools borrowed from Martin Horsch
		  * For the actual shift the halo MUST be present!
		  */
		
		if(_doAlignCentre && !(_simstep % _alignmentInterval))
		{
			_domain->realign(_moleculeContainer);
#ifndef NDEBUG 
#ifndef ENABLE_MPI
			particleNoTest = 0;
			for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) 
			particleNoTest++;
			global_log->info()<<"particles after realign(), halo present:" << particleNoTest<< "\n";
#endif
#endif
		}
		
		// clear halo

		global_log->debug() << "Deleting outer particles / clearing halo." << endl;
		_moleculeContainer->deleteOuterParticles();

		/** @todo For grand canonical ensemble? Sould go into appropriate ensemble class. Needs documentation. */
		if (_simstep >= _initGrandCanonical) {
			_domain->evaluateRho(_moleculeContainer->getNumberOfParticles(),
					_domainDecomposition);
		}

		if (!(_simstep % _collectThermostatDirectedVelocity))
			_domain->calculateThermostatDirectedVelocity(_moleculeContainer);
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
		}
		_longRangeCorrection->calculateLongRange();
		
		/*
		 * radial distribution function
		 */
		if (_simstep >= _initStatistics) {
			if (this->_lmu.size() == 0) {
				this->_domain->record_cv();
			}
		}

		// Inform the integrator about the calculated forces
		global_log->debug() << "Inform the integrator" << endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);

		// calculate the global macroscopic values from the local values
		global_log->debug() << "Calculate macroscopic values" << endl;
		_domain->calculateGlobalValues(_domainDecomposition,
				_moleculeContainer, (!(_simstep % _collectThermostatDirectedVelocity)), Tfactor(
								_simstep));//*/  //TODO: uncomment !!!!!
		
		// scale velocity and angular momentum
		if ( !_domain->NVE() && _temperatureControl == NULL){
		if (_thermostatType ==VELSCALE_THERMOSTAT){
			global_log->debug() << "Velocity scaling" << endl;
			if (_domain->severalThermostats()) {
				_velocityScalingThermostat.enableComponentwise();
				for(unsigned int cid = 0; cid < global_simulation->getEnsemble()->components()->size(); cid++) {
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
		  else if(_thermostatType == ANDERSEN_THERMOSTAT){ //! the Andersen Thermostat
//		    global_log->info() << "Andersen Thermostat" << endl;
		    double nuDt = _nuAndersen * _integrator->getTimestepLength();
//		    global_log->info() << "Timestep length = " << _integrator->getTimestepLength() << " nuDt = " << nuDt << "\n";
		    unsigned numPartThermo = 0; // for testing reasons
		    double tTarget;
		    double stdDevTrans, stdDevRot;
		    if(_domain->severalThermostats()) {
		      for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
			if (_rand.rnd() < nuDt){
			  numPartThermo++;
			  int thermostat = _domain->getThermostat(tM->componentid());
			  tTarget = _domain->getTargetTemperature(thermostat);
			  stdDevTrans = sqrt(tTarget/tM->gMass());
			  for(unsigned short d = 0; d < 3; d++){
			    stdDevRot = sqrt(tTarget*tM->getI(d));
			    tM->setv(d,_rand.gaussDeviate(stdDevTrans));
			    tM->setD(d,_rand.gaussDeviate(stdDevRot));
			  }
			}
		      }
		    }
		    else{
		      tTarget = _domain->getTargetTemperature(0);
		      for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
			if (_rand.rnd() < nuDt){
			  numPartThermo++;
			  // action of the anderson thermostat: mimic a collision by assigning a maxwell distributed velocity
			  stdDevTrans = sqrt(tTarget/tM->gMass());
			  for(unsigned short d = 0; d < 3; d++){
			    stdDevRot = sqrt(tTarget*tM->getI(d));
			    tM->setv(d,_rand.gaussDeviate(stdDevTrans));
			    tM->setD(d,_rand.gaussDeviate(stdDevRot));
			  }
			}
		      }
		    }
//		    global_log->info() << "Andersen Thermostat: n = " << numPartThermo ++ << " particles thermostated\n";
		  }


		// if(_mirror && _applyMirror){
		//  _mirror->VelocityChange(_moleculeContainer, _domain);
		//}


		}
		// mheinen 2015-07-27 --> TEMPERATURE_CONTROL
        else if ( _temperatureControl != NULL)
        {
            _temperatureControl->DoLoopsOverMolecules(_domainDecomposition, _moleculeContainer, _simstep);
        }
        // <-- TEMPERATURE_CONTROL
		

		advanceSimulationTime(_integrator->getTimestepLength());

		/* BEGIN PHYSICAL SECTION:
		 * the system is in a consistent state so we can extract global variables
		 */
		/*ensemble.updateGlobalVariable(NUM_PARTICLES);
		global_log->debug() << "Number of particles in the Ensemble: "
				<< ensemble.N() << endl;
		ensemble.updateGlobalVariable(ENERGY);
		global_log->debug() << "Kinetic energy in the Ensemble: "
				<< ensemble.E() << endl;
		ensemble.updateGlobalVariable(TEMPERATURE);
		global_log->debug() << "Temperature of the Ensemble: " << ensemble.T()
			<< endl;*/
		/* END PHYSICAL SECTION */

		computationTimer.stop();
		perStepIoTimer.start();

		output(_simstep);
		if(_forced_checkpoint_time >= 0 && (loopTimer.get_etime() + ioTimer.get_etime() + perStepIoTimer.get_etime()) >= _forced_checkpoint_time) {
			/* force checkpoint for specified time */
			string cpfile(_outputPrefix + ".timed.restart.xdr");
			global_log->info() << "Writing timed, forced checkpoint to file '" << cpfile << "'" << endl;
			_domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition, _simulationTime);
			_forced_checkpoint_time = -1; /* disable for further timesteps */
		}
		perStepIoTimer.stop();
	}
	loopTimer.stop();
	/***************************************************************************/
	/* END MAIN LOOP                                                           */
	/*****************************//**********************************************/

    ioTimer.start();
    if( _finalCheckpoint ) {
        /* write final checkpoint */
        string cpfile(_outputPrefix + ".restart.xdr");
        global_log->info() << "Writing final checkpoint to file '" << cpfile << "'" << endl;
        _domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition, _simulationTime);
    }
	// finish output
	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		(*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain);
		delete (*outputIter);
	}
	ioTimer.stop();

	global_log->info() << "Computation in main loop took: " << loopTimer.get_etime() << " sec" << endl;
	global_log->info() << "Decomposition took:            " << decompositionTimer.get_etime() << " sec" << endl;
	global_log->info() << "IO in main loop  took:         " << perStepIoTimer.get_etime() << " sec" << endl;
	global_log->info() << "Final IO took:                 " << ioTimer.get_etime() << " sec" << endl;

#if WITH_PAPI
	global_log->info() << "PAPI counter values for loop timer:"  << endl;
	for(int i = 0; i < loopTimer.get_papi_num_counters(); i++) {
		global_log->info() << "  " << papi_event_list[i] << ": " << loopTimer.get_global_papi_counter(i) << endl;
	}
#endif /* WITH_PAPI */

	unsigned long numTimeSteps = _numberOfTimesteps - _initSimulation + 1; // +1 because of <= in loop
	double elapsed_time = loopTimer.get_etime() + decompositionTimer.get_etime();
	if(NULL != _flopCounter) {
		double flop_rate = _flopCounter->getTotalFlopCount() * numTimeSteps / elapsed_time / (1024*1024);
		global_log->info() << "FLOP-Count per Iteration: " << _flopCounter->getTotalFlopCount() << " FLOPs" <<endl;
		global_log->info() << "FLOP-rate: " << flop_rate << " MFLOPS" << endl;
	}
}

void Simulation::output(unsigned long simstep) {

	int mpi_rank = _domainDecomposition->getRank();

	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		OutputBase* output = (*outputIter);
		global_log->debug() << "Output from " << output->getPluginName() << endl;
		output->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep, &(_lmu));
	}

	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileRecordingTimesteps)) {
		_domain->recordProfile(_moleculeContainer, _doRecordVirialProfile);
	}
	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileOutputTimesteps)) {
		_domain->collectProfile(_domainDecomposition, _doRecordVirialProfile);
		if (mpi_rank == 0) {
			ostringstream osstrm;
			osstrm << _profileOutputPrefix << ".";
			osstrm.fill('0');
			osstrm.width(9);
			osstrm << right << simstep;
			//edited by Michaela Heier 
			if(this->_domain->isCylindrical()){
				this->_domain->outputCylProfile(osstrm.str().c_str());
				_domain->outputProfile(osstrm.str().c_str(),_doRecordVirialProfile);
			}
			else{
			_domain->outputProfile(osstrm.str().c_str(), _doRecordVirialProfile);
			}
			osstrm.str("");
			osstrm.clear();
		}
		_domain->resetProfile(_doRecordVirialProfile);
	}

	
	if (_domain->thermostatWarning())
		global_log->warning() << "Thermostat!" << endl;
	/* TODO: thermostat */
	global_log->info() << "Simstep = " << simstep << "\tT = "
			<< _domain->getGlobalCurrentTemperature() << "\tU_pot = "
			<< _domain->getGlobalUpot() << "\tp = "
			<< _domain->getGlobalPressure() << endl;
}

void Simulation::finalize() {
	if (_FMM != NULL) {
		_FMM->printTimers();
		bhfmm::VectorizedLJP2PCellProcessor * temp = dynamic_cast<bhfmm::VectorizedLJP2PCellProcessor*>(_cellProcessor);
		temp->printTimers();
	}

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
	bool forceRebalancing = false;
	_domainDecomposition->balanceAndExchange(forceRebalancing, _moleculeContainer, _domain);
	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateMoleculeCaches();
}

void Simulation::performOverlappingDecompositionAndCellTraversalStep(
		Timer& decompositionTimer, Timer& computationTimer) {

	//_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain);
	bool forceRebalancing = false;

	//TODO: exchange the constructor for a real non-blocking version

#ifdef ENABLE_MPI
#ifdef ENABLE_OVERLAPPING
	NonBlockingMPIHandlerBase* nonBlockingMPIHandler = new NonBlockingMPIMultiStepHandler(&decompositionTimer, &computationTimer,
			static_cast<DomainDecompMPIBase*>(_domainDecomposition), _moleculeContainer, _domain, _cellProcessor);
#else
	NonBlockingMPIHandlerBase* nonBlockingMPIHandler = new NonBlockingMPIHandlerBase(&decompositionTimer, &computationTimer,
				_domainDecomposition, _moleculeContainer, _domain, _cellProcessor);
#endif

	nonBlockingMPIHandler->performOverlappingTasks(forceRebalancing);
#endif
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

	_finalCheckpoint = true;

        // TODO:
#ifndef ENABLE_MPI
	global_log->info() << "Initializing the alibi domain decomposition ... " << endl;
	_domainDecomposition = new DomainDecompBase();
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
	_numberOfTimesteps = 1;
	_outputPrefix = string("mardyn");
	_outputPrefix.append(gettimestring());

	/** @todo the following features should be documented */
	_doRecordProfile = false;
	_doRecordVirialProfile = false;
	_profileRecordingTimesteps = 7;
	_profileOutputTimesteps = 12500;
	_profileOutputPrefix = "out";
	_collectThermostatDirectedVelocity = 100;
	_initCanonical = 5000;
	_initGrandCanonical = 10000000;
	_initStatistics = 20000;
	h = 0.0;

	_thermostatType = VELSCALE_THERMOSTAT;
	_nuAndersen = 0.0;
	_rand.init(8624);

	_doAlignCentre = false;
	_componentSpecificAlignment = false;
	_alignmentInterval = 25;
	_momentumInterval = 1000;
	_wall = NULL;
	_applyWallFun = false;
	_mirror = NULL;
	_applyMirror = false;

	_pressureGradient = new PressureGradient(ownrank);
	global_log->info() << "Constructing domain ..." << endl;
	_domain = new Domain(ownrank, this->_pressureGradient);
	global_log->info() << "Domain construction done." << endl;
	_particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);
	_longRangeCorrection = NULL;
}
