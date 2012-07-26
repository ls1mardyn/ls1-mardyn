/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// Simulation.cpp
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

#define SIMULATION_SRC
#include "Simulation.h"

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/AdaptiveSubCells.h"
#include "parallel/DomainDecompBase.h"
//#include "ParticleInsertion.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#include "parallel/KDDecomposition2.h"
#else
#include "parallel/DomainDecompDummy.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"

#include "io/io.h"
#include "io/xmlreader.h"
#include "io/GeneratorFactory.h"

#include "ensemble/GrandCanonical.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/PressureGradient.h"

#include "RDF.h"

#ifdef STEEREO
#include "utils/SteereoIntegration.h"
#include <steereo/steereoSimSteering.h>
#include <steereo/steereoCouplingSim.h>
#endif

#include "utils/OptionParser.h"
#include "utils/Timer.h"
#include "utils/Logger.h"

#include "io/TcTS.h"

using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;

Simulation* global_simulation;

Simulation::Simulation() :
	_rdf(NULL), _domainDecomposition(NULL) {

	initialize();
}

Simulation::~Simulation() {
	if (_rdf)
		delete _rdf;
	if (_domainDecomposition)
		delete _domainDecomposition;
	if (_pressureGradient)
		delete _pressureGradient;
	if (_domain)
		delete _domain;
	if (_particlePairsHandler)
		delete _particlePairsHandler;
	if (_moleculeContainer)
		delete _moleculeContainer;
	if (_integrator)
		delete _integrator;
	if (_inputReader)
		delete _inputReader;
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
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif
	global_log->info() << "init XML config file: " << inputfilename << endl;
	XMLfileUnits inp(inputfilename);
	XMLReader xmlreader(inputfilename);

	global_log->debug() << "Input XML:" << endl << string(inp) << endl;
	inp.changecurrentnode("/mardyn");

	string version = xmlreader.getVersion();
	global_log->info() << "MarDyn XML config file version: " << version << endl;
	global_log->info() << "Reading simulation parameters ..." << endl;

	global_log->info() << "Reading integrator information..." << endl;
	xmlreader.getIntegrator(_integrator);
	double timestepLength = 0.0;
	timestepLength = _integrator->getTimestepLength();
	global_log->info() << " timestep:             " << timestepLength << endl;
	if (timestepLength == 0) {
		/* TODO: perform checks before simulation start as there are command line parameters, too */
		global_log->error() << "Undefined timestep length." << endl;
	}

	global_log->info() << "Reading parallelization type..." << endl;
	_domainDecompositionType = xmlreader.getDomainDecompositionType();
	switch (_domainDecompositionType) {
	case DUMMY_DECOMPOSITION:
#ifdef ENABLE_MPI
		global_log->warning() << "Dummy domain decomposition does not work in a parallel environment." << endl;
		exit(1);
#else
		_domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
#endif
		break;
#ifdef ENABLE_MPI
		case DOMAIN_DECOMPOSITION:
		_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
		break;
		case KD_DECOMPOSITION:
		/* TODO: remove parameters from KDDecomposition constructor */
		_domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, 1.0, 100);
		break;
#endif
	case UNKNOWN_DECOMPOSITION:
		global_log->error() << "Unknown domain decomposition type." << endl;
		exit(1);
		break;
	}
#ifndef ENABLE_MPI
	if (_domainDecompositionType != DUMMY_DECOMPOSITION) {
		global_log->warning()
				<< "Input demands parallelization, but the current compilation doesn't support parallel execution."
				<< endl;
	}
#endif	

	double simBoxLength[3];
	if (xmlreader.getSimBoxSize(simBoxLength)) {
		global_log->info() << " simulation box size:  (" << simBoxLength[0]
				<< ", " << simBoxLength[1] << ", " << simBoxLength[2] << ")"
				<< endl;

		for (int d = 0; d < 3; d++) {
			_domain->setGlobalLength(d, simBoxLength[d]);
		}
	}

	double temperature = xmlreader.getTemperature();
	_domain->setGlobalTemperature(temperature);
	global_log->info() << "Temperature of the enseble: " << temperature << endl;

	if (inp.changecurrentnode("simulation")) {
		string siminpfile, siminptype;
		if (inp.getNodeValue("input", siminpfile)) {
			global_log->info() << "reading input file:\t" << siminpfile << endl;
			// input type="oldstyle" to include oldstyle input files for backward compatibility - only temporary!!!
			if (inp.getNodeValue("input@type", siminptype) && siminptype
					== "oldstyle") {
				global_log->info() << "         file type:\t" << siminptype
						<< endl;
				initConfigOldstyle(siminpfile);
			}
		}

		if (inp.getNodeValueReduced("cutoffs/radiusLJ", _LJCutoffRadius)) {
			_cutoffRadius = _LJCutoffRadius;
			global_log->info() << "dimensionless LJ cutoff radius:\t"
					<< _LJCutoffRadius << endl;
			global_log->info() << "dimensionless cutoff radius:\t"
					<< _cutoffRadius << endl;
		} else {
			global_log->error() << "Cutoff section missing." << endl;
			exit(1);
		}

		string pspfile;
		if (inp.getNodeValue("ensemble/phasespacepoint/file", pspfile)) {
			pspfile.insert(0, inp.getDir());
			global_log->info() << "phasespacepoint description file:\t"
					<< pspfile << endl;

			string pspfiletype("OldStyle");
			inp.getNodeValue("ensemble/phasespacepoint/file@type", pspfiletype);
			global_log->info() << "       phasespacepoint file type:\t"
					<< pspfiletype << endl;
			if (pspfiletype == "OldStyle") {
				_inputReader = (InputBase*) new InputOldstyle();
				_inputReader->setPhaseSpaceFile(pspfile);
				_inputReader->setPhaseSpaceHeaderFile(pspfile);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else if (pspfiletype == "PartGen") {
				string mode;
				inp.getNodeValue("ensemble/phasespacepoint/file@mode", mode);
				global_log->error() << "NOT IMPLEMENTED YET";
				/*
				 _inputReader = (InputBase*) new PartGen();
				 _inputReader->setPhaseSpaceHeaderFile(pspfile);
				 // PartGen has to modes, a "Homogeneous" mode, where particles
				 // with a homogeneous distribution are created, and a "Cluster" mode,
				 // where droplets are created. Currently, only the cluster mode is supported,
				 // which needs another config line starting with "clusterFile ..." directly
				 // after the config line starting with "phaseSpaceFile"
				 double gasDensity;
				 double fluidDensity;
				 double volPercOfFluid;
				 string clusterFileName;
				 inputfilestream >> token >> gasDensity >> fluidDensity >> volPercOfFluid >> clusterFileName;
				 ((PartGen*) _inputReader)->setClusterFile(gasDensity, fluidDensity, volPercOfFluid, clusterFileName);
				 ((PartGen*) _inputReader)->readPhaseSpaceHeader(_domain, timestepLength);
				 */
			} else if (pspfiletype == "1CLJGen") {
				string mode;
				inp.getNodeValue("ensemble/phasespacepoint/file@mode", mode);
				global_log->error() << "NOT IMPLEMENTED YET";
				/*
				 int N;
				 double T;
				 string line;
				 getline(inputfilestream, line);
				 stringstream lineStream(line);
				 lineStream >> mode >> N >> T;
				 global_log->debug() << "read: mode " << mode << " N " << N << " T " << T << endl;

				 OneCLJGenerator* generator = (OneCLJGenerator*) new OneCLJGenerator(mode, N, T);
				 if (mode == "Homogeneous") {
				 double rho;
				 lineStream >> rho;
				 generator->setHomogeneuosParameter(rho);
				 } else if (mode == "Cluster") {
				 double rho_gas, rho_fluid, fluidVolumePercent, maxSphereVolume, numSphereSizes;
				 lineStream >> rho_gas >> rho_fluid >> fluidVolumePercent >> maxSphereVolume >> numSphereSizes;
				 generator->setClusterParameters(rho_gas, rho_fluid, fluidVolumePercent, maxSphereVolume, numSphereSizes);
				 } else {
				 global_log->error() << "Error in inputfile: OneCLJGenerator option \""<< mode << "\" not supported!" << endl;
				 global_log->error() << " Has to be  \"Homogeneous\"  or \"Cluster\"  " << endl;
				 exit(1);
				 }
				 generator->readPhaseSpaceHeader(_domain, timestepLength);

				 _inputReader = (OneCLJGenerator*) generator;
				 */
			}
			/* TODO: Check this part of the code */
			if (_LJCutoffRadius == 0.0)
				_LJCutoffRadius = _cutoffRadius;
			_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
		}

		if (inp.changecurrentnode("algorithm")) {

			string datastructype;
			if (inp.getNodeValue("datastructure@type", datastructype)) {
				global_log->info() << "datastructure to use:\t"
						<< datastructype << endl;
				if (datastructype == "LinkedCells") {
					int cellsInCutoffRadius = 1;
					inp.getNodeValue("datastructure/cellsInCutoffRadius",
							cellsInCutoffRadius);
					global_log->info()
							<< "LinkedCells cells in cutoff radius:\t"
							<< cellsInCutoffRadius << endl;
					double bBoxMin[3];
					double bBoxMax[3];
					global_log->debug() << "Box dimensions:" << endl;
					for (int d = 0; d < 3; d++) {
						bBoxMin[d] = _domainDecomposition->getBoundingBoxMin(d,
								_domain);
						bBoxMax[d] = _domainDecomposition->getBoundingBoxMax(d,
								_domain);
						global_log->debug() << "[" << d << "] " << bBoxMin[d]
								<< " - " << bBoxMax[d] << endl;
					}
					if (_LJCutoffRadius == 0.0)
						_LJCutoffRadius = _cutoffRadius;
					global_log->debug() << "Cutoff: " << _cutoffRadius
							<< "\nCellsInCutoff: " << cellsInCutoffRadius
							<< endl;
					_moleculeContainer = new LinkedCells(bBoxMin, bBoxMax,
							_cutoffRadius, _LJCutoffRadius,
							_tersoffCutoffRadius, cellsInCutoffRadius);

				} else if (datastructype == "AdaptiveSubCells") {
					double bBoxMin[3];
					double bBoxMax[3];
					for (int i = 0; i < 3; i++) {
						bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i,
								_domain);
						bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i,
								_domain);
					}
					//creates a new Adaptive SubCells datastructure
					if (_LJCutoffRadius == 0.0)
						_LJCutoffRadius = _cutoffRadius;
					_moleculeContainer = new AdaptiveSubCells(bBoxMin, bBoxMax,
							_cutoffRadius, _LJCutoffRadius,
							_tersoffCutoffRadius);
				} else {
					global_log->error() << "UNKOWN DATASTRUCTURE: "
							<< datastructype << endl;
					exit(1);
				}
			}
			inp.changecurrentnode("..");
		} // algo-section
		else {
			global_log->info() << "XML config file " << inputfilename
					<< ": no algo section" << endl;
		}
		if (inp.changecurrentnode("output")) {
			string outputPathAndPrefix;
			if (inp.getNodeValue("Resultwriter", outputPathAndPrefix)) {
				// TODO: add the output frequency to the xml!
				_outputPlugins.push_back(new ResultWriter(1,
						outputPathAndPrefix));
				global_log->debug() << "ResultWriter '" << outputPathAndPrefix
						<< "'.\n";
			}
			inp.changecurrentnode("..");
		} // output-section
		else {
			global_log->info() << "XML config file " << inputfilename
					<< ": no output section" << endl;
		}
		inp.changecurrentnode("..");
	} // simulation-section
	else {
		global_log->error() << "XML config file " << inputfilename
				<< ": no simulation section" << endl;
	}

	// read particle data
	unsigned long maxid = _inputReader->readPhaseSpace(_moleculeContainer,
			&_lmu, _domain, _domainDecomposition);

	if (this->_LJCutoffRadius == 0.0)
		_LJCutoffRadius = this->_cutoffRadius;

	global_log->info() << " reading in components from xml ..." << endl;
	/* TODO: save components into domain:
	 * take care of already existing components
	 * take care about parameter streams
	 */
	std::vector<Component>& dcomponents = _domain->getComponents();
	std::vector<Component> components;
	xmlreader.getComponents(components);
	// 	for( unsigned int j = 0; j < components.size(); j++ ) {
	// 		global_log->info() << components[j] << endl;
	// 	}
	long numComponents = dcomponents.size();
	global_log->info() << " number of components from input: " << numComponents
			<< endl;
	global_log->info() << " number of components in xml (not implemented): "
			<< components.size() << endl;

	// TODO: Check this...
	// 	_domain->initParameterStreams( _cutoffRadius, _LJCutoffRadius );
	_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

	// test new Decomposition
	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	unsigned idi = _lmu.size();
	unsigned j = 0;
	std::list<ChemicalPotential>::iterator cpit;
	for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
		cpit->setIncrement(idi);
		double tmp_molecularMass =
				_domain->getComponents()[cpit->getComponentID()].m();
		cpit->setSystem(_domain->getGlobalLength(0),
				_domain->getGlobalLength(1), _domain->getGlobalLength(2),
				tmp_molecularMass);
		cpit->setGlobalN(_domain->N(cpit->getComponentID()));
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

	double timestepLength;
	unsigned cosetid = 0;
        bool widom = false;

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
				string generatorName; // name of the library to load
				string inputFile; // name of the input file for the generator

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
			if (token=="DomainDecomposition") {
				_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
			}
			else if(token=="KDDecomposition") {
				_domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, 1.0, 100);
			}
			else if(token=="KDDecomposition2") {
				int updateFrequency = 100;
				int fullSearchThreshold = 3;
				string line;
				getline(inputfilestream, line);
				stringstream lineStream(line);
				lineStream >> updateFrequency >> fullSearchThreshold;
				_domainDecomposition = (DomainDecompBase*) new KDDecomposition2(_cutoffRadius, _domain, updateFrequency, fullSearchThreshold);
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
				_particleContainerType = LINKED_CELL;
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
				_moleculeContainer = new LinkedCells(bBoxMin, bBoxMax,
						_cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius,
						cellsInCutoffRadius);
				cout << "Cutoff called with" << _cutoffRadius << " "
						<< _LJCutoffRadius << " " << _tersoffCutoffRadius
						<< endl;
				for (int i = 0; i < 3; i++)
					cout << bBoxMin[i] << " ";
				cout << endl;
				for (int i = 0; i < 3; i++)
					cout << bBoxMax[i] << " ";
				cout << endl;
				cout << _tersoffCutoffRadius << " " << cellsInCutoffRadius
						<< endl;
			} else if (token == "AdaptiveSubCells") {
				_particleContainerType = ADAPTIVE_LINKED_CELL;
				double bBoxMin[3];
				double bBoxMax[3];
				for (int i = 0; i < 3; i++) {
					bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i,
							_domain);
					bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i,
							_domain);
				}
				//creates a new Adaptive SubCells datastructure
				if (_LJCutoffRadius == 0.0)
					_LJCutoffRadius = _cutoffRadius;
				_moleculeContainer = new AdaptiveSubCells(bBoxMin, bBoxMax,
						_cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius);
			} else {
				global_log->error() << "UNKOWN DATASTRUCTURE: " << token
						<< endl;
				exit(1);
			}
		} else if (token == "output") {
			inputfilestream >> token;
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
						outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "XyzWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "CheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new CheckpointWriter(writeFrequency,
						outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "CheckpointWriter " << writeFrequency
						<< " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "PovWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new PovWriter(writeFrequency,
						outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "POVWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "DecompWriter") {
				unsigned long writeFrequency;
				string mode;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> mode
						>> outputPathAndPrefix;
				_outputPlugins.push_back(new DecompWriter(writeFrequency, mode,
						outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "DecompWriter " << writeFrequency
						<< " '" << outputPathAndPrefix << "'.\n";
			} else if ((token == "VisittWriter") || (token == "VISWriter")) {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new VISWriter(writeFrequency,
						outputPathAndPrefix, _numberOfTimesteps, true));
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
				Log::global_log->error() << std::endl << "VKT-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			} else if (token == "VTKGridWriter") {
#ifdef VTK
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;

				if (_particleContainerType == LINKED_CELL) {
					_outputPlugins.push_back(new VTKGridWriter(writeFrequency,
							outputPathAndPrefix,
							static_cast<LinkedCells&> (*_moleculeContainer)));
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
				Log::global_log->error() << std::endl << "VKT-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			} else if (token == "StatisticsWriter") {
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;

				if (_particleContainerType == LINKED_CELL) {
					_outputPlugins.push_back(new StatisticsWriter(
							writeFrequency, outputPathAndPrefix,
							static_cast<LinkedCells&> (*_moleculeContainer)));
					global_log->debug() << "StatisticsWriter "
							<< writeFrequency << " '" << outputPathAndPrefix
							<< "'.\n";
				} else {
					global_log->warning()
							<< "StatisticsWriter only supported with LinkedCells!"
							<< std::endl;
					global_log->warning()
							<< "Generating no statistics output for the grid!"
							<< std::endl;
				}
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
		} else if (token == "profile") {
			unsigned xun, yun, zun;
			inputfilestream >> xun >> yun >> zun;
			_domain->setupProfile(xun, yun, zun);
			_doRecordProfile = true;
		} else if (token == "profileRecordingTimesteps") { /* TODO: subotion of profile */
			inputfilestream >> _profileRecordingTimesteps;
		} else if (token == "profileOutputTimesteps") { /* TODO: subotion of profile */
			inputfilestream >> _profileOutputTimesteps;
		} else if (token == "RDF") {
			double interval;
			unsigned bins;
			inputfilestream >> interval >> bins;
			if (_domain->getComponents().size() <= 0) {
				global_log->error()
						<< "PhaseSpaceFile-Specifiation has to occur befor RDF-Token!"
						<< endl;
				exit(-1);
			}
			_rdf = new RDF(interval, bins, _domain->getComponents());
			//_domain->setupRDF(interval, bins);
		} else if (token == "RDFOutputTimesteps") { /* TODO: subotion of RDF */
			unsigned int RDFOutputTimesteps;
			inputfilestream >> RDFOutputTimesteps;
			_rdf->setOutputTimestep(RDFOutputTimesteps);
		} else if (token == "RDFOutputPrefix") { /* TODO: subotion of RDF */
			std::string RDFOutputPrefix;
			inputfilestream >> RDFOutputPrefix;
			_rdf->setOutputPrefix(RDFOutputPrefix);
		} else if (token == "profiledComponent") { /* TODO: subotion of profile, check if required to enable output in general */
			unsigned cid;
			inputfilestream >> cid;
			cid--;
			_domain->considerComponentInProfile(cid);
		} else if (token == "profileOutputPrefix") { /* TODO: subotion of profile */
			inputfilestream >> _profileOutputPrefix;
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
			tmu.setInstances(intest);
			if (controlVolume)
				tmu.setControlVolume(x0, y0, z0, x1, y1, z1);
			global_log->info() << setprecision(6) << "chemical Potential "
					<< imu << " component " << icid + 1 << " (internally "
					<< icid << ") conduct " << intest << " tests every "
					<< instep << " steps: ";
			_lmu.push_back(tmu);
			global_log->info() << " pushed back." << endl;
		} else if (token == "planckConstant") {
			inputfilestream >> h;
		} else if(token == "Widom") {
                        widom = true;
                } else if (token == "NVE") {
			/* TODO: Documentation, what it does (no "Enerstat" at the moment) */
			_domain->thermostatOff();
		} else if (token == "initCanonical") {
			inputfilestream >> _initCanonical;
		} else if (token == "initGrandCanonical") { /* suboption of chemical potential */
			inputfilestream >> _initGrandCanonical;
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
		double tmp_molecularMass =
				_domain->getComponents()[cpit->getComponentID()].m();
		cpit->setSystem(_domain->getGlobalLength(0),
				_domain->getGlobalLength(1), _domain->getGlobalLength(2),
				tmp_molecularMass);
		cpit->setGlobalN(_domain->N(cpit->getComponentID()));
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

	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();
	global_log->info() << "Updating domain decomposition" << endl;
	updateParticleContainerAndDecomposition();
	global_log->info() << "Performing inital force calculation" << endl;
	_moleculeContainer->traversePairs(_particlePairsHandler);
	// TODO:
	// here we have to call calcFM() manually, otherwise force and moment are not
	// updated inside the molecule (actually this is done in upd_postF)
	// or should we better call the integrator->eventForcesCalculated?
	for (Molecule* tM = _moleculeContainer->begin(); tM
			!= _moleculeContainer->end(); tM = _moleculeContainer->next()) {
		tM->calcFM();
	}

	// clear halo
	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();

	// initialize the radial distribution function
	// TODO this call should not be neccessary!
	//	if (_rdf != NULL)
	//		_rdf->reset();

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
				Ttar =
						_domain->severalThermostats() ? _domain->getTargetTemperature(
								1)
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

#ifdef STEEREO
	_steer = initSteereo (_domainDecomposition->getRank(), _domainDecomposition->getNumProcs());
#ifdef STEEREO_COUPLING
	_coupling = initCoupling(_steer, (long*) &_simstep);
#endif  /* STEEREO_COUPLING */
	registerSteereoCommands (_steer, this);
	startListeningSteereo (_steer);
#endif  /* STEEREO */

	//	if ((_initSimulation > _initStatistics) && this->_rdf != NULL) {
	//		this->_rdf->tickRDF();
	//		this->_particlePairsHandler->setRDF(_rdf);
	//	}

	global_log->info() << "System initialised\n" << endl;
	global_log->info() << "System contains "
			<< _domain->getglobalNumMolecules() << " molecules." << endl;
}

void Simulation::simulate() {

	Molecule* tM;
	global_log->info() << "Started simulation" << endl;

	// added by tijana because of getEnergy(...) in LinkedCells
	_moleculeContainer->setPairsHandler(_particlePairsHandler);

	_simstep = 0;


	// (universal) constant acceleration (number of) timesteps
	unsigned uCAT = _pressureGradient->getUCAT();
	_initSimulation = (unsigned long) (_domain->getCurrentTime()
			/ _integrator->getTimestepLength());

	/* demonstration for the usage of the new ensemble class */
	CanonicalEnsemble ensemble(_moleculeContainer, &(_domain->getComponents()));
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
	// all timers except the ioTimer messure inside the main loop
	Timer loopTimer;
	Timer decompositionTimer;
	Timer perStepIoTimer;
	Timer ioTimer;

#if defined(STEEREO) && defined(STEEREO_COUPLING)
	_simstep = _initSimulation;
	//SteereoLogger::setOutputLevel (4);
	_coupling->waitForConnection();
#endif

	loopTimer.start();
	for (_simstep = _initSimulation; _simstep <= _numberOfTimesteps; _simstep++) {
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
		global_log->debug() << "timestep " << _simstep << endl;

		_integrator->eventNewTimestep(_moleculeContainer, _domain);

#if defined(STEEREO) && defined(STEEREO_COUPLING)
		_steer -> processQueue (1);
		global_log->debug() << "molecules in simulation: " << _moleculeContainer->getNumberOfParticles() << std::endl;
#endif
		// activate RDF sampling
		if ((_simstep >= this->_initStatistics) && this->_rdf != NULL) {
			this->_rdf->tickRDF();
			this->_particlePairsHandler->setRDF(_rdf);
			this->_rdf->accumulateNumberOfMolecules(_domain->getComponents());
		}

		// ensure that all Particles are in the right cells and exchange Particles
		global_log->debug() << "Updating container and decomposition" << endl;
		loopTimer.stop();
		decompositionTimer.start();
		updateParticleContainerAndDecomposition();
		decompositionTimer.stop();
		loopTimer.start();

		// Force calculation
		global_log->debug() << "Traversing pairs" << endl;
		//cout<<"here somehow"<<endl;
		_moleculeContainer->traversePairs(_particlePairsHandler);

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
					_moleculeContainer->grandcanonicalStep(&(*cpit),
							_domain->getGlobalCurrentTemperature(), this->_domain);
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

		// clear halo
		global_log->debug() << "Delete outer particles" << endl;
		_moleculeContainer->deleteOuterParticles();

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

		/*
		 * radial distribution function
		 */
		if (_simstep >= _initStatistics) {
			if (this->_lmu.size() == 0) {
				this->_domain->record_cv();
			}
			//			if (this->_rdf != NULL) {
			//				this->_rdf->tickRDF();
			//
			//			}
		}

		if (_zoscillation) {
			global_log->debug() << "alert z-oscillators" << endl;
			_integrator->zOscillation(_zoscillator, _moleculeContainer);
		}

		// Inform the integrator about the calculated forces
		global_log->debug() << "Inform the integrator" << endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);

		// calculate the global macroscopic values from the local values
		global_log->debug() << "Calculate macroscopic values" << endl;
		_domain->calculateGlobalValues(_domainDecomposition,
				_moleculeContainer, (!(_simstep
						% _collectThermostatDirectedVelocity)), Tfactor(
								_simstep));

		// scale velocity and angular momentum
		if (!_domain->NVE()) {
			global_log->debug() << "Velocity scaling" << endl;
			if (_domain->severalThermostats()) {
				for (tM = _moleculeContainer->begin(); tM
				!= _moleculeContainer->end(); tM
				= _moleculeContainer->next()) {
					int thermostat = _domain->getThermostat(tM->componentid());
					if (0 >= thermostat)
						continue;
					if (_domain->thermostatIsUndirected(thermostat)) {
						/* TODO: thermostat */
						tM->scale_v(_domain->getGlobalBetaTrans(thermostat),
								_domain->getThermostatDirectedVelocity(
										thermostat, 0),
										_domain->getThermostatDirectedVelocity(
												thermostat, 1),
												_domain->getThermostatDirectedVelocity(
														thermostat, 2));
					} else {
						tM->scale_v(_domain->getGlobalBetaTrans(thermostat));
					}
					tM->scale_D(_domain->getGlobalBetaRot(thermostat));
				}
			} else {
				for (tM = _moleculeContainer->begin(); tM
				!= _moleculeContainer->end(); tM
				= _moleculeContainer->next()) {
					tM->scale_v(_domain->getGlobalBetaTrans());
					tM->scale_D(_domain->getGlobalBetaRot());
				}
			}
		}

		_domain->advanceTime(_integrator->getTimestepLength());
#ifdef STEEREO
		_steer -> processQueue (0);
#endif

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
		perStepIoTimer.stop();
		loopTimer.start();

	}
	loopTimer.stop();
	/***************************************************************************/
	/* END MAIN LOOP                                                           */
	/*****************************//**********************************************/

	ioTimer.start();
	/* write final checkpoint */
	string cpfile(_outputPrefix + ".restart.xdr");
	_domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);
	// finish output
	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter
	!= _outputPlugins.end(); outputIter++) {
		(*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition,
				_domain);
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


}

void Simulation::output(unsigned long simstep) {

	int ownrank = _domainDecomposition->getRank();

	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter
			!= _outputPlugins.end(); outputIter++) {
		(*outputIter)->doOutput(_moleculeContainer, _domainDecomposition,
				_domain, simstep, &(_lmu));
	}

	if (_rdf != NULL) {
		_rdf->doOutput(_domainDecomposition, _domain, simstep);
	}

	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep
			% _profileRecordingTimesteps)) {
		_domain->recordProfile(_moleculeContainer);
	}
	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep
			% _profileOutputTimesteps)) {
		_domain->collectProfile(_domainDecomposition);
		if (!ownrank) {
			ostringstream osstrm;
			osstrm << _profileOutputPrefix << ".";
			osstrm.fill('0');
			osstrm.width(9);
			osstrm << right << simstep;
			_domain->outputProfile(osstrm.str().c_str());
			osstrm.str("");
			osstrm.clear();
		}
		_domain->resetProfile();
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
#ifdef ENABLE_MPI
	delete _domainDecomposition;
	_domainDecomposition = NULL;
#endif
}

void Simulation::updateParticleContainerAndDecomposition() {

	// The particles have moved, so the neighbourhood relations have
	// changed and have to be adjusted
	_moleculeContainer->update();
	//_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain->getComponents(), _domain);
	_domainDecomposition->balanceAndExchange(true, _moleculeContainer,
			_domain->getComponents(), _domain);
	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateMoleculeCaches();
}

/* FIXME: we shoud provide a more general way of doing this */
double Simulation::Tfactor(unsigned long simstep) {
	double xi = (double) (simstep - _initSimulation) / (double) (_initCanonical
			- _initSimulation);
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
	_moleculeContainer = NULL;
	_integrator = NULL;
	_inputReader = NULL;

#ifndef ENABLE_MPI
	global_log->info() << "Initializing the alibi domain decomposition ... "
			<< endl;
	_domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
	global_log->info() << "Initialization done" << endl;
#endif

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
	_profileRecordingTimesteps = 7;
	_profileOutputTimesteps = 12500;
	_profileOutputPrefix = "out";
	_collectThermostatDirectedVelocity = 100;
	_zoscillation = false;
	_zoscillator = 512;
	_initCanonical = 5000;
	_initGrandCanonical = 10000000;
	_initStatistics = 20000;
	h = 0.0;

	_pressureGradient = new PressureGradient(ownrank);
	global_log->info() << "Constructing domain ..." << endl;
	_domain = new Domain(ownrank, this->_pressureGradient);
	global_log->info() << "Domain construction done." << endl;
	_particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);
}

void Simulation::mkTcTS(int argc, char** argv) {
	_particleContainerType = LINKED_CELL;

	TcTS(argc, argv, this->_domain, &(this->_domainDecomposition),
			&(this->_integrator), &(this->_moleculeContainer),
			&(this->_outputPlugins), this->_rdf, this);

	if (this->_LJCutoffRadius == 0.0)
		_LJCutoffRadius = this->_cutoffRadius;
	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
	_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

	/*
	 unsigned idi = _lmu.size();
	 unsigned j = 0;
	 std::list<ChemicalPotential>::iterator cpit;
	 double Tcur = _domain->getCurrentTemperature(0);
	 for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
	 cpit->setIncrement(idi);
	 double tmp_molecularMass = _domain->getComponents()[cpit->getComponentID()].m();
	 cpit->setSystem(_domain->getGlobalLength(0), _domain->getGlobalLength(1), _domain->getGlobalLength(2),
	 tmp_molecularMass);
	 cpit->setGlobalN(_domain->N(cpit->getComponentID()));
	 cpit->setNextID(j + (int) (1.001 * (256 + maxid)));
	 cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0), _moleculeContainer->getBoundingBoxMax(0),
	 _moleculeContainer->getBoundingBoxMin(1), _moleculeContainer->getBoundingBoxMax(1),
	 _moleculeContainer->getBoundingBoxMin(2), _moleculeContainer->getBoundingBoxMax(2));
	 double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(0);
	 if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar)) Tcur = Ttar;
	 cpit->submitTemperature(Tcur);
	 if (h != 0.0) cpit->setPlanckConstant(h);
	 j++;
	 }
	 */
}

