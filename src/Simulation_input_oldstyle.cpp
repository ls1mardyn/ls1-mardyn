#include "Simulation.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>

#include "Common.h"
#include "Domain.h"
#include "particleContainer/LinkedCells.h"
#include "parallel/DomainDecompBase.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#else
#include "parallel/DomainDecompBase.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/adapter/FlopCounter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "integrators/ExplicitEuler.h"
#include "molecules/Wall.h"
#include "molecules/Mirror.h"

#include "io/io.h"
#include "io/GeneratorFactory.h"
#include "io/RDF.h"

#include "ensemble/GrandCanonical.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/PressureGradient.h"

#include "thermostats/VelocityScalingThermostat.h"
#include "thermostats/TemperatureControl.h"

#include "utils/Logger.h"
#include "utils/FileUtils.h"

#include "longRange/LongRangeCorrection.h"
#include "longRange/Homogeneous.h"
#include "longRange/Planar.h"
#include "particleContainer/adapter/VectorizationTuner.h"

#include "NEMD/NEMD.h"
#include "NEMD/DriftControl.h"
#include "NEMD/DistControl.h"
#include "NEMD/RegionSampling.h"
#include "NEMD/DensityControl.h"
#include "NEMD/ParticleTracker.h"

using namespace std;

void Simulation::initConfigOldstyle(const string& inputfilename) {
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &ownrank));
#endif

	global_log->info()
			<< "\n# Please address your questions and suggestions to the ls1 mardyn contact point:\n# \n# E-mail: contact@ls1-mardyn.de\n# \n# Phone: +49 631 205 3227\n# Fax: +49 631 205 3835\n# University of Kaiserslautern\n# Computational Molecular Engineering\n# Erwin-Schroedinger-Str. 44\n# D-67663 Kaiserslautern, Germany\n# \n# http://www.ls1-mardyn.de/\n\n";

	global_log->info() << "init oldstyle config file: " << inputfilename << endl;

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
	map<unsigned, double*> cavity_grid;

	// The first line of the config file has to contain the token "MDProjectConfig"
	inputfilestream >> token;
	if ((token != "mardynconfig") && (token != "MDProjectConfig")) {
		global_log->error() << "Not a mardynconfig file! First token: " << token << endl;
		exit(1);
	}

	while (inputfilestream && !inputfilestream.eof()) {
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
				global_log->info() << "phaseSpaceFileFormat is Generator!" << endl;
				string generatorName;  // name of the library to load
				string inputFile;  // name of the input file for the generator

				string line;
				getline(inputfilestream, line);
				stringstream lineStream(line);
				lineStream >> generatorName >> inputFile;
				_inputReader = GeneratorFactory::loadGenerator(generatorName, inputFile);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else if (phaseSpaceFileFormat == "Binary") {
				_inputReader = (InputBase*) new BinaryReader();
				string token;
				inputfilestream >> token;
				_inputReader->setPhaseSpaceHeaderFile(token);
				inputfilestream >> token;
				_inputReader->setPhaseSpaceFile(token);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else if (phaseSpaceFileFormat == "MPI-IO") {
			#ifdef ENABLE_MPI
				_inputReader = (InputBase*) new MPI_IOReader();
				string token;
				inputfilestream >> token;
				_inputReader->setPhaseSpaceHeaderFile(token);
				inputfilestream >> token;
				_inputReader->setPhaseSpaceFile(token);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			#else
				global_log->error() << "MPI not enabled! Program exit..." << endl;
				Simulation::exit(1);
			#endif
			} else {
				global_log->error() << "Don't recognize phasespaceFile reader " << phaseSpaceFileFormat << endl;
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
			string line;
			getline(inputfilestream, line);
			stringstream lineStream(line);
			lineStream >> token;
#else
			string line;
			getline(inputfilestream, line);
			stringstream lineStream(line);
			lineStream >> token;
			if (token == "DomainDecomposition") {
				// default DomainDecomposition is already set in initialize();
				//_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
				if (line.find("direct") != string::npos) {
					dynamic_cast<DomainDecompMPIBase*>(_domainDecomposition)->setCommunicationScheme("direct");
				}
			} else if (token == "KDDecomposition") {
				delete _domainDecomposition;
				int updateFrequency = 100;
				int fullSearchThreshold = 3;
				lineStream >> updateFrequency >> fullSearchThreshold;
				bool hetero = false, cutsmaller = false, forceRatio = false;
				if (line.find("hetero") != string::npos) {
					hetero = true;
				}
				if (line.find("cutSmaller") != string::npos) {
					cutsmaller = true;  // allow domain to be split not only along biggest side
				}
				if (line.find("forceRatio") != string::npos) {
					forceRatio = true; // never do a search for the best Partitioning, always force the ratio (should be similar to fullSeachThreshold=0)
				}
				_domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, updateFrequency,
						fullSearchThreshold, hetero, cutsmaller, forceRatio);
				double bBoxMin[3];
				double bBoxMax[3];
				_domainDecomposition->getBoundingBoxMinMax(_domain, bBoxMin, bBoxMax);
				if (_moleculeContainer != NULL) {
					_moleculeContainer->rebuild(bBoxMin, bBoxMax);
				}

				if (line.find("direct") != string::npos) {
					dynamic_cast<DomainDecompMPIBase*>(_domainDecomposition)->setCommunicationScheme("direct"); // never do a search for the best Partitioning, always force the ratio (should be similar to fullSeachThreshold=0)
				}
			}
#endif
		} else if (token == "datastructure") {

			if (_domainDecomposition == nullptr) {
				global_log->error()
						<< "_domainDecomposition is NULL! Probably you compiled for MPI, but didn't specify line \"parallelization\" before line \"datastructure\"!"
						<< endl;
				exit(1);
			}

			inputfilestream >> token;
			if (token == "LinkedCells") {
				int cellsInCutoffRadius;
				inputfilestream >> cellsInCutoffRadius;
//				cellsInCutoffRadius - Not used anymore. Read it for backwards compatibility with
				// the to-be-removed .cfg files
				double bBoxMin[3];
				double bBoxMax[3];
				for (int i = 0; i < 3; i++) {
					bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
					bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
				}
				if (this->_LJCutoffRadius == 0.0)
					_LJCutoffRadius = this->_cutoffRadius;
				_moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, _cutoffRadius);
			} else if (token == "AdaptiveSubCells") {
				global_log->error() << "AdaptiveSubCells no longer supported." << std::endl;
				global_simulation->exit(-1);
			} else {
				global_log->error() << "UNKNOWN DATASTRUCTURE: " << token << endl;
				exit(1);
			}
		} else if (token == "output") {
			token.clear();
			inputfilestream >> token;
			if (token == "ResultWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new ResultWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "ResultWriter '" << outputPathAndPrefix << "'.\n";
			} else if (token == "XyzWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new XyzWriter(writeFrequency, outputPathAndPrefix, true));
				global_log->debug() << "XyzWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "CavityWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new CavityWriter(writeFrequency, outputPathAndPrefix, true));
				global_log->debug() << "CavityWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "CheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new CheckpointWriter(writeFrequency, outputPathAndPrefix, true));
				if (writeFrequency == 0) {
					global_log->error() << "Write frequency must be a positive nonzero integer, but is "
							<< writeFrequency << endl;
					exit(-1);
				}
				global_log->debug() << "CheckpointWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "BinaryCheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new BinaryCheckpointWriter(writeFrequency, outputPathAndPrefix, true));
				global_log->debug() << "BinaryCheckpointWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
#ifdef ENABLE_MPI
			} else if (token == "MPI-IOCheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new MPI_IOCheckpointWriter(writeFrequency, outputPathAndPrefix, true));
				global_log->debug() << "MPI-IOCheckpointWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
#endif
			} else if (token == "PovWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new PovWriter(writeFrequency, outputPathAndPrefix, true));
				global_log->debug() << "POVWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "DecompWriter") {
				unsigned long writeFrequency;
				string mode;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> mode >> outputPathAndPrefix;
				_outputPlugins.push_back(new DecompWriter(writeFrequency, mode, outputPathAndPrefix, true));
				global_log->debug() << "DecompWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if ((token == "VisittWriter") || (token == "VISWriter")) {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new VISWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "VISWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "MmspdBinWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new MmspdBinWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "MmspdBinWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "MmpldWriter") {
				std::string strSpheres;
				uint64_t startTimestep;
				uint64_t writeFrequency;
				uint64_t stopTimestep;
				uint64_t numFramesPerFile;
				std::string strWriteControl;
				std::string outputPathAndPrefix;
				std::string strInitSphereData;
				std::string strSphereDataFilename = "unknown";
				uint8_t bInitSphereData = ISD_USE_DEFAULT;
				inputfilestream >> strSpheres >> strWriteControl >> outputPathAndPrefix >> strInitSphereData;

				// tokenize write control parameters
				std::vector<string> fields;
				fields = split( fields, strWriteControl, ":", split_type::no_empties );
				startTimestep    = atoi( fields.at(0).c_str() );
				writeFrequency   = atoi( fields.at(1).c_str() );
				stopTimestep     = atoi( fields.at(2).c_str() );
				numFramesPerFile = atoi( fields.at(3).c_str() );

				if("file" == strInitSphereData)
				{
					inputfilestream >> strSphereDataFilename;
					bInitSphereData = ISD_READ_FROM_FILE;
				}
				else if("default" == strInitSphereData)
					bInitSphereData = ISD_USE_DEFAULT;
				else
				{
					global_log->error() << "MmpldWriter: wrong statement, expected default|file. Program exit... " << endl;
					exit(-1);
				}

				MmpldWriter* mmpldWriter = NULL;
				if("simple" == strSpheres)
					mmpldWriter = new MmpldWriterSimpleSphere(startTimestep, writeFrequency, stopTimestep, numFramesPerFile, outputPathAndPrefix);
				else if("multi" == strSpheres)
					mmpldWriter = new MmpldWriterMultiSphere (startTimestep, writeFrequency, stopTimestep, numFramesPerFile, outputPathAndPrefix);
				else
				{
					global_log->error() << "MmpldWriter: wrong statement, expected simple|multi. Program exit... " << endl;
					exit(-1);
				}
				if(NULL != mmpldWriter)
				{
					mmpldWriter->SetInitSphereDataParameters(bInitSphereData, strSphereDataFilename);
					_outputPlugins.push_back(mmpldWriter);
				}
				global_log->debug() << "MmpldWriter " << strSpheres << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "VTKWriter") {
#ifdef VTK

				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new VTKMoleculeWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "VTKWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
#else
				Log::global_log->warning() << std::endl << "VTK-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			} else if (token == "VTKGridWriter") {
#ifdef VTK
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;

				if (dynamic_cast<LinkedCells*>(_moleculeContainer)) {
					_outputPlugins.push_back(new VTKGridWriter(writeFrequency, outputPathAndPrefix));
					global_log->debug() << "VTKGridWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
				} else {
					global_log->warning() << "VTKGridWriter only supported with LinkedCells!" << std::endl;
					global_log->warning() << "Generating no VTK output for the grid!" << std::endl;
				}
#else
				Log::global_log->error() << std::endl << "VTK-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			}
			// by Stefan Becker <stefan.becker@mv.uni-kl.de>
			// output for the MegaMol Simple Particle Data File Format (*.mmspd)
			else if (token == "MmspdWriter") {
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
				_outputPlugins.push_back(new MPICheckpointWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "MPICheckpointWriter " << writeFrequency << " '" << outputPathAndPrefix
						<< "'.\n";
			} else if (token == "GammaWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new GammaWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "GammaWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "VectorizationTuner") {
				string outputPathAndPrefix;
				unsigned int minMoleculeCnt, maxMoleculeCnt;
				int type;
				MoleculeCntIncreaseTypeEnum moleculeCntIncreaseType;
				inputfilestream >> outputPathAndPrefix >> minMoleculeCnt >> maxMoleculeCnt >> type;
				moleculeCntIncreaseType = static_cast<MoleculeCntIncreaseTypeEnum>(type);
				_outputPlugins.push_back(
						new VectorizationTuner(outputPathAndPrefix, minMoleculeCnt, maxMoleculeCnt,
								moleculeCntIncreaseType, _cutoffRadius, _LJCutoffRadius, &_cellProcessor));
				global_log->debug() << "VectorizationTuner " << outputPathAndPrefix << "'.\n";
			} else {
				global_log->warning() << "Unknown output plugin " << token << endl;
			}
		} else if (token == "accelerate") {
			cosetid++;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "comp") {
				global_log->error() << "Expected 'comp' instead of '" << token << "'.\n";
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
				global_log->error() << "Expected 'towards' instead of '" << token << "'.\n";
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
				global_log->error() << "Expected 'within' instead of '" << token << "'.\n";
				exit(1);
			}
			double tau;
			inputfilestream >> tau;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "from") {
				global_log->error() << "Expected 'from' instead of '" << token << "'.\n";
				exit(1);
			}
			double ainit[3];
			for (unsigned d = 0; d < 3; d++)
				inputfilestream >> ainit[3];
			if (timestepLength == 0.0) {
				global_log->error() << "timestep missing." << endl;
				exit(1);
			}
			_pressureGradient->specifyComponentSet(cosetid, dir, tau, ainit, timestepLength);
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

		}
		//by Stefan Becker
		else if (token == "yOffset") {
			double yOffset;
			inputfilestream >> yOffset;
			_domain->sYOffset(yOffset);
		} else if (token == "profileVirial") {
			_doRecordVirialProfile = true;
		} else if (token == "profileRecordingTimesteps") { /* TODO: suboption of profile */
			inputfilestream >> _profileRecordingTimesteps;
		} else if (token == "profileOutputTimesteps") { /* TODO: suboption of profile */
			inputfilestream >> _profileOutputTimesteps;
		} else if (token == "RDF") {
			double interval;
			unsigned bins;
			inputfilestream >> interval >> bins;
			if (global_simulation->getEnsemble()->getComponents()->size() <= 0) {
				global_log->error() << "PhaseSpaceFile-Specification has to occur before RDF-Token!" << endl;
				exit(-1);
			}
			_rdf = new RDF(interval, bins, global_simulation->getEnsemble()->getComponents());
			_outputPlugins.push_back(_rdf);
		} else if (token == "RDFOutputTimesteps") { /* TODO: suboption of RDF */
			unsigned int RDFOutputTimesteps;
			inputfilestream >> RDFOutputTimesteps;
			_rdf->setOutputTimestep(RDFOutputTimesteps);
		} else if (token == "RDFOutputPrefix") { /* TODO: suboption of RDF */
			std::string RDFOutputPrefix;
			inputfilestream >> RDFOutputPrefix;
			_rdf->setOutputPrefix(RDFOutputPrefix);
		} else if (token == "profiledComponent") { /* TODO: suboption of profile, check if required to enable output in general */
			unsigned cid;
			inputfilestream >> cid;
			cid--;
			_domain->considerComponentInProfile(cid);
		}
		// by Stefan Becker <stefan.becker@mv.uni-kl.de>, token determining the coordinate system of the density profile
		else if (token == "SessileDrop") {
			this->_domain->sesDrop();
		} else if (token == "profileOutputPrefix") { /* TODO: suboption of profile */
			inputfilestream >> _profileOutputPrefix;
		} else if (token == "collectThermostatDirectedVelocity") { /* suboption of the thermostat replace with direct thermostat */
			inputfilestream >> _collectThermostatDirectedVelocity;
		}

		else if (token == "thermostat") {
			inputfilestream >> _thermostatType;
			if (_thermostatType == ANDERSEN_THERMOSTAT) {
				inputfilestream >> _nuAndersen;
				global_log->info() << "Using the Andersen Thermostat with nu = " << _nuAndersen << "\n";
			}
		} else if (token == "zOscillator") {
			global_log->error() << "zOscillator was used for the Tersoff potential, which is no longer supported."
					<< std::endl;
			global_simulation->exit(-1);
		}
		// by Stefan Becker
		else if (token == "AlignCentre") {
			_doAlignCentre = true;
			inputfilestream >> _alignmentInterval >> _alignmentCorrection;
		} else if (token == "ComponentForYShift") {
			_componentSpecificAlignment = true;
			unsigned cidMin, cidMax;
			inputfilestream >> cidMin >> cidMax;
			cidMin--; // since internally the component number is reduced by one, i.e. cid == 1 in the input file corresponds to the internal cid == 0
			cidMax--;
			_domain->considerComponentForYShift(cidMin, cidMax);
		}
		// chemicalPotential <mu> component <cid> [control <x0> <y0> <z0>
		// to <x1> <y1> <z1>] conduct <ntest> tests every <nstep> steps
		else if (token == "chemicalPotential") {
			double imu;
			inputfilestream >> imu;
			inputfilestream >> token;
			if (token != "component") {
				global_log->error() << "Expected 'component' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
						<< endl;
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
					global_log->error() << "Expected 'to' instead of '" << token << "'.\n";
					global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
							<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
							<< "conduct <ntest> tests every <nstep> steps" << endl;
					exit(1);
				}
				inputfilestream >> x1 >> y1 >> z1;
				inputfilestream >> token;
			}
			if (token != "conduct") {
				global_log->error() << "Expected 'conduct' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
						<< endl;
				exit(1);
			}
			unsigned intest;
			inputfilestream >> intest;
			inputfilestream >> token;
			if (token != "tests") {
				global_log->error() << "Expected 'tests' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
						<< endl;
				exit(1);
			}
			inputfilestream >> token;
			if (token != "every") {
				global_log->error() << "Expected 'every' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
						<< endl;
				exit(1);
			}
			unsigned instep;
			inputfilestream >> instep;
			inputfilestream >> token;
			if (token != "steps") {
				global_log->error() << "Expected 'steps' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
						<< "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
						<< endl;
				exit(1);
			}
			ChemicalPotential tmu = ChemicalPotential();
			tmu.setMu(icid, imu);
			tmu.setInterval(instep);
			tmu.setInstances(intest);
			if (controlVolume)
				tmu.setControlVolume(x0, y0, z0, x1, y1, z1);
			global_log->info() << setprecision(6) << "chemical Potential " << imu << " component " << icid + 1
					<< " (internally " << icid << ") conduct " << intest << " tests every " << instep << " steps: "
					<< endl;
			global_log->info() << flush;
			_lmu.push_back(tmu);
			global_log->info() << " pushed back." << endl;
		} else if (token == "planckConstant") {
			inputfilestream >> h;
		} else if (token == "Widom") {
			widom = true;
		} else if (token == "cavity") {
			unsigned cavity_cid;
			inputfilestream >> cavity_cid;
			cavity_cid--;
			inputfilestream >> token;
			if (token != "radius") {
				global_log->error() << "Expected 'radius' instead of '" << token << "'.\n";
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
				global_log->error() << "Expected 'grid' instead of '" << token << "'.\n";
				global_log->debug()
						<< "Syntax: cavity <cid> radius <R> <coord_max> grid <N_x> <N_y> <N_z> every <n_step> steps\n";
				exit(1);
			}
			unsigned gridx, gridy, gridz;
			inputfilestream >> gridx >> gridy >> gridz;
			inputfilestream >> token;
			if (token != "every") {
				global_log->error() << "Expected 'every' instead of '" << token << "'.\n";
				global_log->debug()
						<< "Syntax: cavity <cid> radius <R> <coord_max> grid <N_x> <N_y> <N_z> every <n_step> steps\n";
				exit(1);
			}
			unsigned cavity_steps;
			inputfilestream >> cavity_steps;
			inputfilestream >> token;
			if (token != "steps") {
				global_log->error() << "Expected 'steps' instead of '" << token << "'.\n";
				global_log->debug()
						<< "Syntax: cavity <cid> radius <R> <coord_max> grid <N_x> <N_y> <N_z> every <n_step> steps\n";
				exit(1);
			}
			this->_mcav[cavity_cid] = CavityEnsemble();
			this->_mcav[cavity_cid].setInterval(cavity_steps);
			this->_mcav[cavity_cid].setMaxNeighbours(neighbours, cavity_radius * cavity_radius);
			cavity_grid[cavity_cid] = new double[3];
			(cavity_grid[cavity_cid])[0] = gridx;
			(cavity_grid[cavity_cid])[1] = gridy;
			(cavity_grid[cavity_cid])[2] = gridz;
		} else if (token == "NVE") {
			/* TODO: Documentation, what it does (no "Enerstat" at the moment) */
			_domain->thermostatOff();
			global_log->error() << "NVE ensemble not implemented." << endl;
			this->exit(1);
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
			global_log->error() << "tersoff no longer supported." << std::endl;
			global_simulation->exit(-1);
		} else if (token == "WallFun_LJ_9_3") {
			double rho_w, sig_w, eps_w, y_off, y_cut;
			unsigned numComponents;
			_applyWallFun_LJ_9_3 = true;
			inputfilestream >> numComponents >> rho_w >> sig_w >> eps_w >> y_off >> y_cut;
			double *xi_sf = new double[numComponents];
			double *eta_sf = new double[numComponents];
			for (unsigned nc = 0; nc < numComponents; nc++) {
				inputfilestream >> xi_sf[nc] >> eta_sf[nc];
			}

			std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
			_wall = new Wall();
			_wall->initializeLJ93(components, rho_w, sig_w, eps_w, xi_sf, eta_sf, y_off, y_cut);
			delete[] xi_sf;
			delete[] eta_sf;

		} else if (token == "WallFun_LJ_10_4") {
			double rho_w, sig_w, eps_w, y_off, y_cut, Delta;
			unsigned numComponents;
			_applyWallFun_LJ_10_4 = true;
			inputfilestream >> numComponents >> rho_w >> sig_w >> eps_w >> y_off >> y_cut >> Delta;
			double *xi_sf = new double[numComponents];
			double *eta_sf = new double[numComponents];
			for (unsigned nc = 0; nc < numComponents; nc++) {
				inputfilestream >> xi_sf[nc] >> eta_sf[nc];
			}

			std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
			_wall = new Wall();
			_wall->initializeLJ104(components, rho_w, sig_w, eps_w, xi_sf, eta_sf, y_off, y_cut, Delta);
			delete[] xi_sf;
			delete[] eta_sf;

		}

		else if (token == "Mirror") {
			double yMirr, forceConstant;
			_applyMirror = true;
			inputfilestream >> yMirr >> forceConstant;
			std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
			_mirror = new Mirror();
			_mirror->initialize(components, yMirr, forceConstant);
		} else if (token == "slabsLRC") {
			double slabs;
			inputfilestream >> slabs;
			_longRangeCorrection = new Planar(_cutoffRadius, _LJCutoffRadius, _domain, _domainDecomposition,
					_moleculeContainer, slabs, global_simulation);
			_longRangeCorrection->init();
		} else if (token == "NumberOfFluidComponents") {
			double numFluidComp;
			inputfilestream >> numFluidComp;
			_domain->setNumFluidComponents(numFluidComp);
		} else if (token == "flagsNEMD") {
			std::string strFlag;
			inputfilestream >> strFlag;
			if (strFlag == "CHANGE_COMPONENT_AC_TO_N2")
				_flagsNEMD = (_flagsNEMD | CHANGE_COMPONENT_AC_TO_N2);
		}

		/** mheinen 2015-07-27 --> TEMPERATURE_CONTROL
		 *
		 * Temperature Control (Slab Thermostat)
		 *
		 * - applicable to different regions
		 * - using different target temperatures
		 * - using different number of slabs
		 * - componentwise if required
		 * - thermostating only selected directions: x, y, z, xy, xz, yz or xyz
		 *
		 **/
		else if (token == "TemperatureControl" || token == "Temperaturecontrol" || token == "temperatureControl") {

			            string strToken;

			            inputfilestream >> strToken;

			            if(strToken == "param")
			            {

			                unsigned long nControlFreq;
			                unsigned long nStart;
			                unsigned long nStop;

			                inputfilestream >> nControlFreq;
			                inputfilestream >> nStart;
			                inputfilestream >> nStop;

			                if(_temperatureControl == NULL)
			                {
			                    _temperatureControl = new TemperatureControl(_domain, _domainDecomposition, nControlFreq, nStart, nStop);

			                    // turn off explosion heuristics
			                    _domain->SetExplosionHeuristics(false);
			                }
			                else
			                {
			                  global_log->error() << "TemperatureControl object allready exist, programm exit..." << endl;
			                  exit(-1);
			                }
			            }
			            else if(strToken == "heat")
			            {
							paramLineHeat paramLine;
							inputfilestream >> paramLine.nWriteFreq;
							inputfilestream >> paramLine.slabsKey;
							inputfilestream >> paramLine.slabsVal;
							inputfilestream >> paramLine.nWriteFreqRegions;

			                if(_temperatureControl != NULL)
			                {
			                    // set heat supply measurement parameters
			                    _temperatureControl->SetDeltaEkinParameters(paramLine);
			                }
			                else
			                {
			                    global_log->error() << "TemperatureControl object doesnt exist, programm exit..." << endl;
			                    exit(-1);
			                }
			            }
			            else if (strToken == "region")
			            {
			                double dLowerCorner[3];
			                double dUpperCorner[3];
			                string strSubdivisionKey;
			                string strSubdivisionVal;
			                unsigned int nComp;

			                // read lower corner
			                for(unsigned short d=0; d<3; d++)
			                {
			                    inputfilestream >> dLowerCorner[d];
			                }

			                // read upper corner
			                for(unsigned short d=0; d<3; d++)
			                {
			                	string strToken;
			                	inputfilestream >> strToken;
			                	if(strToken == "box")
			                		dUpperCorner[d] = _domain->getGlobalLength(d);
			                	else
			                		dUpperCorner[d] = atof(strToken.c_str() );

			//                	cout << "uc[" << d << "] = " << dUpperCorner[d] << endl;
			                }

			                // target component / temperature
			                inputfilestream >> strSubdivisionKey;
			                inputfilestream >> strSubdivisionVal;
			                inputfilestream >> nComp;

			            	string strToken;
			            	int nTemperatureControlType = TCT_UNKNOWN;
			                double dTargetTemperature[2];
			                double dTemperatureExponent;
			                string strTransDirections;

			                unsigned long nStartAdjust = 0;
			                unsigned long nStopAdjust  = 0;
			                unsigned long nAdjustFreq  = 0;

			            	inputfilestream >> strToken;

			            	if (strToken == "constant")
			            		nTemperatureControlType = TCT_CONSTANT_TEMPERATURE;
			            	else if (strToken == "gradient")
			            		nTemperatureControlType = TCT_TEMPERATURE_GRADIENT;
			            	else if (strToken == "gradient_raise")
			            		nTemperatureControlType = TCT_TEMPERATURE_GRADIENT_RAISE;
			            	else if (strToken == "gradient_lower")
			            		nTemperatureControlType = TCT_TEMPERATURE_GRADIENT_LOWER;
			            	else if(strToken == "adjust")
			            		nTemperatureControlType = TCT_TEMPERATURE_ADJUST;
			            	else
							{
								cout << "TemperatureControl: No valid TemperatureControl type! Program exit...  ";
								exit(-1);
							}

							inputfilestream >> dTargetTemperature[0];

							if( TCT_TEMPERATURE_GRADIENT       == nTemperatureControlType ||
								TCT_TEMPERATURE_GRADIENT_LOWER == nTemperatureControlType ||
								TCT_TEMPERATURE_GRADIENT_RAISE == nTemperatureControlType ||
								TCT_TEMPERATURE_ADJUST         == nTemperatureControlType )
							{
								inputfilestream >> strToken;
								if(strToken != "to")
								{
									cout << "TemperatureControl: Wrong statement in cfg-file, expected 'to'! Program exit...  ";
									exit(-1);
								}
								inputfilestream >> dTargetTemperature[1];
							}

							if( TCT_TEMPERATURE_GRADIENT_LOWER == nTemperatureControlType ||
								TCT_TEMPERATURE_GRADIENT_RAISE == nTemperatureControlType ||
								TCT_TEMPERATURE_ADJUST         == nTemperatureControlType )
							{
								inputfilestream >> strToken;

								char * cstr = new char [strToken.length()+1];
								std::strcpy (cstr, strToken.c_str());
								char* pch;

								pch = strtok(cstr, ":");
								nStartAdjust = atoi(pch);
			//                    cout << "nStartAdjust = " << nStartAdjust << endl;

								pch = strtok (NULL, ":");
								nAdjustFreq = atoi(pch);
			//                    cout << "nAdjustFreq = " << nAdjustFreq << endl;

								pch = strtok (NULL, ":");
								nStopAdjust = atoi(pch);
			//                    cout << "nStopAdjust = " << nStopAdjust << endl;

								delete[] cstr;
							}

							inputfilestream >> dTemperatureExponent;
							inputfilestream >> strTransDirections;

			                if( strTransDirections != "x"  && strTransDirections != "y"  && strTransDirections != "z"  &&
			                    strTransDirections != "xy" && strTransDirections != "xz" && strTransDirections != "yz" &&
			                    strTransDirections != "xyz")
			                {
			                      global_log->error() << "TemperatureControl: Wrong statement! Expected x, y, z, xy, xz, yz or xyz!" << endl;
			                      exit(-1);
			                }

			                if(_temperatureControl == NULL)
			                {
			                      global_log->error() << "TemperatureControl object doesnt exist, programm exit..." << endl;
			                      exit(-1);
			                }
			                else
			                {
			                    // add regions
			                	tec::ControlRegion* region = new tec::ControlRegion( _temperatureControl, dLowerCorner, dUpperCorner, nComp, dTargetTemperature, dTemperatureExponent, strTransDirections,
			                			                                     nTemperatureControlType, nStartAdjust, nStopAdjust, nAdjustFreq );
			                	_temperatureControl->AddRegion(region);

			                    unsigned short refCoords[6];

			                    for(unsigned short d=0; d<6; d++)
			                        inputfilestream >> refCoords[d];

			                    region->PrepareAsObserver(refCoords);

			                    if(_distControl != NULL)
			                    	_distControl->registerObserver(region);

			                    // subdivision of region
			                    if("num" == strSubdivisionKey)
			                    {
			                    	unsigned int nNumSlabs = atoi(strSubdivisionVal.c_str() );
			                    	region->SetSubdivision(nNumSlabs);
			                    }
			                    else if("width" == strSubdivisionKey)
			                    {
			                    	double dSlabWidth = atof(strSubdivisionVal.c_str() );
			                    	region->SetSubdivision(dSlabWidth);
			                    }
			                    else
			                    {
			                    	global_log->error() << "TemperatureControl: wrong subdivision statement. Expected 'num' or 'width' keyword! Programm exit..." << endl;
			                    	exit(-1);
			                	}
			                }

			            }  // if (strToken == "region")
			            else
			            {
			                global_log->error() << "TemperatureControl: Wrong statement in cfg, programm exit..." << endl;
			                exit(-1);
			            }

			        // <-- TEMPERATURE_CONTROL


			        // mheinen 2015-05-27 --> DRIFT_CONTROL
			         } else if (token == "DriftControl" || token == "Driftcontrol" || token == "driftControl") {

			             string strToken;

			             inputfilestream >> strToken;

			             if(strToken == "param")
			             {
			                 unsigned long nControlFreq;
			                 unsigned long nStart;
			                 unsigned long nStop;

			                 inputfilestream >> nControlFreq;
			                 inputfilestream >> nStart;
			                 inputfilestream >> nStop;

			                 if(_driftControl == NULL)
			                 {
			                     _driftControl = new DriftControl(_domain, _domainDecomposition, nControlFreq, nStart, nStop);
			                 }
			                 else
			                 {
			                   global_log->error() << "DriftControl object allready exist, programm exit..." << endl;
			                   exit(-1);
			                 }
			             }
			             else if (strToken == "region")
			             {
			                 double dLowerCorner[3];
			                 double dUpperCorner[3];
			                 unsigned int nTargetComponent;
			                 double dDirection[3];
			                 double dDriftVeloTargetVal;

			                 // read lower corner
			                 for(unsigned short d=0; d<3; d++)
			                 {
			                     inputfilestream >> dLowerCorner[d];
			                 }

			                 // read upper corner
			                 for(unsigned short d=0; d<3; d++)
			                 {
			                 	string strToken;
			                 	inputfilestream >> strToken;
			                 	if(strToken == "box")
			                 		dUpperCorner[d] = _domain->getGlobalLength(d);
			                 	else
			                 		dUpperCorner[d] = atof(strToken.c_str() );

			 //                	cout << "uc[" << d << "] = " << dUpperCorner[d] << endl;
			                 }

			                 // target component ID
			                 inputfilestream >> nTargetComponent;

			                 // drift velocity target direction
			                 for(unsigned short d=0; d<3; d++)
			                 {
			                     inputfilestream >> dDirection[d];
			                 }

			                 // drift velocity target value
			                 inputfilestream >> dDriftVeloTargetVal;

			                 if(_driftControl == NULL)
			                 {
			                       global_log->error() << "DriftControl object doesnt exist, programm exit..." << endl;
			                       exit(-1);
			                 }
			                 else
			                 {
			                     // add regions
			                	 drc::ControlRegion* region = new drc::ControlRegion(_driftControl, dLowerCorner, dUpperCorner, nTargetComponent, dDirection, dDriftVeloTargetVal);
			                     _driftControl->AddRegion(region);

			                     unsigned short refCoords[6];

			                     for(unsigned short d=0; d<6; d++)
			                         inputfilestream >> refCoords[d];

			                     region->PrepareAsObserver(refCoords);

			                     if(_distControl != NULL)
			                     	_distControl->registerObserver(region);
			                 }
			             }
			             else
			             {
			                 global_log->error() << "DriftControl: Wrong statement in cfg, programm exit..." << endl;
			                 exit(-1);
			             }
			         // <-- DRIFT_CONTROL


			         // mheinen 2015-05-29 --> DENSITY_CONTROL
			         } else if (token == "DensityControl" || token == "Densitycontrol" || token == "densityControl") {

			             string strToken;

			             inputfilestream >> strToken;

			             if(strToken == "param")
			             {
			                 unsigned long nControlFreq;
			                 unsigned long nStart;
			                 unsigned long nStop;

			                 inputfilestream >> nControlFreq;
			                 inputfilestream >> nStart;
			                 inputfilestream >> nStop;

			                 if(_densityControl == NULL)
			                 {
			                     _densityControl = new DensityControl(_domainDecomposition, _domain, nControlFreq, nStart, nStop);
			                 }
			                 else
			                 {
			                   global_log->error() << "DensityControl object allready exist, programm exit..." << endl;
			                   exit(-1);
			                 }
			             }
			             else if (strToken == "region")
			             {
			                 double dLowerCorner[3];
			                 double dUpperCorner[3];
			                 unsigned int nTargetComponent;
			                 double dTargetDensity;

			                 // read lower corner
			                 for(unsigned short d=0; d<3; d++)
			                 {
			                     inputfilestream >> dLowerCorner[d];
			                 }

			                 // read upper corner
			                 for(unsigned short d=0; d<3; d++)
			                 {
			                 	string strToken;
			                 	inputfilestream >> strToken;
			                 	if(strToken == "box")
			                 		dUpperCorner[d] = _domain->getGlobalLength(d);
			                 	else
			                 		dUpperCorner[d] = atof(strToken.c_str() );

			 //                	cout << "uc[" << d << "] = " << dUpperCorner[d] << endl;
			                 }

			                 // target component ID
			                 inputfilestream >> nTargetComponent;

			                 // target density
			                 inputfilestream >> dTargetDensity;

			                 if(_densityControl == NULL)
			                 {
			                       global_log->error() << "DensityControl object doesnt exist, programm exit..." << endl;
			                       exit(-1);
			                 }
			                 else
			                 {
			                     // add regions
			                     dec::ControlRegion* region = new dec::ControlRegion(_densityControl, dLowerCorner, dUpperCorner, nTargetComponent, dTargetDensity);
			                     _densityControl->AddRegion(region);

			                     unsigned short refCoords[6];

			                     for(unsigned short d=0; d<6; d++)
			                    	 inputfilestream >> refCoords[d];

			                      region->PrepareAsObserver(refCoords);

			                      if(_distControl != NULL)
			                    	  _distControl->registerObserver(region);
			                 }
			             }
			             else
			             {
			                 global_log->error() << "DensityControl: Wrong statement in cfg, programm exit..." << endl;
			                 exit(-1);
			             }
			         // <-- DENSITY_CONTROL


			         // mheinen 2015-03-16 --> DISTANCE_CONTROL
			         } else if (token == "DistControl" || token == "distControl" ) {

						string strToken;
						inputfilestream >> strToken;

						if(strToken == "param")
						{
							unsigned long nUpdateFreq = 1000;
							unsigned int nWriteFreqProfiles = 1000;
							string strSubdivisionKey;
			                string strSubdivisionVal;

							inputfilestream >> nUpdateFreq >> nWriteFreqProfiles >> strSubdivisionKey >> strSubdivisionVal;

							if(_distControl == NULL)
							{
								_distControl = new DistControl(_domainDecomposition, _domain, nUpdateFreq, nWriteFreqProfiles);

								// subdivision of system to sample profiles
			                    if("num" == strSubdivisionKey)
			                    {
			                    	unsigned int nNumSlabs = atoi(strSubdivisionVal.c_str() );
			                    	_distControl->SetSubdivision(nNumSlabs);
			                    }
			                    else if("width" == strSubdivisionKey)
			                    {
			                    	double dSlabWidth = atof(strSubdivisionVal.c_str() );
			                    	_distControl->SetSubdivision(dSlabWidth);
			                    }
			                    else
			                    {
			                    	global_log->error() << "DistControl: wrong subdivision statement. Expected 'num' or 'width' keyword! Programm exit..." << endl;
			                    	exit(-1);
			                	}
			                    _distControl->PrepareSubdivision();
			                    _distControl->PrepareDataStructures();
							}
							else
							{
								global_log->error() << "DistControl object allready exist, programm exit..." << endl;
								exit(-1);
							}
						}
						else if(strToken == "method")
						{
							if(_distControl == NULL)
							{
								global_log->error() << "DistControl object does not exist, programm exit..." << endl;
								exit(-1);
							}
							std::string strMethod;
							int nMethod;
							unsigned short nVal1, nVal2, nVal3;
							double dVal;

							inputfilestream >> strMethod;
							if(strMethod == "density")
							{
								inputfilestream >> dVal >> nVal3;
								nMethod = DCUM_DENSITY_PROFILE;
							}
							else if(strMethod == "denderiv")
							{
								inputfilestream >> nVal1 >> nVal2 >> nVal3;
								nMethod = DCUM_DENSITY_PROFILE_DERIVATION;
							}
							_distControl->SetUpdateMethod(nMethod, nVal1, nVal2, nVal3, dVal);
						}
						else if(strToken == "init")
						{
							if(_distControl == NULL)
							{
								global_log->error() << "DistControl object does not exist, programm exit..." << endl;
								exit(-1);
							}
							std::string strMethod;
							int nMethod;
							double dVal1 = 0.;
							double dVal2 = 0.;
							std::string strVal = "unknown";
							unsigned long nVal = 0;

							inputfilestream >> strMethod;
							if(strMethod == "startconfig")
							{
								nMethod = DCIM_START_CONFIGURATION;
							}
							else if(strMethod == "values")
							{
								inputfilestream >> dVal1 >> dVal2;
								nMethod = DCIM_MIDPOINT_VALUES;
							}
							else if(strMethod == "file")
							{
								inputfilestream >> strVal >> nVal;
								nMethod = DCIM_READ_FROM_FILE;
							}
							else
								nMethod = DCIM_UNKNOWN;
							_distControl->SetInitMethod(nMethod, dVal1, dVal2, strVal, nVal);
						}
			         // <-- DISTANCE_CONTROL


			         // mheinen 2015-03-18 --> REGION_SAMPLING
			         } else if (token == "RegionSampling" || token == "regionSampling" ) {

			             if(_regionSampling == NULL)
			             {
			                 _regionSampling = new RegionSampling(_domain, _domainDecomposition);
			             }

			             inputfilestream >> token;

			             if(token == "region" || token == "Region")
			             {
			                 double dLowerCorner[3];
			                 double dUpperCorner[3];

			                 // read lower corner
			                 for(unsigned short d=0; d<3; d++)
			                 {
			                     inputfilestream >> dLowerCorner[d];
			                 }

			                 // read upper corner
			                 for(unsigned short d=0; d<3; d++)
			                 {
			                 	string strToken;
			                 	inputfilestream >> strToken;
			                 	if(strToken == "box")
			                 		dUpperCorner[d] = _domain->getGlobalLength(d);
			                 	else
			                 		dUpperCorner[d] = atof(strToken.c_str() );

			 //                	cout << "uc[" << d << "] = " << dUpperCorner[d] << endl;
			                 }

			                 SampleRegion* region = new SampleRegion(_regionSampling, dLowerCorner, dUpperCorner);
			                 _regionSampling->AddRegion(region);

			                 unsigned short refCoords[6];

			                 for(unsigned short d=0; d<6; d++)
			                     inputfilestream >> refCoords[d];

			                 region->PrepareAsObserver(refCoords);

			                 if(_distControl != NULL)
			                 	_distControl->registerObserver(region);
			             }
			             else if(token == "Trho" || token == "tRho")
			             {
							unsigned short nID;
							unsigned long initSamplingProfiles;
							unsigned long writeFrequencyProfiles;
							string strSubdivisionKey;
							string strSubdivisionVal;

							inputfilestream >> nID >> initSamplingProfiles >> writeFrequencyProfiles >> strSubdivisionKey >> strSubdivisionVal;

							SampleRegion* region = _regionSampling->GetSampleRegion(nID);
							region->SetParamProfiles( initSamplingProfiles, writeFrequencyProfiles);

							// subdivision of region
							if("num" == strSubdivisionKey)
							{
								unsigned int nNumSlabs = atoi(strSubdivisionVal.c_str() );
								region->SetSubdivisionProfiles(nNumSlabs);
							}
							else if("width" == strSubdivisionKey)
							{
								double dSlabWidth = atof(strSubdivisionVal.c_str() );
								region->SetSubdivisionProfiles(dSlabWidth);
							}
							else
							{
								global_log->error() << "RegionSampling: wrong subdivision statement. Expected 'num' or 'width' keyword! Programm exit..." << endl;
								exit(-1);
							}

						}
						else if(token == "vdf" || token == "VDF")
						{
							unsigned short nID;
							unsigned long initSamplingVDF;
							unsigned long writeFrequencyVDF;
							string strSubdivisionKey;
							string strSubdivisionVal;
							unsigned int nNumDiscreteStepsVDF;
							double dVeloMax;

							inputfilestream >> nID >> initSamplingVDF >> writeFrequencyVDF;
							inputfilestream >> strSubdivisionKey >> strSubdivisionVal >> nNumDiscreteStepsVDF >> dVeloMax;

							SampleRegion* region = _regionSampling->GetSampleRegion(nID);
							region->SetParamVDF( initSamplingVDF, writeFrequencyVDF, nNumDiscreteStepsVDF, dVeloMax);

							// subdivision of region
							if("num" == strSubdivisionKey)
							{
								unsigned int nNumSlabs = atoi(strSubdivisionVal.c_str() );
								region->SetSubdivisionVDF(nNumSlabs);
							}
							else if("width" == strSubdivisionKey)
							{
								double dSlabWidth = atof(strSubdivisionVal.c_str() );
								region->SetSubdivisionVDF(dSlabWidth);
							}
							else
							{
								global_log->error() << "RegionSampling: wrong subdivision statement. Expected 'num' or 'width' keyword! Programm exit..." << endl;
								exit(-1);
							}
						}
						else
						{
							global_log->error() << "RegionSampling options using wrong statements, programm exit..." << endl;
							exit(-1);
						}

						// <-- REGION_SAMPLING

			// mheinen 2017-03-01 --> PARTICLE_TRACKER
			} else if (token == "Particletracker" || token == "ParticleTracker" ) {
				unsigned long simstepStart, simstepStop;
				inputfilestream >> simstepStart >> simstepStop;
				if(NULL == _particleTracker) {
					_particleTracker = new ParticleTracker(_domainDecomposition, simstepStart, simstepStop);
				}
			// <-- PARTICLE_TRACKER

		} else if (token == "finalCheckpoint") {
			std::string strKey;
			inputfilestream >> strKey;
			if (strKey == "binary")
			{
				_finalCheckpoint = true;
				_finalCheckpointBinary = true;
			}
		} else {
			if (token != "")
				global_log->warning() << "Did not process unknown token " << token << endl;
		}
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

	// read particle data
	_maxMoleculeId = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain, _domainDecomposition);

	if (this->_LJCutoffRadius == 0.0)
		_LJCutoffRadius = this->_cutoffRadius;
	//_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

	// @todo comment
#ifndef MARDYN_WR
	_integrator = new Leapfrog(timestepLength);
#else
	_integrator = new ExplicitEuler(timestepLength);
#endif

	// test new Decomposition
	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	unsigned idi = _lmu.size();
	unsigned j = 0;
	std::list<ChemicalPotential>::iterator cpit;
	for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
		if (widom)
			cpit->enableWidom();
		cpit->setIncrement(idi);
		double tmp_molecularMass = global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->m();
		cpit->setSystem(_domain->getGlobalLength(0), _domain->getGlobalLength(1), _domain->getGlobalLength(2),
				tmp_molecularMass);
		cpit->setGlobalN(global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->getNumMolecules());
		cpit->setNextID(j + (int) (1.001 * (256 + _maxMoleculeId)));

		cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0), _moleculeContainer->getBoundingBoxMax(0),
				_moleculeContainer->getBoundingBoxMin(1), _moleculeContainer->getBoundingBoxMax(1),
				_moleculeContainer->getBoundingBoxMin(2), _moleculeContainer->getBoundingBoxMax(2));
		/* TODO: thermostat */
		double Tcur = _domain->getCurrentTemperature(0);
		/* FIXME: target temperature from thermostat ID 0 or 1?  */
		double Ttar =
				_domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;
		cpit->submitTemperature(Tcur);
		if (h != 0.0)
			cpit->setPlanckConstant(h);

		j++;
	}

	unsigned long Nbasis = 3 * (this->_domain->getglobalNumMolecules() + 3 * this->_lmu.size());
	map<unsigned, CavityEnsemble>::iterator ceit;
	for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {
		ceit->second.setSystem(_domain->getGlobalLength(0), _domain->getGlobalLength(1), _domain->getGlobalLength(2));
		ceit->second.setIdOffset((1 + ceit->first) * Nbasis);

		ceit->second.setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0),
				_moleculeContainer->getBoundingBoxMax(0), _moleculeContainer->getBoundingBoxMin(1),
				_moleculeContainer->getBoundingBoxMax(1), _moleculeContainer->getBoundingBoxMin(2),
				_moleculeContainer->getBoundingBoxMax(2));
		double Tcur = _domain->getCurrentTemperature(0);
		double Ttar =
				_domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(0);
		if ((Tcur < 0.667 * Ttar) || (Tcur > 1.5 * Ttar))
			Tcur = Ttar;
		ceit->second.submitTemperature(Tcur);

		ceit->second.init(global_simulation->getEnsemble()->getComponent(ceit->first), (cavity_grid[ceit->first])[0],
				(cavity_grid[ceit->first])[1], (cavity_grid[ceit->first])[2]);
		ceit->second.communicateNumCavities(this->_domainDecomposition);
	}
}
