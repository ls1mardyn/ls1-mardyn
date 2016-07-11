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

#include "longRange/LongRangeCorrection.h"
#include "longRange/Homogeneous.h"
#include "longRange/Planar.h"
#include "particleContainer/adapter/VectorizationTuner.h"

using namespace std;


void Simulation::initConfigOldstyle(const string& inputfilename) {
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif

        global_log->info() << "\n# Please address your questions and suggestions to the ls1 mardyn contact point:\n# \n# E-mail: contact@ls1-mardyn.de\n# \n# Phone: +49 631 205 3227\n# Fax: +49 631 205 3835\n# University of Kaiserslautern\n# Computational Molecular Engineering\n# Erwin-Schroedinger-Str. 44\n# D-67663 Kaiserslautern, Germany\n# \n# http://www.ls1-mardyn.de/\n\n";

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
				bool hetero=false, cutsmaller=false, forceRatio=false;
				if(line.find("hetero") != string::npos){
					hetero=true;
				}
				if(line.find("cutSmaller") != string::npos){
					cutsmaller=true;  // allow domain to be split not only along biggest side
				}
				if(line.find("forceRatio") != string::npos){
					forceRatio=true;  // allow domain to be split not only along biggest side
				}
				_domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, updateFrequency,
						fullSearchThreshold, hetero, cutsmaller, forceRatio);
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
				global_log->error() << "AdaptiveSubCells no longer supported." << std::endl;
				global_simulation->exit(-1);
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
						outputPathAndPrefix, true));
				global_log->debug() << "XyzWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "CavityWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new CavityWriter(writeFrequency,
						outputPathAndPrefix, true));
				global_log->debug() << "CavityWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "CheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new CheckpointWriter(writeFrequency,
						outputPathAndPrefix, true));
				if(writeFrequency == 0 ){
					global_log->error() << "Write frequency must be a positive nonzero integer, but is " << writeFrequency << endl;
					exit(-1);
				}
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
						outputPathAndPrefix));
				global_log->debug() << "VISWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			} else if (token == "MmspdBinWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new MmspdBinWriter(writeFrequency,
						outputPathAndPrefix));
				global_log->debug() << "MmspdBinWriter " << writeFrequency << " '"
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
				Log::global_log->warning() << std::endl << "VTK-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			} else if (token == "VTKGridWriter") {
#ifdef VTK
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;

				if (dynamic_cast<LinkedCells*>(_moleculeContainer)) {
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
						outputPathAndPrefix));
				global_log->debug() << "MPICheckpointWriter " << writeFrequency
						<< " '" << outputPathAndPrefix << "'.\n";
			}
			else if (token == "GammaWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new GammaWriter(writeFrequency,
						outputPathAndPrefix));
				global_log->debug() << "GammaWriter " << writeFrequency << " '"
						<< outputPathAndPrefix << "'.\n";
			}
			else if (token == "VectorizationTuner") {
				string outputPathAndPrefix;
				unsigned int minMoleculeCnt, maxMoleculeCnt;
				int type;
				MoleculeCntIncreaseTypeEnum moleculeCntIncreaseType;
				inputfilestream >> outputPathAndPrefix >> minMoleculeCnt >> maxMoleculeCnt >> type;
				moleculeCntIncreaseType = static_cast<MoleculeCntIncreaseTypeEnum>(type);
				_outputPlugins.push_back(
						new VectorizationTuner(outputPathAndPrefix, minMoleculeCnt, maxMoleculeCnt,
								moleculeCntIncreaseType, _cutoffRadius,_LJCutoffRadius,&_cellProcessor)
				);
				global_log->debug() << "VectorizationTuner " << outputPathAndPrefix << "'.\n";
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
		}
		else if (token == "profileVirial") {
                        _doRecordVirialProfile = true;
                } else if (token == "profileRecordingTimesteps") { /* TODO: suboption of profile */
			inputfilestream >> _profileRecordingTimesteps;
		} else if (token == "profileOutputTimesteps") { /* TODO: suboption of profile */
			inputfilestream >> _profileOutputTimesteps;
		} else if (token == "RDF") {
			double interval;
			unsigned bins;
			inputfilestream >> interval >> bins;
			if (global_simulation->getEnsemble()->components()->size() <= 0) {
				global_log->error()
						<< "PhaseSpaceFile-Specification has to occur before RDF-Token!"
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
		} else if (token == "profiledComponent") { /* TODO: suboption of profile, check if required to enable output in general */
			unsigned cid;
			inputfilestream >> cid;
			cid--;
			_domain->considerComponentInProfile(cid);
		}
		// by Stefan Becker <stefan.becker@mv.uni-kl.de>, token determining the coordinate system of the density profile
		else if (token == "SessileDrop"){
			this->_domain->sesDrop();
		}
		else if (token == "profileOutputPrefix") { /* TODO: suboption of profile */
			inputfilestream >> _profileOutputPrefix;
		} else if (token == "collectThermostatDirectedVelocity") { /* suboption of the thermostat replace with direct thermostat */
			inputfilestream >> _collectThermostatDirectedVelocity;
		} 
			
		else if (token == "thermostat"){
			inputfilestream >> _thermostatType;
			if(_thermostatType == ANDERSEN_THERMOSTAT){
			  inputfilestream >>_nuAndersen;
			  global_log->info() << "Using the Andersen Thermostat with nu = " << _nuAndersen << "\n";
			}
		}
		else if (token == "zOscillator") {
			global_log->error()
					<< "zOscillator was used for the Tersoff potential, which is no longer supported."
					<< std::endl;
			global_simulation->exit(-1);
		}
		// by Stefan Becker
		else if(token == "AlignCentre"){	
			_doAlignCentre = true;
			inputfilestream >> _alignmentInterval >> _alignmentCorrection;
		}
		else if(token == "ComponentForYShift"){
		    _componentSpecificAlignment = true;
		    unsigned cidMin, cidMax;
		    inputfilestream  >> cidMin >> cidMax;
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
							<< "conduct <ntest> tests every <nstep> steps" << endl;
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
						<< "conduct <ntest> tests every <nstep> steps" << endl;
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
					<< instep << " steps: " << endl;
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
		} 
		else if (token == "WallFun_LJ_9_3"){
		  double rho_w, sig_w, eps_w, y_off, y_cut;
		  unsigned numComponents;
		  _applyWallFun = true;
		  inputfilestream  >> numComponents >> rho_w >> sig_w >> eps_w >> y_off >> y_cut;
		  double *xi_sf = new double[numComponents];
		  double *eta_sf = new double[numComponents];
		  for(unsigned nc = 0; nc < numComponents; nc++ ){
		    inputfilestream >> xi_sf[nc] >> eta_sf[nc];
		  }
		 
		  std::vector<Component>* components = global_simulation->getEnsemble()->components();
		  _wall=new Wall();
		  _wall->initialize(components,  rho_w, sig_w, eps_w, xi_sf, eta_sf, y_off, y_cut);
		  delete[] xi_sf;
		  delete[] eta_sf;

		}		
		
		else if (token == "Mirror"){
			double yMirr, forceConstant;
			_applyMirror=true;
			inputfilestream >> yMirr >> forceConstant;
			std::vector<Component>* components = global_simulation->getEnsemble()->components();
			_mirror=new Mirror();
			_mirror->initialize(components, yMirr, forceConstant);
		} else if (token == "slabsLRC") {
			double slabs;
			inputfilestream >> slabs;
			_longRangeCorrection = new Planar(_cutoffRadius,_LJCutoffRadius,_domain,_domainDecomposition,_moleculeContainer,slabs,global_simulation);
		} 
		else if (token == "NumberOfFluidComponents"){
		    double numFluidComp;
		    inputfilestream >> numFluidComp;
		    _domain->setNumFluidComponents(numFluidComp);

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

        } else if (token == "TemperatureControl" || token == "Temperaturecontrol" || token == "temperatureControl") {

            string strToken;

            inputfilestream >> strToken;

            if(strToken == "param")
            {
                unsigned long nControlFreq;
                unsigned long nStart;
                unsigned long nStop;
                bool bUseExplosionHeuristics;

                inputfilestream >> nControlFreq;
                inputfilestream >> nStart;
                inputfilestream >> nStop;
                inputfilestream >> bUseExplosionHeuristics;

                if(_temperatureControl == NULL)
                {
                    _temperatureControl = new TemperatureControl(nControlFreq, nStart, nStop);

                    // turn off explosion heuristics
                    _domain->SetExplosionHeuristics(bUseExplosionHeuristics);
                }
                else
                {
                  global_log->error() << "TemperatureControl object already exist, program exit..." << endl;
                  exit(-1);
                }
            }
            else if (strToken == "region")
            {
                double dLowerCorner[3];
                double dUpperCorner[3];
                unsigned int nNumSlabs;
                unsigned int nComp;
                double dTargetTemperature;
                double dTemperatureExponent;
                string strTransDirections;

                // read lower corner
                for(unsigned short d=0; d<3; d++)
                {
                    inputfilestream >> dLowerCorner[d];
                }

                // read upper corner
                for(unsigned short d=0; d<3; d++)
                {
                    inputfilestream >> dUpperCorner[d];
                }

                // target component / temperature
                inputfilestream >> nNumSlabs;
                inputfilestream >> nComp;
                inputfilestream >> dTargetTemperature;
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
                    _temperatureControl->AddRegion(dLowerCorner, dUpperCorner, nNumSlabs, nComp, dTargetTemperature, dTemperatureExponent, strTransDirections);
                }
            }
            else
            {
                global_log->error() << "TemperatureControl: Wrong statement in cfg, programm exit..." << endl;
                exit(-1);
            }

        // <-- TEMPERATURE_CONTROL
		} else {
			if (token != "")
				global_log->warning() << "Did not process unknown token " << token << endl;
		}
	}
	/*if(inputfilestream.fail()){
		global_log->error() << "Reading input data ended with failure (failbit or badbit). rdstate:"
				<< (inputfilestream.bad() ? "badbit" : "failbit") << endl << "last token was: " << token << endl;
		exit(-1);
	}*/

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
	_maxMoleculeId = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain,
			_domainDecomposition);

	if (this->_LJCutoffRadius == 0.0)
		_LJCutoffRadius = this->_cutoffRadius;
	//_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

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
		cpit->setNextID(j + (int) (1.001 * (256 + _maxMoleculeId)));

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
