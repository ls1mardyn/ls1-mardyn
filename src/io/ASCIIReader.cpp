#include "io/ASCIIReader.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <climits>
#include <string>

#include "Domain.h"
#include "ensemble/BoxDomain.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include "molecules/Molecule.h"
#include "molecules/mixingrules/MixingRuleBase.h"
#include "molecules/mixingrules/LorentzBerthelot.h"

#ifdef ENABLE_MPI
#include "parallel/ParticleData.h"
#include "parallel/DomainDecompBase.h"
#endif

#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"


ASCIIReader::ASCIIReader() {}

void ASCIIReader::setPhaseSpaceFile(std::string filename) {
	_phaseSpaceFile = filename;
}

void ASCIIReader::setPhaseSpaceHeaderFile(std::string filename) {
	_phaseSpaceHeaderFile = filename;
}

void ASCIIReader::readXML(XMLfileUnits& xmlconfig) {
	std::string pspfile;
	if(xmlconfig.getNodeValue(".", pspfile)) {
		pspfile = string_utils::trim(pspfile);
		// only prefix xml dir if path is not absolute
		if (pspfile[0] != '/') {
			pspfile.insert(0, xmlconfig.getDir());
		}
		Log::global_log->info() << "phasespacepoint description file:\t" << pspfile << std::endl;
	}
	setPhaseSpaceFile(pspfile);
}

void ASCIIReader::readPhaseSpaceHeader(Domain* domain, double timestep) {
	std::string token;

	Log::global_log->info() << "Opening phase space header file " << _phaseSpaceHeaderFile << std::endl;
	_phaseSpaceHeaderFileStream.open(_phaseSpaceHeaderFile.c_str());
	_phaseSpaceHeaderFileStream >> token;
	if(token != "mardyn") {
		Log::global_log->error() << _phaseSpaceHeaderFile << " not a valid mardyn input file." << std::endl;
		mardyn_exit(1);
	}

	std::string inputversion;
	_phaseSpaceHeaderFileStream >> token >> inputversion;
	// FIXME: remove tag trunk from file specification?
	if(token != "trunk") {
		Log::global_log->error() << "Wrong input file specifier (\'" << token << "\' instead of \'trunk\')." << std::endl;
		mardyn_exit(1);
	}

	if(std::stoi(inputversion) < 20080701) {
		Log::global_log->error() << "Input version too old (" << inputversion << ")" << std::endl;
		mardyn_exit(1);
	}

	Log::global_log->info() << "Reading phase space header from file " << _phaseSpaceHeaderFile << std::endl;

	std::vector<Component>& dcomponents = *(_simulation.getEnsemble()->getComponents());
	bool header = true; // When the last header element is reached, "header" is set to false

	while(header) {
		char c;
		_phaseSpaceHeaderFileStream >> c;
		if(c == '#') {
			// comment line
			_phaseSpaceHeaderFileStream.ignore(INT_MAX, '\n');
			continue;
		}
		_phaseSpaceHeaderFileStream.putback(c);

		token.clear();
		_phaseSpaceHeaderFileStream >> token;
		Log::global_log->info() << "{{" << token << "}}" << std::endl;

		if((token == "currentTime") || (token == "t")) {
			// set current simulation time
			_phaseSpaceHeaderFileStream >> token;
			_simulation.setSimulationTime(strtod(token.c_str(), NULL));
		} else if((token == "Temperature") || (token == "T")) {
			// set global thermostat temperature
			domain->disableComponentwiseThermostat(); //disable component wise thermostats
			double targetT;
			_phaseSpaceHeaderFileStream >> targetT;
			domain->setGlobalTemperature(targetT);
		} else if((token == "ThermostatTemperature") || (token == "ThT") || (token == "h")) {
			// set up a new thermostat
			int thermostat_id;
			double targetT;
			_phaseSpaceHeaderFileStream >> thermostat_id;
			_phaseSpaceHeaderFileStream >> targetT;
			Log::global_log->info() << "Thermostat number " << thermostat_id << " has T = " << targetT << ".\n";
			domain->setTargetTemperature(thermostat_id, targetT);
		} else if((token == "ComponentThermostat") || (token == "CT") || (token == "o")) {
			// specify a thermostat for a component
			if(!domain->severalThermostats())
				domain->enableComponentwiseThermostat();
			int component_id;
			int thermostat_id;
			_phaseSpaceHeaderFileStream >> component_id >> thermostat_id;
			Log::global_log->info() << "Component " << component_id << " (internally: " << component_id - 1
							   << ") is regulated by thermostat number " << thermostat_id << ".\n";
			component_id--; // FIXME thermostat IDs start with 0 in the program but not in the config file?!
			if(thermostat_id < 0) // thermostat IDs start with 0
				continue;
			domain->setComponentThermostat(component_id, thermostat_id);
		} else if((token == "Undirected") || (token == "U")) {
			// set undirected thermostat
			int thermostat_id;
			_phaseSpaceHeaderFileStream >> thermostat_id;
			domain->enableUndirectedThermostat(thermostat_id);
		} else if((token == "Length") || (token == "L")) {
			// simulation box dimensions
			double globalLength[3];
			_phaseSpaceHeaderFileStream >> globalLength[0] >> globalLength[1] >> globalLength[2];
			delete _simulation.getEnsemble()->domain();
			_simulation.getEnsemble()->domain() = new BoxDomain();
			for(int d = 0; d < 3; d++) {
				static_cast<BoxDomain*>(_simulation.getEnsemble()->domain())->setLength(d, globalLength[d]);
				domain->setGlobalLength(d, _simulation.getEnsemble()->domain()->length(d));
			}
		} else if((token == "HeatCapacity") || (token == "cv") || (token == "I")) {
			unsigned N;
			double U, UU;
			_phaseSpaceFileStream >> N >> U >> UU;
			domain->init_cv(N, U, UU);
		} else if((token == "NumberOfComponents") || (token == "C")) {
			// read in component definitions and
			// read in mixing coefficients

			// components:
			unsigned int numcomponents = 0;
			_phaseSpaceHeaderFileStream >> numcomponents;
			Log::global_log->info() << "Reading " << numcomponents << " components" << std::endl;
			dcomponents.resize(numcomponents);
			for(unsigned int i = 0; i < numcomponents; i++) {
				Log::global_log->info() << "comp. i = " << i << ": " << std::endl;
				dcomponents[i].setID(i);
				unsigned int numljcenters = 0;
				unsigned int numcharges = 0;
				unsigned int numdipoles = 0;
				unsigned int numquadrupoles = 0;
				unsigned int numtersoff = 0; //previously tersoff
				_phaseSpaceHeaderFileStream >> numljcenters >> numcharges >> numdipoles
											>> numquadrupoles >> numtersoff;
				if(numtersoff != 0) {
					Log::global_log->error() << "tersoff no longer supported."
										<< std::endl;
					mardyn_exit(-1);
				}
				double x, y, z, m;
				for(unsigned int j = 0; j < numljcenters; j++) {
					double eps, sigma, tcutoff, do_shift;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> m >> eps >> sigma >> tcutoff >> do_shift;
					dcomponents[i].addLJcenter(x, y, z, m, eps, sigma, tcutoff, (do_shift != 0));
					Log::global_log->info() << "LJ at [" << x << " " << y << " " << z << "], mass: " << m << ", epsilon: "
									   << eps << ", sigma: " << sigma << std::endl;
				}
				for(unsigned int j = 0; j < numcharges; j++) {
					double q;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> m >> q;
					dcomponents[i].addCharge(x, y, z, m, q);
					Log::global_log->info() << "charge at [" << x << " " << y << " " << z << "], mass: " << m << ", q: " << q
									   << std::endl;
				}
				for(unsigned int j = 0; j < numdipoles; j++) {
					double eMyx, eMyy, eMyz, absMy;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> eMyx >> eMyy >> eMyz >> absMy;
					dcomponents[i].addDipole(x, y, z, eMyx, eMyy, eMyz, absMy);
					Log::global_log->info() << "dipole at [" << x << " " << y << " " << z << "] " << std::endl;
				}
				for(unsigned int j = 0; j < numquadrupoles; j++) {
					double eQx, eQy, eQz, absQ;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> eQx >> eQy >> eQz >> absQ;
					dcomponents[i].addQuadrupole(x, y, z, eQx, eQy, eQz, absQ);
					Log::global_log->info() << "quad at [" << x << " " << y << " " << z << "] " << std::endl;
				}
				double IDummy1, IDummy2, IDummy3;
				// FIXME! Was soll das hier? Was ist mit der Initialisierung im Fall I <= 0.
				_phaseSpaceHeaderFileStream >> IDummy1 >> IDummy2 >> IDummy3;
				if(IDummy1 > 0.) dcomponents[i].setI11(IDummy1);
				if(IDummy2 > 0.) dcomponents[i].setI22(IDummy2);
				if(IDummy3 > 0.) dcomponents[i].setI33(IDummy3);
				domain->setProfiledComponentMass(dcomponents[i].m());
				Log::global_log->info() << std::endl;
			}

#ifndef NDEBUG
			for(unsigned int i = 0; i < numcomponents; i++) {
				Log::global_log->debug() << "Component " << (i + 1) << " of " << numcomponents << std::endl;
				Log::global_log->debug() << dcomponents[i] << std::endl;
			}
#endif

			// Mixing coefficients
			for(unsigned int cidi = 0; cidi < numcomponents-1; cidi++) {
				for(unsigned int cidj = cidi + 1; cidj <= numcomponents-1; cidj++) {
					double xi, eta;
					_phaseSpaceHeaderFileStream >> xi >> eta;
#ifndef NDEBUG
					Log::global_log->debug() << "Mixing: " << cidi+1 << " + " << cidj+1
											 << " : xi=" << xi << " eta=" << eta << std::endl;
#endif
					// Only LB mixing rule is supported for now
					LorentzBerthelotMixingRule* mixingrule = new LorentzBerthelotMixingRule();
					mixingrule->setCid1(cidi);
					mixingrule->setCid2(cidj);
					mixingrule->setEta(eta);
					mixingrule->setXi(xi);
					_simulation.getEnsemble()->setMixingrule(mixingrule);
				}
			}

			// read in global factor \epsilon_{RF}
			// FIXME: Maybe this should go better to a seperate token?!
			_phaseSpaceHeaderFileStream >> token;
			domain->setepsilonRF(strtod(token.c_str(), NULL));
			long int fpos;
			if(_phaseSpaceFile == _phaseSpaceHeaderFile) {
				// in the case of a single phase space header + phase space file
				// find out the actual position, because the phase space definition will follow
				// FIXME: is there a more elegant way?
				fpos = _phaseSpaceHeaderFileStream.tellg();
				_phaseSpaceFileStream.seekg(fpos, std::ios::beg);
			}
			// FIXME: Is there a better solution than skipping the rest of the file?
			header = false;
		} else if((token == "NumberOfMolecules") || (token == "N")) {
			// set number of Molecules
			// FIXME: Is this part called in any case as the token is handled in the readPhaseSpace method?
			_phaseSpaceHeaderFileStream >> token;
			domain->setglobalNumMolecules( strtoul(token.c_str(),NULL,0) );
		}
		// LOCATION OF OLD PRESSURE GRADIENT TOKENS
		else {
			Log::global_log->error() << "Invalid token \'" << token << "\' found. Skipping rest of the header." << std::endl;
			header = false;
		}
	}

	_simulation.getEnsemble()->setComponentLookUpIDs();

	_phaseSpaceHeaderFileStream.close();
}

unsigned long
ASCIIReader::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {

	global_simulation->timers()->start("INPUT_OLDSTYLE_INPUT");

#ifdef ENABLE_MPI
	if (domainDecomp->getRank() == 0)
	{ // Rank 0 only
#endif
	Log::global_log->info() << "Opening phase space file " << _phaseSpaceFile << std::endl;
	_phaseSpaceFileStream.open(_phaseSpaceFile.c_str());
	if(!_phaseSpaceFileStream.is_open()) {
		Log::global_log->error() << "Could not open phaseSpaceFile " << _phaseSpaceFile << std::endl;
		mardyn_exit(1);
	}
	Log::global_log->info() << "Reading phase space file " << _phaseSpaceFile << std::endl;
#ifdef ENABLE_MPI
	} // Rank 0 only
#endif

	std::string token;
	std::vector<Component>& dcomponents = *(_simulation.getEnsemble()->getComponents());
	unsigned int numcomponents = dcomponents.size();
	unsigned long nummolecules = 0;
	unsigned long maxid = 0; // stores the highest molecule ID found in the phase space file
	std::string ntypestring("ICRVQD");
	enum class Ndatatype {
		ICRVQDV, ICRVQD, IRV, ICRV
	} ntype = Ndatatype::ICRVQD;

#ifdef ENABLE_MPI
	if (domainDecomp->getRank() == 0)
	{ // Rank 0 only
#endif
	while(_phaseSpaceFileStream && (token != "NumberOfMolecules") && (token != "N")) {
		_phaseSpaceFileStream >> token;
	}
	if((token != "NumberOfMolecules") && (token != "N")) {
		Log::global_log->error() << "Expected the token 'NumberOfMolecules (N)' instead of '" << token << "'" << std::endl;
		mardyn_exit(1);
	}
	_phaseSpaceFileStream >> nummolecules;
#ifdef ENABLE_MPI
	} // Rank 0 only
	// TODO: Better do the following in setGlobalNumMolecules?!
	MPI_Bcast(&nummolecules, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
#endif
	Log::global_log->info() << " number of molecules: " << nummolecules << std::endl;


#ifdef ENABLE_MPI
	if (domainDecomp->getRank() == 0)
	{ // Rank 0 only
#endif
	std::streampos spos = _phaseSpaceFileStream.tellg();
	_phaseSpaceFileStream >> token;
	if((token == "MoleculeFormat") || (token == "M")) {
		_phaseSpaceFileStream >> ntypestring;
		ntypestring.erase(ntypestring.find_last_not_of(" \t\n") + 1);
		ntypestring.erase(0, ntypestring.find_first_not_of(" \t\n"));
		if(ntypestring == "ICRVQDV") ntype = Ndatatype::ICRVQDV;
		else if(ntypestring == "ICRVQD") ntype = Ndatatype::ICRVQD;
		else if(ntypestring == "ICRV") ntype = Ndatatype::ICRV;
		else if(ntypestring == "IRV") ntype = Ndatatype::IRV;
		else {
			Log::global_log->error() << "Unknown molecule format '" << ntypestring << "'" << std::endl;
			mardyn_exit(1);
		}
	} else {
		_phaseSpaceFileStream.seekg(spos);
	}
	Log::global_log->info() << " molecule format: " << ntypestring << std::endl;
	if(numcomponents < 1) {
		Log::global_log->warning() << "No components defined! Setting up single one-centered LJ" << std::endl;
		numcomponents = 1;
		dcomponents.resize(numcomponents);
		dcomponents[0].setID(0);
		dcomponents[0].addLJcenter(0., 0., 0., 1., 1., 1., 6., false);
	}

#ifdef ENABLE_MPI
	} // Rank 0 only
#endif

#ifdef ENABLE_MPI
#define PARTICLE_BUFFER_SIZE  (16*1024)
	ParticleData particle_buff[PARTICLE_BUFFER_SIZE];
	int particle_buff_pos = 0;
	MPI_Datatype mpi_Particle;
	ParticleData::getMPIType(mpi_Particle);

	int size;
	MPI_CHECK(MPI_Type_size(mpi_Particle, &size));
	Log::global_log->debug() << "size of custom datatype is " << size << std::endl;

#endif

	double x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz, Vix, Viy, Viz;
	unsigned long id = 0ul;
	unsigned int componentid = 1;  // Default componentID when using IRV format

	x = y = z = vx = vy = vz = q1 = q2 = q3 = Dx = Dy = Dz = Vix = Viy = Viz = 0.;
	q0 = 1.;

	for(unsigned long i = 0; i < nummolecules; i++) {

#ifdef ENABLE_MPI
		if (domainDecomp->getRank() == 0) { // Rank 0 only
#endif
		switch (ntype) {
			case Ndatatype::ICRVQDV:
				_phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
									  >> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz >> Vix >> Viy >> Viz;
				break;
			case Ndatatype::ICRVQD:
				_phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
									  >> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz;
				break;
			case Ndatatype::ICRV :
				_phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz;
				break;
			case Ndatatype::IRV :
				_phaseSpaceFileStream >> id >> x >> y >> z >> vx >> vy >> vz;
				componentid = 1;
				break;
			default:
				Log::global_log->error() << "[ASCIIReader.cpp] Unknown ntype" << std::endl;
				mardyn_exit(1);
		}
		if((x < 0.0 || x >= domain->getGlobalLength(0))
		   || (y < 0.0 || y >= domain->getGlobalLength(1))
		   || (z < 0.0 || z >= domain->getGlobalLength(2))) {
			Log::global_log->warning() << "Molecule " << id << " out of box: " << x << ";" << y << ";" << z << std::endl;
		}

		if(componentid > numcomponents) {
			Log::global_log->error() << "Molecule id " << id
								<< " has a component ID greater than the existing number of components: "
								<< componentid
								<< ">"
								<< numcomponents << std::endl;
			mardyn_exit(1);
		}
		// ComponentIDs are used as array IDs, hence need to start at 0.
		// In the input files they always start with 1 so we need to adapt that all the time.
		componentid--;
		// store only those molecules within the domain of this process
		// The necessary check is performed in the particleContainer addParticle method
		// FIXME: Datastructures? Pass pointer instead of object, so that we do not need to copy?!
		Molecule m1 = Molecule(id, &dcomponents[componentid], x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz);
#ifdef ENABLE_MPI
		ParticleData::MoleculeToParticleData(particle_buff[particle_buff_pos], m1);
	} // Rank 0 only

	particle_buff_pos++;
	if ((particle_buff_pos >= PARTICLE_BUFFER_SIZE) || (i == nummolecules - 1)) {
		//MPI_Bcast(&particle_buff_pos, 1, MPI_INT, 0, MPI_COMM_WORLD);
		Log::global_log->debug() << "broadcasting(sending/receiving) particles with buffer_position " << particle_buff_pos << std::endl;
		MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, 0, MPI_COMM_WORLD); // TODO: MPI_COMM_WORLD
		for (int j = 0; j < particle_buff_pos; j++) {
			Molecule m;
			ParticleData::ParticleDataToMolecule(particle_buff[j], m);
			// only add particle if it is inside of the own domain!
			if(particleContainer->isInBoundingBox(m.r_arr().data())) {
				particleContainer->addParticle(m, true, false);
			}
			componentid=m.componentid();

			dcomponents[componentid].incNumMolecules();
			domain->setglobalRotDOF(dcomponents[componentid].getRotationalDegreesOfFreedom() + domain->getglobalRotDOF());

			if(m.getID() > maxid) maxid = m.getID();

			// Only called inside GrandCanonical
			global_simulation->getEnsemble()->storeSample(&m, componentid);
		}
		Log::global_log->debug() << "broadcasting(sending/receiving) complete" << particle_buff_pos << std::endl;
		particle_buff_pos = 0;
	}
#else
		if(particleContainer->isInBoundingBox(m1.r_arr().data())) {
			particleContainer->addParticle(m1, true, false);
		}

		componentid = m1.componentid();
		// TODO: The following should be done by the addPartice method.
		dcomponents[componentid].incNumMolecules();
		domain->setglobalRotDOF(dcomponents[componentid].getRotationalDegreesOfFreedom() + domain->getglobalRotDOF());

		if(id > maxid) maxid = id;

		// Only called inside GrandCanonical
		global_simulation->getEnsemble()->storeSample(&m1, componentid);

#endif
	}
	Log::global_log->info() << "Reading Molecules done" << std::endl;


#ifdef ENABLE_MPI
	if (domainDecomp->getRank() == 0)
	{ // Rank 0 only
#endif
	_phaseSpaceFileStream.close();
#ifdef ENABLE_MPI
	} // Rank 0 only
#endif

	global_simulation->timers()->stop("INPUT_OLDSTYLE_INPUT");
	global_simulation->timers()->setOutputString("INPUT_OLDSTYLE_INPUT", "Initial IO took:                 ");
	global_simulation->timers()->print("INPUT_OLDSTYLE_INPUT");
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Type_free(&mpi_Particle) );
#endif
	domain->setglobalNumMolecules(nummolecules);
	domain->setglobalRho(nummolecules / domain->getGlobalVolume());
	return maxid;
}
