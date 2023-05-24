/*
 * MPI_IOReader.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: andal
 */

#include "MPI_IOReader.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <string>
#include <sstream>
#include <climits>

#include "Domain.h"
#include "ensemble/BoxDomain.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include "molecules/Molecule.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecompBase.h"
#include "parallel/ParticleData.h"
#endif

#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/Timer.h"

//#include <time.h>


MPI_IOReader::MPI_IOReader() {
	// TODO Auto-generated constructor stub

}

MPI_IOReader::~MPI_IOReader() {
	// TODO Auto-generated destructor stub
}

void MPI_IOReader::setPhaseSpaceFile(std::string filename) {
	_phaseSpaceFile = filename;
}

void MPI_IOReader::setPhaseSpaceHeaderFile(std::string filename) {
	_phaseSpaceHeaderFile = filename;
}

void MPI_IOReader::readPhaseSpaceHeader(Domain* domain, double timestep) {
	std::string token, token2;

	global_log->info() << "Opening phase space header file " << _phaseSpaceHeaderFile << std::endl;
	_phaseSpaceHeaderFileStream.open(_phaseSpaceHeaderFile.c_str());
	_phaseSpaceHeaderFileStream >> token;
	if(token != "mardyn") {
		global_log->error() << _phaseSpaceHeaderFile << " not a valid mardyn input file." << std::endl;
		Simulation::exit(1);
	}

	std::string inputversion;
	_phaseSpaceHeaderFileStream >> token >> inputversion;
	// FIXME: remove tag trunk from file specification?
	if(token != "trunk") {
		global_log->error() << "Wrong input file specifier (\'" << token << "\' instead of \'trunk\')." << std::endl;
		Simulation::exit(1);
	}

	if(strtoul(inputversion.c_str(), NULL, 0) < 20080701) {
		global_log->error() << "Input version tool old (" << inputversion << ")" << std::endl;
		Simulation::exit(1);
	}

	global_log->info() << "Reading phase space header from file " << _phaseSpaceHeaderFile << std::endl;


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
		global_log->info() << "{{" << token << "}}" << std::endl;

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
		} else if((token == "MoleculeFormat") || (token == "M")) {

			std::string ntypestring("ICRVQD");

			_phaseSpaceFileStream >> ntypestring;
			ntypestring.erase(ntypestring.find_last_not_of(" \t\n") + 1);
			ntypestring.erase(0, ntypestring.find_first_not_of(" \t\n"));

			if(!(ntypestring == "ICRVQD" || ntypestring == "ICRV"
				 || ntypestring == "IRV")) {
				global_log->error() << "Unknown molecule format: '"
									<< ntypestring << "'" << std::endl;
				Simulation::exit(1);
			}
			_moleculeFormat = ntypestring;
			global_log->info() << " molecule format: " << ntypestring << std::endl;
			header = false;
		} else if((token == "ThermostatTemperature") || (token == "ThT") || (token == "h")) {
			// set up a new thermostat
			int thermostat_id;
			double targetT;
			_phaseSpaceHeaderFileStream >> thermostat_id;
			_phaseSpaceHeaderFileStream >> targetT;
			global_log->info() << "Thermostat number " << thermostat_id << " has T = " << targetT << ".\n";
			domain->setTargetTemperature(thermostat_id, targetT);
		} else if((token == "ComponentThermostat") || (token == "CT") || (token == "o")) {
			// specify a thermostat for a component
			if(!domain->severalThermostats())
				domain->enableComponentwiseThermostat();
			int component_id;
			int thermostat_id;
			_phaseSpaceHeaderFileStream >> component_id >> thermostat_id;
			global_log->info() << "Component " << component_id << " (internally: "
							   << component_id - 1 << ") is regulated by thermostat number " << thermostat_id << ".\n";
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
			_simulation.getEnsemble()->domain() = new BoxDomain();
			for(int d = 0; d < 3; d++) {
				static_cast<BoxDomain*>(_simulation.getEnsemble()->domain())->setLength(d, globalLength[d]);
				domain->setGlobalLength(d, _simulation.getEnsemble()->domain()->length(d));
			}
		} else if((token == "HeatCapacity") || (token == "cv") || (token == "I")) {
			unsigned N;
			double U, UU;
			_phaseSpaceHeaderFileStream >> N >> U >> UU;
			domain->init_cv(N, U, UU);
		} else if((token == "NumberOfComponents") || (token == "C")) {
			// read in component definitions and
			// read in mixing coefficients

			// components:
			unsigned int numcomponents = 0;
			_phaseSpaceHeaderFileStream >> numcomponents;
			global_log->info() << "Reading " << numcomponents << " components" << std::endl;
			dcomponents.resize(numcomponents);
			for(unsigned int i = 0; i < numcomponents; i++) {
				global_log->info() << "comp. i = " << i << ": " << std::endl;
				dcomponents[i].setID(i);
				unsigned int numljcenters = 0;
				unsigned int numcharges = 0;
				unsigned int numdipoles = 0;
				unsigned int numquadrupoles = 0;
				unsigned int numtersoff = 0; //previously tersoff
				_phaseSpaceHeaderFileStream >> numljcenters >> numcharges >> numdipoles
											>> numquadrupoles >> numtersoff;
				if(numtersoff != 0) {
					global_log->error() << "tersoff no longer supported." << std::endl;
					Simulation::exit(-1);
				}
				double x, y, z, m;
				for(unsigned int j = 0; j < numljcenters; j++) {
					double eps, sigma, tcutoff, do_shift;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> m >> eps >> sigma >> tcutoff >> do_shift;
					dcomponents[i].addLJcenter(x, y, z, m, eps, sigma, tcutoff, (do_shift != 0));
					global_log->info() << "LJ at [" << x << " " << y << " " << z
									   << "], mass: " << m << ", epsilon: " << eps << ", sigma: " << sigma << std::endl;
				}
				for(unsigned int j = 0; j < numcharges; j++) {
					double q;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> m >> q;
					dcomponents[i].addCharge(x, y, z, m, q);
					global_log->info() << "charge at [" << x << " " << y << " " << z
									   << "], mass: " << m << ", q: " << q << std::endl;
				}
				for(unsigned int j = 0; j < numdipoles; j++) {
					double eMyx, eMyy, eMyz, absMy;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> eMyx >> eMyy >> eMyz >> absMy;
					dcomponents[i].addDipole(x, y, z, eMyx, eMyy, eMyz, absMy);
					global_log->info() << "dipole at [" << x << " " << y << " " << z << "] " << std::endl;
				}
				for(unsigned int j = 0; j < numquadrupoles; j++) {
					double eQx, eQy, eQz, absQ;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> eQx >> eQy >> eQz >> absQ;
					dcomponents[i].addQuadrupole(x, y, z, eQx, eQy, eQz, absQ);
					global_log->info() << "quad at [" << x << " " << y << " " << z << "] " << std::endl;
				}
				double IDummy1, IDummy2, IDummy3;
				// FIXME! Was soll das hier? Was ist mit der Initialisierung im Fall I <= 0.
				_phaseSpaceHeaderFileStream >> IDummy1 >> IDummy2 >> IDummy3;
				if(IDummy1 > 0.)
					dcomponents[i].setI11(IDummy1);
				if(IDummy2 > 0.)
					dcomponents[i].setI22(IDummy2);
				if(IDummy3 > 0.)
					dcomponents[i].setI33(IDummy3);
				domain->setProfiledComponentMass(dcomponents[i].m());
				global_log->info() << std::endl;
			}

#ifndef NDEBUG
			for(unsigned int i = 0; i < numcomponents; i++) {
				global_log->debug() << "Component " << (i + 1) << " of " << numcomponents << std::endl;
				global_log->debug() << dcomponents[i] << std::endl;
			}
#endif

			// Mixing coefficients
			std::vector<double>& dmixcoeff = domain->getmixcoeff();
			dmixcoeff.clear();
			for(unsigned int i = 1; i < numcomponents; i++) {
				for(unsigned int j = i + 1; j <= numcomponents; j++) {
					double xi, eta;
					_phaseSpaceHeaderFileStream >> xi >> eta;
					dmixcoeff.push_back(xi);
					dmixcoeff.push_back(eta);
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
			// This is not the last line of the header. The last line is 'MoleculeFormat'
			//header = false;
		} else if((token == "NumberOfMolecules") || (token == "N")) {
			// set number of Molecules
			// FIXME: Is this part called in any case as the token is handled in the readPhaseSpace method?
			// Yes, now it is needed, because we do not skip this line anymore.
			_phaseSpaceHeaderFileStream >> token;
			domain->setglobalNumMolecules( strtoul(token.c_str(),NULL,0) );
		}
		// LOCATION OF OLD PRESSURE GRADIENT TOKENS
		else {
			global_log->error() << "Invalid token \'" << token << "\' found. Skipping rest of the header." << std::endl;
			header = false;
		}
	}

	_simulation.getEnsemble()->setComponentLookUpIDs();

	_phaseSpaceHeaderFileStream.close();
}

unsigned long
MPI_IOReader::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {
#ifdef ENABLE_MPI
	Timer inputTimer;
	inputTimer.start();

	std::string token;
	std::vector<Component>& dcomponents = *(_simulation.getEnsemble()->getComponents());
	unsigned int numcomponents = dcomponents.size();
	unsigned long localMaxid = 0; // stores the highest molecule ID found in the phase space file

	if (numcomponents < 1) {
		global_log->warning()
				<< "No components defined! Setting up single one-centered LJ"
				<< std::endl;
		numcomponents = 1;
		dcomponents.resize(numcomponents);
		dcomponents[0].setID(0);
		dcomponents[0].addLJcenter(0., 0., 0., 1., 1., 1., 6., false);
	}

	if (domainDecomp->getRank() == 0) {
		global_log->info() << "Opening phase space file " << _phaseSpaceFile
				<< std::endl;
	}

	const char * fileName = _phaseSpaceFile.c_str();
	int ret, size;
	MPI_Status status;
	MPI_Offset header_offset = 0;

	MPI_File fh;
	MPI_Info info;

	MPI_Info_create(&info);

	/*
	timeval timer1, timer2;
	double timeDiff;

	if (domainDecomp->getRank() == 0) {
		gettimeofday(&timer1, NULL);
	}
	*/

	ret = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, info, &fh);
	if (ret != MPI_SUCCESS) {
		std::cerr << "Could not open phaseSpaceFile " << _phaseSpaceFile
				<< std::endl;
		handle_error(ret);
	}

	//numCellsAndMolecules[0] - number of cells
	//numCellsAndMolecules[1] - number of molecules
	long numCellsAndMolecules[2];
	double oldCellLength[3];

	//read the number of cells and the number of molecules
	if (domainDecomp->getRank() == 0) {
		ret = MPI_File_seek(fh, 0, MPI_SEEK_SET);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}
		ret = MPI_File_read(fh, numCellsAndMolecules, 2, MPI_LONG, &status);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}
	}

	ret = MPI_Bcast(numCellsAndMolecules, 2, MPI_LONG, 0, MPI_COMM_WORLD);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}

	domain->setglobalNumMolecules(numCellsAndMolecules[1]);

	//read the old cutoffRadius used for domain decomposition in the file

	ret = MPI_Type_size(MPI_LONG, &size);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}
	header_offset = 2 * size;

	if (domainDecomp->getRank() == 0) {
		//set file pointer after the first two long-variables
		ret = MPI_File_seek(fh, header_offset, MPI_SEEK_SET);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}

		ret = MPI_File_read(fh, oldCellLength, 3, MPI_DOUBLE, &status);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}
	}

	ret = MPI_Bcast(oldCellLength, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}

	//read the number of particles which each cell contains
	std::vector<int> globalNumParticlesPerCell (numCellsAndMolecules[0]);

	ret = MPI_Type_size(MPI_DOUBLE, &size);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}
	header_offset += (3 * size);

	if (domainDecomp->getRank() == 0) {
		ret = MPI_File_seek(fh, header_offset, MPI_SEEK_SET);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}

		ret = MPI_File_read(fh, globalNumParticlesPerCell.data(),
				numCellsAndMolecules[0], MPI_INT, &status);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}
	}

	ret = MPI_Bcast(globalNumParticlesPerCell.data(), numCellsAndMolecules[0],
			MPI_INT, 0, MPI_COMM_WORLD);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}

	//from now on the header_offset gives us the position after the header
	ret = MPI_Type_size(MPI_INT, &size);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}
	header_offset += (numCellsAndMolecules[0] * size);

	/*
	if (domainDecomp->getRank() == 0) {
		gettimeofday(&timer2, NULL);
		timeDiff = timer2.tv_sec - timer1.tv_sec + (timer2.tv_usec
				- timer1.tv_usec) / 1.E6;
		global_log->info() << "Das Lesen des MPI-IO Headers hat " << timeDiff
				<< " Sekunden benötigt" << std::endl;
	}
	*/

	//begin reading the cells

	//compute the length in cells. The cell size is given through oldCellLength
	int lengthInCells[3];
	for (unsigned short i = 0; i < 3; i++) {
		lengthInCells[i] = floor(
				domain->getGlobalLength(i) / oldCellLength[i]);

		mardyn_assert(lengthInCells[i] >= 0);
	}

	//compute bounding boxes in cells. The cellsize is determinded through the old
	//cellLength, which was used for the cellsize in the file
	int boundingBoxMinInCells[3];
	int boundingBoxMaxInCells[3];

	for (unsigned short i = 0; i < 3; i++) {
		boundingBoxMinInCells[i]
				= floor(
						domainDecomp->getBoundingBoxMin(i, domain)
								/ oldCellLength[i]);
		boundingBoxMaxInCells[i] = ceil(
				domainDecomp->getBoundingBoxMax(i, domain) / oldCellLength[i]
						- 1);

		mardyn_assert(boundingBoxMinInCells[i] >= 0);
		mardyn_assert(boundingBoxMaxInCells[i] >= 0);
		mardyn_assert(boundingBoxMinInCells[i] <= boundingBoxMaxInCells[i]);
	}

	int maxNumCellsLocal = (boundingBoxMaxInCells[2] - boundingBoxMinInCells[2]
			+ 1) * (boundingBoxMaxInCells[1] - boundingBoxMinInCells[1] + 1)
			* (boundingBoxMaxInCells[0] - boundingBoxMinInCells[0] + 1);

	//setup the elementary Datatype
	//it's the same, which is used for sending particles
	MPI_Datatype mpiParticleData;
	ParticleData::getMPIType(mpiParticleData);

	ret = MPI_Type_size(mpiParticleData, &size);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}


	//gettimeofday(&timer1, NULL);

	//setup the iterators
	int n = boundingBoxMinInCells[0];
	int k = boundingBoxMinInCells[1];
	int j = boundingBoxMinInCells[2];

	for (int i = 0; i < maxNumCellsLocal; i++) {

		int index = (j * lengthInCells[1] + k) * lengthInCells[0] + n;

		if (index >= (lengthInCells[2] * lengthInCells[1] * lengthInCells[0])) {
			break;
		}

		long offset = header_offset;

		for (int m = 0; m < index; m++) {
			offset += (globalNumParticlesPerCell[m] * size);
		}

		ret = MPI_File_seek(fh, offset, MPI_SEEK_SET);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}

		std::vector<ParticleData> data(globalNumParticlesPerCell[index]);

		ret = MPI_File_read(fh, data.data(), globalNumParticlesPerCell[index],
				mpiParticleData, &status);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}

		if (globalNumParticlesPerCell[index] > 0) {

			Molecule m;
			for (int l = 0; l < globalNumParticlesPerCell[index]; l++) {

				ParticleData::ParticleDataToMolecule(data[l], m);

				// only add particle if it is inside of the own domain!
				if(particleContainer->isInBoundingBox(m.r_arr().data())) {
					particleContainer->addParticle(m, true, false);
				}
				mardyn_assert(m.mass() != 0);

				// TODO: The following should be done by the addPartice method.
				if (m.r(0) >= domainDecomp->getBoundingBoxMin(0, domain)
						&& m.r(1)
								>= domainDecomp->getBoundingBoxMin(1, domain)
						&& m.r(2)
								>= domainDecomp->getBoundingBoxMin(2, domain)
						&& m.r(0)
								<= domainDecomp->getBoundingBoxMax(0, domain)
						&& m.r(1)
								<= domainDecomp->getBoundingBoxMax(1, domain)
						&& m.r(2)
								<= domainDecomp->getBoundingBoxMax(2, domain)) {
					dcomponents[m.componentid()].incNumMolecules();
				}

				if (m.getID() > localMaxid) {
					localMaxid = m.getID();
				}

				// Only called inside GrandCanonical
				global_simulation->getEnsemble()->storeSample(&m, m.componentid());
			}
		}

		//iterate over the boundary box of process

		//iterate over the x dimension
		if (n <= boundingBoxMaxInCells[0]) {
			if (n == boundingBoxMaxInCells[0]) {
				n = boundingBoxMinInCells[0];
			} else {
				n++;
			}
		}
		//iterate over the y dimension
		if (k <= boundingBoxMaxInCells[1] && n == boundingBoxMinInCells[0]) {
			if (k == boundingBoxMaxInCells[1]) {
				k = boundingBoxMinInCells[1];
			} else {
				k++;
			}
		}
		//iterate the y dimension
		if (n == boundingBoxMinInCells[0] && k == boundingBoxMinInCells[1]) {
			j++;
		}
	}

	/*
	gettimeofday(&timer2, NULL);
	timeDiff = timer2.tv_sec - timer1.tv_sec
			+ (timer2.tv_usec - timer1.tv_usec) / 1.E6;
	double timeDiffGlobal = 0;
	MPI_Reduce(&timeDiff, &timeDiffGlobal, 1, MPI_DOUBLE, MPI_MAX, 0,
				MPI_COMM_WORLD);
	if (domainDecomp->getRank() == 0) {

		global_log->info() << "Das Lesen der Zellen hat " << timeDiffGlobal
				<< " Sekunden benötigt" << std::endl;
	}
	*/

	ret = MPI_File_close(&fh);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}

	ret = MPI_Type_free(&mpiParticleData);
	if (ret != MPI_SUCCESS) {
		handle_error(ret);
	}

	/*
	if (domainDecomp->getRank() == 0) {
		gettimeofday(&timer1, NULL);
	}
	*/

	std::vector<long> numComponentMoleculesLocal(numcomponents);
	std::vector<long> numComponentMoleculesGlobal(numcomponents);

	for (unsigned int i = 0; i < numcomponents; i++) {
		numComponentMoleculesLocal[i] = dcomponents[i].getNumMolecules();
		numComponentMoleculesGlobal[i] = 0;
	}

	MPI_Allreduce(numComponentMoleculesLocal.data(), numComponentMoleculesGlobal.data(),
			numcomponents, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	for (unsigned int i = 0; i < numcomponents; i++) {
		dcomponents[i].setNumMolecules(numComponentMoleculesGlobal[i]);
		//for each molecule add dcomponents[i].getRotationalDegreesOfFreedom() to the globalRotDOF
		domain->setglobalRotDOF(
				dcomponents[i].getRotationalDegreesOfFreedom()
						* numComponentMoleculesGlobal[i]);
	}

	/*
	if (domainDecomp->getRank() == 0) {
		gettimeofday(&timer2, NULL);
		timeDiff = timer2.tv_sec - timer1.tv_sec + (timer2.tv_usec
				- timer1.tv_usec) / 1.E6;
		global_log->info() << "NumMolsInComps und globalRotDOF hat "
				<< timeDiff << " Sekunden benötigt" << std::endl;
	}
	*/

	global_log->info() << "Finished reading molecules: 100%" << std::endl;
	global_log->info() << "Reading Molecules done" << std::endl;

	// TODO: Shouldn't we always calculate this?
	if (domain->getglobalRho() < 1e-5) {
		domain->setglobalRho(
				domain->getglobalNumMolecules(true, particleContainer, domainDecomp) / domain->getGlobalVolume());
		global_log->info() << "Calculated Rho_global = "
				<< domain->getglobalRho() << std::endl;
	}
	//get maximum I/O time of each process and output it
	inputTimer.stop();

	double ioTime = inputTimer.get_etime();
	double maxIOTime = 1;

	MPI_Reduce(&ioTime, &maxIOTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (domainDecomp->getRank() == 0) {
		global_log->info() << "Initial IO took:                 " << maxIOTime
				<< " sec" << std::endl;
	}

	//get the global maximumID
	unsigned long globalMaxid = numCellsAndMolecules[1] - 1;
	MPI_Reduce(&localMaxid, &globalMaxid, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0,
			MPI_COMM_WORLD);

	return globalMaxid;
#else
	return 0;
#endif
}

void MPI_IOReader::handle_error(int i) {
#ifdef ENABLE_MPI
	char error_string[BUFSIZ];
	int length_of_error_string;

	MPI_Error_string(i, error_string, &length_of_error_string);

	global_log->error() << "Writing of file was not successfull " << " , " << i
			<< " , " << error_string << std::endl;
	Simulation::exit(1);
#endif
}
