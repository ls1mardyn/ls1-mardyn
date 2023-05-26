/*
 * BinaryReader.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: andal
 */

#include "BinaryReader.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "Domain.h"
#include "ensemble/BoxDomain.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include "molecules/Molecule.h"

#ifdef ENABLE_MPI
#include "parallel/ParticleData.h"
#include "parallel/DomainDecompBase.h"
#endif

#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/Timer.h"
#include "utils/xmlfileUnits.h"

#include <climits>
#include <cstdint>
#include <string>


enum MoleculeFormat : uint32_t {
	ICRVQD, IRV, ICRV
};

BinaryReader::BinaryReader()
		: _nMoleculeFormat(ICRVQD) {
	// TODO Auto-generated constructor stub
}

BinaryReader::~BinaryReader() = default;

void BinaryReader::readXML(XMLfileUnits& xmlconfig) {
	std::string pspfile;
	std::string pspheaderfile;
	xmlconfig.getNodeValue("header", pspheaderfile);
    pspheaderfile = string_utils::trim(pspheaderfile);
    if (pspheaderfile[0] != '/') {
      pspheaderfile.insert(0, xmlconfig.getDir());
    }
	Log::global_log->info() << "phase space header file: " << pspheaderfile << std::endl;
	xmlconfig.getNodeValue("data", pspfile);
	pspfile = string_utils::trim(pspfile);
	// only prefix xml dir if path is not absolute
	if (pspfile[0] != '/') {
	  pspfile.insert(0, xmlconfig.getDir());
	}
	Log::global_log->info() << "phase space data file: " << pspfile << std::endl;
	setPhaseSpaceHeaderFile(pspheaderfile);
	setPhaseSpaceFile(pspfile);
}

void BinaryReader::setPhaseSpaceFile(std::string filename) {
	_phaseSpaceFile = filename;
}

void BinaryReader::setPhaseSpaceHeaderFile(std::string filename) {
	_phaseSpaceHeaderFile = filename;
}

void BinaryReader::readPhaseSpaceHeader(Domain* domain, double timestep) {
	XMLfileUnits inp(_phaseSpaceHeaderFile);

	if(not inp.changecurrentnode("/mardyn")) {
		Log::global_log->error() << "Could not find root node /mardyn in XML input file." << std::endl;
		Log::global_log->fatal() << "Not a valid MarDyn XML input file." << std::endl;
		Simulation::exit(1);
	}

	bool bInputOk = true;
	double dCurrentTime = 0.;
	double dBoxLength[3] = {0., 0., 0.};
	uint64_t numMolecules = 0;
	std::string strMoleculeFormat;
	bInputOk = bInputOk && inp.changecurrentnode("headerinfo");
	bInputOk = bInputOk && inp.getNodeValue("time", dCurrentTime);
	bInputOk = bInputOk && inp.getNodeValue("length/x", dBoxLength[0]);
	bInputOk = bInputOk && inp.getNodeValue("length/y", dBoxLength[1]);
	bInputOk = bInputOk && inp.getNodeValue("length/z", dBoxLength[2]);
	bInputOk = bInputOk && inp.getNodeValue("number", numMolecules);
	bInputOk = bInputOk && inp.getNodeValue("format@type", strMoleculeFormat);

	if(not bInputOk) {
		Log::global_log->error() << "Content of file: '" << _phaseSpaceHeaderFile << "' corrupted! Program exit ..." << std::endl;
		Simulation::exit(1);
	}

	if("ICRVQD" == strMoleculeFormat)
		_nMoleculeFormat = ICRVQD;
	else if("IRV" == strMoleculeFormat)
		_nMoleculeFormat = IRV;
	else if("ICRV" == strMoleculeFormat)
		_nMoleculeFormat = ICRV;
	else {
		Log::global_log->error() << "Not a valid molecule format: " << strMoleculeFormat << ", program exit ..." << std::endl;
		Simulation::exit(1);
	}

	// Set parameters of Domain and Simulation class
	_simulation.setSimulationTime(dCurrentTime);
	for(uint8_t d = 0; d < 3; ++d) {
		domain->setGlobalLength(d, dBoxLength[d]);
	}
	domain->setglobalNumMolecules(numMolecules);
}

unsigned long
BinaryReader::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {

	Timer inputTimer;
	inputTimer.start();

#ifdef ENABLE_MPI
	if (domainDecomp->getRank() == 0) { // Rank 0 only
#endif
	Log::global_log->info() << "Opening phase space file " << _phaseSpaceFile
					   << std::endl;
	_phaseSpaceFileStream.open(_phaseSpaceFile.c_str(),
							   std::ios::binary | std::ios::in);
	if(!_phaseSpaceFileStream.is_open()) {
		Log::global_log->error() << "Could not open phaseSpaceFile "
							<< _phaseSpaceFile << std::endl;
		Simulation::exit(1);
	}
	Log::global_log->info() << "Reading phase space file " << _phaseSpaceFile
					   << std::endl;
#ifdef ENABLE_MPI
	} // Rank 0 only
#endif

	std::string token;
	std::vector<Component>& dcomponents = *(_simulation.getEnsemble()->getComponents());
	size_t numcomponents = dcomponents.size();
	unsigned long maxid = 0; // stores the highest molecule ID found in the phase space file

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

	double x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz;
	uint64_t id;
	uint32_t componentid = 0;

	x = y = z = vx = vy = vz = q1 = q2 = q3 = Dx = Dy = Dz = 0.;
	q0 = 1.;

	// Global number of particles must not be updated as this would result in numMolecules = 0
	uint64_t numMolecules = domain->getglobalNumMolecules(false);

	for(uint64_t i = 0; i < numMolecules; i++) {

#ifdef ENABLE_MPI
		if (domainDecomp->getRank() == 0) { // Rank 0 only
#endif
		if(_phaseSpaceFileStream.eof()) {
			Log::global_log->error() << "End of file was hit before all " << numMolecules << " expected molecules were read."
				<< std::endl;
			Simulation::exit(1);
        }
		_phaseSpaceFileStream.read(reinterpret_cast<char*> (&id), 8);
		switch (_nMoleculeFormat) {
			case ICRVQD:
				_phaseSpaceFileStream.read(
						reinterpret_cast<char*> (&componentid), 4);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&x), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&y), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&z), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vx), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vy), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vz), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&q0), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&q1), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&q2), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&q3), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&Dx), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&Dy), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&Dz), 8);
				break;
			case ICRV:
				_phaseSpaceFileStream.read(
						reinterpret_cast<char*> (&componentid), 4);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&x), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&y), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&z), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vx), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vy), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vz), 8);
				break;
			case IRV:
				componentid = 1;
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&x), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&y), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&z), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vx), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vy), 8);
				_phaseSpaceFileStream.read(reinterpret_cast<char*> (&vz), 8);
				break;
			default:
				Log::global_log->error() << "BinaryReader: Unknown phase space format: " << _nMoleculeFormat << std::endl
									<< "Aborting simulation." << std::endl;
				Simulation::exit(12);
		}
		if ((x < 0.0 || x >= domain->getGlobalLength(0)) || (y < 0.0 || y >= domain->getGlobalLength(1)) ||
			(z < 0.0 || z >= domain->getGlobalLength(2))) {
			Log::global_log->warning() << "Molecule " << id << " out of box: " << x << ";" << y << ";" << z << std::endl;
		}

		if(componentid > numcomponents) {
			Log::global_log->error() << "Molecule id " << id
								<< " has a component ID greater than the existing number of components: "
								<< componentid
								<< ">"
								<< numcomponents << std::endl;
			Simulation::exit(1);
		}
		if(componentid == 0) {
			Log::global_log->error() << "Molecule id " << id
								<< " has componentID == 0." << std::endl;
			Simulation::exit(1);
		}
		// ComponentIDs are used as array IDs, hence need to start at 0.
		// In the input files they always start with 1 so we need to adapt that all the time.
		componentid--;

		// store only those molecules within the domain of this process
		// The neccessary check is performed in the particleContainer addPartice method
		// FIXME: Datastructures? Pass pointer instead of object, so that we do not need to copy?!
		Molecule m1 = Molecule(id, &dcomponents[componentid], x, y, z, vx,
							   vy, vz, q0, q1, q2, q3, Dx, Dy, Dz);
#ifdef ENABLE_MPI
		ParticleData::MoleculeToParticleData(
				particle_buff[particle_buff_pos], m1);
	} // Rank 0 only

	particle_buff_pos++;
	if ((particle_buff_pos >= PARTICLE_BUFFER_SIZE) || (i
			== numMolecules - 1)) {
		MPI_Bcast(&particle_buff_pos, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, 0,
				MPI_COMM_WORLD); // TODO: MPI_COMM_WORLD
		for (int j = 0; j < particle_buff_pos; j++) {
			Molecule m;
			ParticleData::ParticleDataToMolecule(particle_buff[j], m);
			// only add particle if it is inside of the own domain!
			if(particleContainer->isInBoundingBox(m.r_arr().data())) {
				particleContainer->addParticle(m, true, false);
			}

			dcomponents[m.componentid()].incNumMolecules();
			domain->setglobalRotDOF(
					dcomponents[m.componentid()].getRotationalDegreesOfFreedom()
							+ domain->getglobalRotDOF());

			if (m.getID() > maxid)
				maxid = m.getID();

			// Only called inside GrandCanonical
			global_simulation->getEnsemble()->storeSample(&m, componentid);
		}
		particle_buff_pos = 0;
	}
#else
		if (particleContainer->isInBoundingBox(m1.r_arr().data())) {
			particleContainer->addParticle(m1, true, false);
		}

		// TODO: The following should be done by the addPartice method.
		dcomponents[componentid].incNumMolecules();
		domain->setglobalRotDOF(
				dcomponents[componentid].getRotationalDegreesOfFreedom()
				+ domain->getglobalRotDOF());

		if(id > maxid)
			maxid = id;

		// Only called inside GrandCanonical
		global_simulation->getEnsemble()->storeSample(&m1, componentid);
#endif

		// Print status message
		unsigned long iph = numMolecules / 100;
		if(iph != 0 && (i % iph) == 0)
			Log::global_log->info() << "Finished reading molecules: " << i / iph
							   << "%\r" << std::flush;
	}

	Log::global_log->info() << "Finished reading molecules: 100%" << std::endl;
	Log::global_log->info() << "Reading Molecules done" << std::endl;

	// TODO: Shouldn't we always calculate this?
	if (domain->getglobalRho() < 1e-5) {
		domain->setglobalRho(
				domain->getglobalNumMolecules(true, particleContainer, domainDecomp) / domain->getGlobalVolume());
		Log::global_log->info() << "Calculated Rho_global = "
						   << domain->getglobalRho() << std::endl;
	}

#ifdef ENABLE_MPI
	if (domainDecomp->getRank() == 0) { // Rank 0 only
#endif
	_phaseSpaceFileStream.close();
#ifdef ENABLE_MPI
	} // Rank 0 only
#endif

	inputTimer.stop();
	Log::global_log->info() << "Initial IO took:                 "
					   << inputTimer.get_etime() << " sec" << std::endl;
#ifdef ENABLE_MPI
	MPI_CHECK(MPI_Type_free(&mpi_Particle));
#endif
	return maxid;
}
