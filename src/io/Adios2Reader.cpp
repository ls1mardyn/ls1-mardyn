/*
 * Adios2Reader.cpp
 *
 *  Created on: April 2021
 *      Author: Tobias Rau, Matthias Heinen, Patrick Gralka, Christoph Niethammer, Simon Homes
 */

#ifdef ENABLE_ADIOS2

#include "Adios2Reader.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "parallel/DomainDecompBase.h"
#ifdef ENABLE_MPI
#include "parallel/ParticleData.h"
#endif
#include "utils/xmlfile.h"

using Log::global_log;

void Adios2Reader::init(ParticleContainer* particleContainer,
        DomainDecompBase* domainDecomp, Domain* domain) {
};

void Adios2Reader::readXML(XMLfileUnits& xmlconfig) {
	_mode = "rootOnly";
	xmlconfig.getNodeValue("mode", _mode);
	global_log->info() << "[Adios2Reader] Input mode: " << _mode << endl;
	_inputfile = "mardyn.bp";
	xmlconfig.getNodeValue("filename", _inputfile);
	global_log->info() << "[Adios2Reader] Inputfile: " << _inputfile << endl;
	_adios2enginetype = "BP4";
	xmlconfig.getNodeValue("adios2enginetype", _adios2enginetype);
	global_log->info() << "[Adios2Reader] Adios2 engine type: " << _adios2enginetype << endl;
	_step = -1;
	xmlconfig.getNodeValue("adios2Step", _step);
	global_log->info() << "[Adios2Reader] step to load from input file: " << _step << endl;

	if (!mainInstance) initAdios2();
};

void Adios2Reader::testInit(std::string infile, int step, std::string adios2enginetype, std::string mode) {
	using std::endl;
	_inputfile = infile;
	global_log->info() << "[Adios2Reader] Inputfile: " << _inputfile << endl;
	_adios2enginetype = adios2enginetype;
	global_log->info() << "[Adios2Reader] Adios2 engine type: " << _adios2enginetype << endl;
	_step = step;
	global_log->info() << "[Adios2Reader] step to load from input file: " << _step << endl;
	_mode = mode;
	global_log->info() << "[Adios2Reader] Input mode: " << _mode << endl;

	if (!mainInstance) initAdios2();
}

void Adios2Reader::readPhaseSpaceHeader(Domain* domain, double timestep) {
//EMPTY
};


unsigned long Adios2Reader::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain,
										   DomainDecompBase* domainDecomp) {
	auto variables = io->AvailableVariables();

	auto total_steps = std::stoi(variables["simulationtime"]["AvailableStepsCount"]);
	global_log->info() << "[Adios2Reader] Total Steps in adios file " << total_steps << std::endl;

	if (_step == -1) {
		_step = total_steps - 1;
	}
	if (_step > total_steps) {
		global_log->error() << "[Adios2Reader] Specified step is out of scope" << std::endl;
	}
	if (_step < 0) {
		_step = total_steps + (_step + 1);
	}
	particle_count = std::stoi(variables["rx"]["Shape"]);
	global_log->info() << "[Adios2Reader] Particle count: " << particle_count << std::endl;

	Timer inputTimer;
	inputTimer.start();

	if (_mode == "rootOnly") {
		rootOnlyRead(particleContainer, domain, domainDecomp);
#ifdef ENABLE_MPI
		MPI_Bcast(&_simtime, 1, MPI_DOUBLE, 0, domainDecomp->getCommunicator());
#endif
	} else if (_mode == "parallelRead") {
		parallelRead(particleContainer, domain, domainDecomp);
	} else {
		global_log->error() << "[Adios2Reader] Unknown _mode '" << _mode << "'" << std::endl;
	}

	_simulation.setSimulationTime(_simtime);
	global_log->info() << "[Adios2Reader] simulation time is: " << _simtime << std::endl;

	global_log->info() << "[Adios2Reader] Finished reading molecules: 100%" << std::endl;
	global_log->info() << "[Adios2Reader] Reading Molecules done" << std::endl;

	inputTimer.stop();
	global_log->info() << "[Adios2Reader] Initial IO took: " << inputTimer.get_etime() << " sec" << std::endl;

	if (domain->getglobalRho() == 0.) {
		domain->setglobalRho(domain->getglobalNumMolecules(true, particleContainer, domainDecomp) / domain->getGlobalVolume());
		global_log->info() << "[Adios2Reader] Calculated Rho_global = " << domain->getglobalRho() << endl;
	}

	engine->Close();
	global_log->info() << "[Adios2Reader] finish." << std::endl;

	return particle_count;
};

void Adios2Reader::rootOnlyRead(ParticleContainer* particleContainer, Domain* domain,
								  DomainDecompBase* domainDecomp) {

	std::variant<std::vector<float>, std::vector<double>> rx, ry, rz, vx, vy, vz, qw, qx, qy, qz, Lx, Ly, Lz;
	std::vector<uint64_t> mol_id, comp_id;
	// Root sends bufferSize particles in one MPI_Bcast to all nodes
	uint64_t bufferSize = 1024;
#ifdef ENABLE_MPI
	MPI_Datatype mpi_Particle;
	ParticleData::getMPIType(mpi_Particle);
#endif

	auto num_reads = particle_count / bufferSize;
	if (particle_count % bufferSize != 0) num_reads += 1;
	global_log->info() << "[Adios2Reader] Input is divided into " << num_reads << " sequential reads." << endl;

	auto variables = io->AvailableVariables();

	for (const auto &var : variables) {
		if (var.first == "rx") {
			if (var.second.at("Type") != "double") {
				global_log->info() << "[Adios2Reader] Detected single precision" << endl;
				_single_precision = true;
				rx = std::vector<float>();
				ry = std::vector<float>();
				rz = std::vector<float>();
				vx = std::vector<float>();
				vy = std::vector<float>();
				vz = std::vector<float>();
				qw = std::vector<float>();
				qx = std::vector<float>();
				qy = std::vector<float>();
				qz = std::vector<float>();
				Lx = std::vector<float>();
				Ly = std::vector<float>();
				Lz = std::vector<float>();
			} else {
				global_log->info() << "[Adios2Reader] Detected double precision" << endl;
				rx = std::vector<double>();
				ry = std::vector<double>();
				rz = std::vector<double>();
				vx = std::vector<double>();
				vy = std::vector<double>();
				vz = std::vector<double>();
				qw = std::vector<double>();
				qx = std::vector<double>();
				qy = std::vector<double>();
				qz = std::vector<double>();
				Lx = std::vector<double>();
				Ly = std::vector<double>();
				Lz = std::vector<double>();
			}
		}
	}

	for (int read = 0; read < num_reads; read++) {
		global_log->info() << "[Adios2Reader] Performing read " << read << endl;
		const uint64_t offset = read * bufferSize;
		if (read == num_reads - 1) bufferSize = particle_count % bufferSize;
		if (domainDecomp->getRank() == 0) {
			if (_single_precision) {
				performInquire(variables, bufferSize, offset, std::get<std::vector<float>>(rx),
							   std::get<std::vector<float>>(ry), std::get<std::vector<float>>(rz),
							   std::get<std::vector<float>>(vx), std::get<std::vector<float>>(vy),
							   std::get<std::vector<float>>(vz), std::get<std::vector<float>>(qw),
							   std::get<std::vector<float>>(qx), std::get<std::vector<float>>(qy),
							   std::get<std::vector<float>>(qz), std::get<std::vector<float>>(Lx),
							   std::get<std::vector<float>>(Ly), std::get<std::vector<float>>(Lz), mol_id, comp_id);
			} else {
				performInquire(variables, bufferSize, offset, std::get<std::vector<double>>(rx),
							   std::get<std::vector<double>>(ry), std::get<std::vector<double>>(rz),
							   std::get<std::vector<double>>(vx), std::get<std::vector<double>>(vy),
							   std::get<std::vector<double>>(vz), std::get<std::vector<double>>(qw),
							   std::get<std::vector<double>>(qx), std::get<std::vector<double>>(qy),
							   std::get<std::vector<double>>(qz), std::get<std::vector<double>>(Lx),
							   std::get<std::vector<double>>(Ly), std::get<std::vector<double>>(Lz), mol_id, comp_id);
			}
		}

		engine->PerformGets();
		global_log->info() << "[Adios2Reader] Read " << read << " done." << endl;

		if (_simulation.getEnsemble()->getComponents()->empty()) {
			auto attributes = io->AvailableAttributes();
			auto comp_id_copy = comp_id;
			std::sort(comp_id_copy.begin(), comp_id_copy.end());
			auto u_iter = std::unique(comp_id_copy.begin(), comp_id_copy.end());
			comp_id_copy.resize(std::distance(comp_id_copy.begin(), u_iter));
			_dcomponents.resize(comp_id_copy.size());

			for (auto comp : comp_id_copy) {
				std::string name = "component_" + std::to_string(comp);
				auto centers_var = io->InquireAttribute<double>(name + "_centers");
				auto centers = centers_var.Data();

				auto sigma_var = io->InquireAttribute<double>(name + "_sigma");
				auto sigma = sigma_var.Data();

				auto mass_var = io->InquireAttribute<double>(name + "_mass");
				auto mass = mass_var.Data();

				auto eps_var = io->InquireAttribute<double>(name + "_epsilon");
				auto eps = eps_var.Data();

				auto cname_var = io->InquireAttribute<std::string>(name + "_name");
				auto cname = cname_var.Data();

				_dcomponents[comp].addLJcenter(centers[0], centers[0], centers[0], mass[0], eps[0], sigma[0]);
				_dcomponents[comp].setName(cname[0]);
			}
		} else {
			_dcomponents = *(_simulation.getEnsemble()->getComponents());
		}
		global_log->info() << "[Adios2Reader] Gathered components." << std::endl;

#ifdef ENABLE_MPI
		std::vector<ParticleData> particle_buff(bufferSize);
		if (domainDecomp->getRank() == 0) {
			for (int i = 0; i < bufferSize; i++) {
				Molecule m1;
				if (_single_precision) {
					m1 = fillMolecule(i, mol_id, comp_id, std::get<std::vector<float>>(rx),
									  std::get<std::vector<float>>(ry), std::get<std::vector<float>>(rz),
									  std::get<std::vector<float>>(vx), std::get<std::vector<float>>(vy),
									  std::get<std::vector<float>>(vz), std::get<std::vector<float>>(qw),
									  std::get<std::vector<float>>(qx), std::get<std::vector<float>>(qy),
									  std::get<std::vector<float>>(qz), std::get<std::vector<float>>(Lx),
									  std::get<std::vector<float>>(Ly), std::get<std::vector<float>>(Lz));
				} else {
					m1 = fillMolecule(i, mol_id, comp_id, std::get<std::vector<double>>(rx),
									  std::get<std::vector<double>>(ry), std::get<std::vector<double>>(rz),
									  std::get<std::vector<double>>(vx), std::get<std::vector<double>>(vy),
									  std::get<std::vector<double>>(vz), std::get<std::vector<double>>(qw),
									  std::get<std::vector<double>>(qx), std::get<std::vector<double>>(qy),
									  std::get<std::vector<double>>(qz), std::get<std::vector<double>>(Lx),
									  std::get<std::vector<double>>(Ly), std::get<std::vector<double>>(Lz));
				}
				ParticleData::MoleculeToParticleData(particle_buff[i], m1);
			}
		}
		MPI_Bcast(particle_buff.data(), bufferSize, mpi_Particle, 0, domainDecomp->getCommunicator());

		for (int j = 0; j < bufferSize; j++) {
			Molecule m =
				Molecule(particle_buff[j].id, &_dcomponents[particle_buff[j].cid], particle_buff[j].r[0],
						 particle_buff[j].r[1], particle_buff[j].r[2], particle_buff[j].v[0], particle_buff[j].v[1],
						 particle_buff[j].v[2], particle_buff[j].q[0], particle_buff[j].q[1], particle_buff[j].q[2],
						 particle_buff[j].q[3], particle_buff[j].D[0], particle_buff[j].D[1], particle_buff[j].D[2]);

			// only add particle if it is inside of the own domain!
			if (particleContainer->isInBoundingBox(m.r_arr().data())) {
				particleContainer->addParticle(m, true, false);
			}

			_dcomponents[m.componentid()].incNumMolecules();
			domain->setglobalRotDOF(_dcomponents[m.componentid()].getRotationalDegreesOfFreedom() +
									domain->getglobalRotDOF());
		}
#else
		for (int i = 0; i < bufferSize; i++) {
			global_log->debug() << "[Adios2Reader] Processing particle " << offset + i << std::endl;
			Molecule m;
			if (_single_precision) {
				m = fillMolecule(i, mol_id, comp_id, std::get<std::vector<float>>(rx), std::get<std::vector<float>>(ry),
								 std::get<std::vector<float>>(rz), std::get<std::vector<float>>(vx),
								 std::get<std::vector<float>>(vy), std::get<std::vector<float>>(vz),
								 std::get<std::vector<float>>(qw), std::get<std::vector<float>>(qx),
								 std::get<std::vector<float>>(qy), std::get<std::vector<float>>(qz),
								 std::get<std::vector<float>>(Lx), std::get<std::vector<float>>(Ly),
								 std::get<std::vector<float>>(Lz));
			} else {
				m = fillMolecule(i, mol_id, comp_id, std::get<std::vector<double>>(rx),
								 std::get<std::vector<double>>(ry), std::get<std::vector<double>>(rz),
								 std::get<std::vector<double>>(vx), std::get<std::vector<double>>(vy),
								 std::get<std::vector<double>>(vz), std::get<std::vector<double>>(qw),
								 std::get<std::vector<double>>(qx), std::get<std::vector<double>>(qy),
								 std::get<std::vector<double>>(qz), std::get<std::vector<double>>(Lx),
								 std::get<std::vector<double>>(Ly), std::get<std::vector<double>>(Lz));
			}

			// only add particle if it is inside of the own domain!
			if (particleContainer->isInBoundingBox(m.r_arr().data())) {
				particleContainer->addParticle(m, true, false);
			}

			_dcomponents[m.componentid()].incNumMolecules();
			domain->setglobalRotDOF(_dcomponents[m.componentid()].getRotationalDegreesOfFreedom() +
									domain->getglobalRotDOF());
		}
#endif

		// Print status message
		unsigned long iph = num_reads / 100;
		if (iph != 0 && (read % iph) == 0)
			global_log->info() << "[Adios2Read] Finished reading molecules: " << read / iph << "%\r" << std::flush;
	}
}

void Adios2Reader::parallelRead(ParticleContainer* particleContainer, Domain* domain,
								  DomainDecompBase* domainDecomp) {

    std::variant<std::vector<float>, std::vector<double>> rx, ry, rz, vx, vy, vz, qw, qx, qy, qz, Lx, Ly, Lz;
	std::vector<uint64_t> mol_id, comp_id;
	uint64_t bufferSize = particle_count / domainDecomp->getNumProcs();
    auto variables = io->AvailableVariables();

	for (const auto &var : variables) {
		if (var.first == "rx") {
			if (var.second.at("Type") != "double") {
				global_log->info() << "[Adios2Reader] Detected single precision" << endl;
				_single_precision = true;
				rx = std::vector<float>();
				ry = std::vector<float>();
				rz = std::vector<float>();
				vx = std::vector<float>();
				vy = std::vector<float>();
				vz = std::vector<float>();
				qw = std::vector<float>();
				qx = std::vector<float>();
				qy = std::vector<float>();
				qz = std::vector<float>();
				Lx = std::vector<float>();
				Ly = std::vector<float>();
				Lz = std::vector<float>();
			} else {
				global_log->info() << "[Adios2Reader] Detected double precision" << endl;
				rx = std::vector<double>();
				ry = std::vector<double>();
				rz = std::vector<double>();
				vx = std::vector<double>();
				vy = std::vector<double>();
				vz = std::vector<double>();
				qw = std::vector<double>();
				qx = std::vector<double>();
				qy = std::vector<double>();
				qz = std::vector<double>();
				Lx = std::vector<double>();
				Ly = std::vector<double>();
				Lz = std::vector<double>();
			}
		}
	}

    uint64_t offset = domainDecomp->getRank() * bufferSize;
	if (domainDecomp->getRank() == domainDecomp->getNumProcs() - 1) {
		bufferSize = particle_count % bufferSize == 0 ? bufferSize : bufferSize + particle_count % bufferSize;
	}

	if (_single_precision) {
		performInquire(
			variables, bufferSize, offset, std::get<std::vector<float>>(rx), std::get<std::vector<float>>(ry),
			std::get<std::vector<float>>(rz), std::get<std::vector<float>>(vx), std::get<std::vector<float>>(vy),
			std::get<std::vector<float>>(vz), std::get<std::vector<float>>(qw), std::get<std::vector<float>>(qx),
			std::get<std::vector<float>>(qy), std::get<std::vector<float>>(qz), std::get<std::vector<float>>(Lx),
			std::get<std::vector<float>>(Ly), std::get<std::vector<float>>(Lz), mol_id, comp_id);
	} else {
		performInquire(
			variables, bufferSize, offset, std::get<std::vector<double>>(rx), std::get<std::vector<double>>(ry),
			std::get<std::vector<double>>(rz), std::get<std::vector<double>>(vx), std::get<std::vector<double>>(vy),
			std::get<std::vector<double>>(vz), std::get<std::vector<double>>(qw), std::get<std::vector<double>>(qx),
			std::get<std::vector<double>>(qy), std::get<std::vector<double>>(qz), std::get<std::vector<double>>(Lx),
			std::get<std::vector<double>>(Ly), std::get<std::vector<double>>(Lz), mol_id, comp_id);
	}

    engine->PerformGets();
	global_log->debug() << "[Adios2Reader] Processed gets." << endl;

	std::vector<uint64_t> global_component_ids{};
#ifdef ENABLE_MPI
	MPI_Datatype mpi_Particle;
	ParticleData::getMPIType(mpi_Particle);

	std::vector<int> counts_array(domainDecomp->getNumProcs(), 0);
	std::vector<int> displacements_array(domainDecomp->getNumProcs(), 0);
	counts_array[domainDecomp->getRank()] = bufferSize;
	displacements_array[domainDecomp->getRank()] = offset;
	std::vector<int> recvBuffer(domainDecomp->getNumProcs(), 0);

	MPI_Allgather(&counts_array[domainDecomp->getRank()], 1, MPI_INT, recvBuffer.data(), 1, MPI_INT, domainDecomp->getCommunicator());
	counts_array = recvBuffer;
	MPI_Allgather(&displacements_array[domainDecomp->getRank()], 1, MPI_INT, recvBuffer.data(), 1, MPI_INT,
				  domainDecomp->getCommunicator());
	displacements_array = recvBuffer;

	global_component_ids.resize(particle_count);
	MPI_Allgatherv(comp_id.data(), bufferSize, MPI_UINT64_T, global_component_ids.data(),
				   counts_array.data(),
				   displacements_array.data(), MPI_UINT64_T, domainDecomp->getCommunicator());
#endif

    if (_simulation.getEnsemble()->getComponents()->empty()) {
		global_log->debug() << "[Adios2Reader] getComponents is empty. Reading components from Adios file ..." << std::endl;
        auto attributes = io->AvailableAttributes();
		auto comp_id_copy = global_component_ids.empty() ? comp_id : global_component_ids;
		std::sort(comp_id_copy.begin(), comp_id_copy.end());
		auto u_iter = std::unique(comp_id_copy.begin(), comp_id_copy.end());
		comp_id_copy.resize(std::distance(comp_id_copy.begin(), u_iter));
		_dcomponents.resize(comp_id_copy.size());

    	for (auto comp : comp_id_copy) {
			std::string name = "component_" + std::to_string(comp);
			auto centers_var = io->InquireAttribute<double>(name + "_centers");
    		auto centers = centers_var.Data();

			auto sigma_var = io->InquireAttribute<double>(name + "_sigma");
			auto sigma = sigma_var.Data();

			auto mass_var = io->InquireAttribute<double>(name + "_mass");
			auto mass = mass_var.Data();

			auto eps_var = io->InquireAttribute<double>(name + "_epsilon");
			auto eps = eps_var.Data();

			auto cname_var = io->InquireAttribute<std::string>(name + "_name");
			auto cname = cname_var.Data();

    		_dcomponents[comp].addLJcenter(centers[0], centers[0], centers[0], mass[0], eps[0], sigma[0]);
			_dcomponents[comp].setName(cname[0]);
    	}
    } else {
		_dcomponents = *(_simulation.getEnsemble()->getComponents());
    }
	global_log->debug() << "[Adios2Reader] Gathered components." << std::endl;


#ifdef ENABLE_MPI
	std::vector<ParticleData> particle_buff(particle_count);
	for (int i = 0; i < bufferSize; i++) {
    	Molecule m1;
		if (_single_precision) {
			m1 = fillMolecule(
				i, mol_id, comp_id, std::get<std::vector<float>>(rx), std::get<std::vector<float>>(ry),
				std::get<std::vector<float>>(rz), std::get<std::vector<float>>(vx), std::get<std::vector<float>>(vy),
				std::get<std::vector<float>>(vz), std::get<std::vector<float>>(qw), std::get<std::vector<float>>(qx),
				std::get<std::vector<float>>(qy), std::get<std::vector<float>>(qz), std::get<std::vector<float>>(Lx),
				std::get<std::vector<float>>(Ly), std::get<std::vector<float>>(Lz));
		} else {
			m1 = fillMolecule(
				i, mol_id, comp_id, std::get<std::vector<double>>(rx), std::get<std::vector<double>>(ry),
				std::get<std::vector<double>>(rz), std::get<std::vector<double>>(vx), std::get<std::vector<double>>(vy),
				std::get<std::vector<double>>(vz), std::get<std::vector<double>>(qw), std::get<std::vector<double>>(qx),
				std::get<std::vector<double>>(qy), std::get<std::vector<double>>(qz), std::get<std::vector<double>>(Lx),
				std::get<std::vector<double>>(Ly), std::get<std::vector<double>>(Lz));
		}
		ParticleData::MoleculeToParticleData(particle_buff[i + offset], m1);
    }

	global_log->debug() << "[Adios2Reader] performing allgather" << std::endl;
	MPI_Allgatherv(&particle_buff[offset], bufferSize, mpi_Particle, particle_buff.data(),
				  counts_array.data(), displacements_array.data(),
				  mpi_Particle,
				  domainDecomp->getCommunicator());
	for (int j = 0; j < particle_buff.size(); j++) {
		Molecule m = Molecule(particle_buff[j].id, &_dcomponents[particle_buff[j].cid], particle_buff[j].r[0],
							  particle_buff[j].r[1], particle_buff[j].r[2], particle_buff[j].v[0], particle_buff[j].v[1],
							  particle_buff[j].v[2], particle_buff[j].q[0], particle_buff[j].q[1], particle_buff[j].q[2],
							  particle_buff[j].q[3], particle_buff[j].D[0], particle_buff[j].D[1], particle_buff[j].D[2]);

        // only add particle if it is inside of the own domain!
        if(particleContainer->isInBoundingBox(m.r_arr().data())) {
            particleContainer->addParticle(m, true, false);
        }

        _dcomponents[m.componentid()].incNumMolecules();
        domain->setglobalRotDOF(
        _dcomponents[m.componentid()].getRotationalDegreesOfFreedom()
                + domain->getglobalRotDOF());
    }
#else
	for (int i = 0; i < particle_count; i++) {
		Molecule m;
		if (_single_precision) {
			m = fillMolecule(
				i, mol_id, comp_id, std::get<std::vector<float>>(rx), std::get<std::vector<float>>(ry),
				std::get<std::vector<float>>(rz), std::get<std::vector<float>>(vx), std::get<std::vector<float>>(vy),
				std::get<std::vector<float>>(vz), std::get<std::vector<float>>(qw), std::get<std::vector<float>>(qx),
				std::get<std::vector<float>>(qy), std::get<std::vector<float>>(qz), std::get<std::vector<float>>(Lx),
				std::get<std::vector<float>>(Ly), std::get<std::vector<float>>(Lz));
		} else {
			m = fillMolecule(
				i, mol_id, comp_id, std::get<std::vector<double>>(rx), std::get<std::vector<double>>(ry),
				std::get<std::vector<double>>(rz), std::get<std::vector<double>>(vx), std::get<std::vector<double>>(vy),
				std::get<std::vector<double>>(vz), std::get<std::vector<double>>(qw), std::get<std::vector<double>>(qx),
				std::get<std::vector<double>>(qy), std::get<std::vector<double>>(qz), std::get<std::vector<double>>(Lx),
				std::get<std::vector<double>>(Ly), std::get<std::vector<double>>(Lz));
		}

		// only add particle if it is inside of the own domain!
		if (particleContainer->isInBoundingBox(m.r_arr().data())) {
			particleContainer->addParticle(m, true, false);
		}

		_dcomponents[m.componentid()].incNumMolecules();
		domain->setglobalRotDOF(_dcomponents[m.componentid()].getRotationalDegreesOfFreedom() +
								domain->getglobalRotDOF());
	}
#endif


}

void Adios2Reader::initAdios2() {
  auto& domainDecomp = _simulation.domainDecomposition();

  try {
    //get adios2 instance
#ifdef ENABLE_MPI
		mainInstance = std::make_shared<adios2::ADIOS>(domainDecomp.getCommunicator());
#else
		mainInstance = std::make_shared<adios2::ADIOS>();
#endif
		io = std::make_shared<adios2::IO>(mainInstance->DeclareIO("Input"));

        io->SetEngine(_adios2enginetype);

        if (!engine) {
            global_log->info() << "[Adios2Reader] Opening File for reading: " << _inputfile.c_str() << std::endl;
			engine = std::make_shared<adios2::Engine>(io->Open(_inputfile, adios2::Mode::Read));
        }
  }
    catch (std::invalid_argument& e) {
        global_log->fatal()
                << "[Adios2Reader] Invalid argument exception, STOPPING PROGRAM from rank"
                << domainDecomp.getRank()
                << ": " << e.what() << std::endl;
	mardyn_exit(1);
    } catch (std::ios_base::failure& e) {
        global_log->fatal()
                << "[Adios2Reader] IO System base failure exception, STOPPING PROGRAM from rank "
                << domainDecomp.getRank()
                << ": " << e.what() << std::endl;
	mardyn_exit(1);
    } catch (std::exception& e) {
        global_log->fatal()
                << "[Adios2Reader] Exception, STOPPING PROGRAM from rank"
                << domainDecomp.getRank()
                << ": " << e.what() << std::endl;
	mardyn_exit(1);
    }
    global_log->info() << "[Adios2Reader] Init complete." << std::endl;
};

#endif // ENABLE_ADIOS2
