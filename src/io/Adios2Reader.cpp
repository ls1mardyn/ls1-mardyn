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

	if (!inst) initAdios2();
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
	
	if (!inst) initAdios2();
}

void Adios2Reader::readPhaseSpaceHeader(Domain* domain, double timestep) {
//EMPTY
};


unsigned long Adios2Reader::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {

  auto variables = io->AvailableVariables();

  auto total_steps = std::stoi(variables["simulationtime"]["AvailableStepsCount"]);
  global_log->info() << "[Adios2Reader]: TOTAL STEPS " << total_steps << std::endl;

  if (_step == -1 ) {
      _step = total_steps - 1;
  }
  if (_step > total_steps) {
      global_log->error() << "[Adios2Reader]: Specified step is out of scope" << std::endl;
  } 
  if (_step < 0) {
      _step = total_steps + (_step + 1);
  }
  particle_count = std::stoi(variables["rx"]["Shape"]);
  global_log->info() << "[Adios2Reader]: Particle count: " << particle_count << std::endl;

  if (_mode == "rootOnly") {
      rootOnlyRead(particleContainer, domain, domainDecomp);
  } else if (_mode == "equalRanks") {
      // TODO: Implement
      equalRanksRead(particleContainer, domain, domainDecomp);
  } else {
      global_log->error() << "[Adios2Reader]: Unkown _mode '" << _mode << "'" << std::endl;
  }

    engine->Close();
    global_log->info() << "[Adios2Reader]: finish." << std::endl;

    return particle_count;
};

void Adios2Reader::rootOnlyRead(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {
	Timer inputTimer;
	inputTimer.start();

    std::vector<double> rx ,ry, rz, vx, vy, vz, qw, qx, qy, qz, Lx, Ly, Lz;
	std::vector<uint64_t> mol_id, comp_id;
    double simtime = 0.0;
    uint64_t buffer = 1000;
#ifdef ENABLE_MPI
    MPI_Datatype mpi_Particle;
	ParticleData::getMPIType(mpi_Particle);
#endif

    auto num_reads = particle_count / buffer;
    if (particle_count % buffer != 0) num_reads += 1;

    auto variables = io->AvailableVariables();

    for (int read = 0; read < num_reads; read++) {
        uint64_t offset = read * buffer;
        if (read == num_reads -1) buffer = particle_count % buffer;
        for (auto var : variables) {
            if(domainDecomp->getRank() == 0) {
                if (var.first == "rx") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, rx, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Positions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "ry") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, ry, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Positions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "rz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, rz, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Positions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "vx") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, vx, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Velocities should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "vy") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, vy, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Velocities should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "vz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, vz, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Velocities should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qw") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qw, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qx") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qx, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qy") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qy, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qz, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "Lx") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, Lx, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Angular momentum should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "Ly") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, Ly, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Angular momentum should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "Lz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, Lz, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Angular momentum should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "component_id") {
                    if (var.second["Type"] == "uint64_t") {
                        doTheRead(var.first, comp_id, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Component ids should be uint32_t (for now)." << std::endl;
                    }
                }
                if (var.first == "molecule_id") {
                    if (var.second["Type"] == "uint64_t") {
                        doTheRead(var.first, mol_id, buffer, offset);
                    } else {
                        global_log->error() << "[Adios2Reader] Molecule ids should be uint64_t (for now)." << std::endl;
                    }
                }
	            if (var.first == "simulationtime") {
	                if (var.second["Type"] == "double") {
	                    doTheRead(var.first, simtime);
	                } else {
	                    global_log->error() << "[Adios2Reader] Simulation time should be double (for now)." << std::endl;
	                }
	            }
			}
        }

        engine->PerformGets();

    	std::vector<Component> dcomponents;
    	if (_simulation.getEnsemble()->getComponents()->empty()) {
            auto attributes = io->AvailableAttributes();
    		auto comp_id_copy = comp_id;
			std::sort(comp_id_copy.begin(), comp_id_copy.end());
			auto u_iter = std::unique(comp_id_copy.begin(), comp_id_copy.end());
			comp_id_copy.resize(std::distance(comp_id_copy.begin(), u_iter));
			dcomponents.resize(comp_id_copy.size());
    		
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

    			dcomponents[comp].addLJcenter(centers[0], centers[0], centers[0], mass[0], eps[0], sigma[0]);
				dcomponents[comp].setName(cname[0]);
    		}
    	} else {
			dcomponents = *(_simulation.getEnsemble()->getComponents());
    	}


#ifdef ENABLE_MPI
        std::vector<ParticleData> particle_buff(buffer);
        if(domainDecomp->getRank() == 0) {
            for (int i = 0; i < buffer; i++) {

                    Molecule m1 = Molecule(mol_id[i], &dcomponents[comp_id[i]], rx[i], ry[i], rz[i], vx[i],
                                            vy[i], vz[i], qw[i], qx[i], qy[i], qz[i], Lx[i], Ly[i], Lz[i]);
                    ParticleData::MoleculeToParticleData(particle_buff[i], m1);
            }
        }
        MPI_Bcast(particle_buff.data(), buffer, mpi_Particle, 0, MPI_COMM_WORLD);
    	auto asd = 0;
        for (int j = 0; j < buffer; j++) {
			Molecule m = Molecule(particle_buff[j].id, &dcomponents[particle_buff[j].cid], particle_buff[j].r[0], particle_buff[j].r[1],
								particle_buff[j].r[2], particle_buff[j].v[0], particle_buff[j].v[1], particle_buff[j].v[2],
								particle_buff[j].q[0], particle_buff[j].q[1], particle_buff[j].q[2], particle_buff[j].q[3],
								particle_buff[j].D[0], particle_buff[j].D[1], particle_buff[j].D[2]);
            //ParticleData::ParticleDataToMolecule(particle_buff[j], m);
            // only add particle if it is inside of the own domain!
            if(particleContainer->isInBoundingBox(m.r_arr().data())) {
                particleContainer->addParticle(m, true, false);
            }

            dcomponents[m.componentid()].incNumMolecules();
            domain->setglobalRotDOF(
            dcomponents[m.componentid()].getRotationalDegreesOfFreedom()
                    + domain->getglobalRotDOF());

            // Only called inside GrandCanonical
            // TODO
            //global_simulation->getEnsemble()->storeSample(&m, componentid);
        }
#else
		for (int i = 0; i < particle_count; i++) {
			Molecule m = Molecule(mol_id[i], &dcomponents[comp_id[i]], rx[i], ry[i], rz[i], vx[i], vy[i], vz[i], qw[i],
								   qx[i], qy[i], qz[i], Lx[i], Ly[i], Lz[i]);

			// only add particle if it is inside of the own domain!
			if (particleContainer->isInBoundingBox(m.r_arr().data())) {
				particleContainer->addParticle(m, true, false);
			}

			dcomponents[m.componentid()].incNumMolecules();
			domain->setglobalRotDOF(dcomponents[m.componentid()].getRotationalDegreesOfFreedom() +
									domain->getglobalRotDOF());
		}
#endif
    	
        // Print status message
		unsigned long iph = num_reads / 100;
		if(iph != 0 && (read % iph) == 0)
			global_log->info() << "Finished reading molecules: " << read / iph
							   << "%\r" << std::flush;
    }

    global_log->info() << "Finished reading molecules: 100%" << std::endl;
	global_log->info() << "Reading Molecules done" << std::endl;

    inputTimer.stop();
	global_log->info() << "Initial IO took:                 "
					   << inputTimer.get_etime() << " sec" << std::endl;
    
	if(domain->getglobalRho() == 0.) {
		domain->setglobalRho(
				domain->getglobalNumMolecules() / domain->getGlobalVolume());
		global_log->info() << "Calculated Rho_global = "
						   << domain->getglobalRho() << endl;
	}

    _simulation.setSimulationTime(simtime);
    global_log->info() << "[Adios2Reader] simulation time is: " << simtime << std::endl;
};

void Adios2Reader::initAdios2() {
  auto& domainDecomp = _simulation.domainDecomposition();

  try {
    //get adios2 instance
#ifdef ENABLE_MPI
        inst = std::make_shared<adios2::ADIOS>((MPI_Comm) MPI_COMM_WORLD);
#else
		inst = std::make_shared<adios2::ADIOS>();
#endif
		io = std::make_shared<adios2::IO>(inst->DeclareIO("Input"));

        io->SetEngine(_adios2enginetype);

        if (!engine) {
            global_log->info() << "[Adios2Reader]: Opening File for reading." << _inputfile.c_str() << std::endl;
			engine = std::make_shared<adios2::Engine>(io->Open(_inputfile, adios2::Mode::Read));
        }
  }
    catch (std::invalid_argument& e) {
        global_log->fatal() 
                << "Invalid argument exception, STOPPING PROGRAM from rank" 
                << domainDecomp.getRank() 
                << ": " << e.what() << std::endl;
	mardyn_exit(1);
    } catch (std::ios_base::failure& e) {
        global_log->fatal() 
                << "IO System base failure exception, STOPPING PROGRAM from rank " 
                << domainDecomp.getRank() 
                << ": " << e.what() << std::endl;
	mardyn_exit(1);
    } catch (std::exception& e) {
        global_log->fatal() 
                << "Exception, STOPPING PROGRAM from rank" 
                << domainDecomp.getRank() 
                << ": " << e.what() << std::endl;
	mardyn_exit(1);
    }
    global_log->info() << "[Adios2Reader]: Init complete." << std::endl;
};

#endif // ENABLE_ADIOS2
