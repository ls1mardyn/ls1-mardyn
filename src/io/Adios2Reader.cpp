/*
 * Adios2Reader.cpp
 *
 *  Created on: April 2019
 *      Author: Tobias Rau
 */

#ifdef ENABLE_ADIOS2

#include "Adios2Reader.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#ifdef ENABLE_MPI
#include "parallel/ParticleData.h"
#include "parallel/DomainDecompBase.h"
#endif

using Log::global_log;

void Adios2Reader::init(ParticleContainer* particleContainer,
        DomainDecompBase* domainDecomp, Domain* domain) {
};

void Adios2Reader::readXML(XMLfileUnits& xmlconfig) {
  mode = "rootOnly";
  fname = "test.bp";
  xmlconfig.getNodeValue("filename", fname);
  xmlconfig.getNodeValue("mode", mode);
  // Adios2 step (frame)
  xmlconfig.getNodeValue("adios2Step", step);

  global_log->info() << "    [Adios2Reader]: readXML." << std::endl;

  if (!inst) initAdios2();
};

void Adios2Reader::readPhaseSpaceHeader(Domain* domain, double timestep) {
//EMPTY
};


unsigned long Adios2Reader::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {

  auto variables = io->AvailableVariables();

  auto total_steps = std::stoi(variables["simulationtime"]["AvailableStepsCount"]);
  global_log->info() << "[Adios2Reader]: TOTAL STEPS " << total_steps << std::endl;


  if (step > total_steps) {
      global_log->error() << "[Adios2Reader]: Specified step is out of scope" << std::endl;
  } 
  if (step < 0) {
      step = total_steps + (step + 1);
  }
  particle_count = std::stoi(variables["rx"]["Shape"]);
  global_log->info() << "    [Adios2Reader]: Particle count: " << particle_count << std::endl;

  if (mode == "rootOnly") {
      rootOnlyRead(particleContainer, domain, domainDecomp);
  } else if (mode == "equalRanks") {
      // TODO: Implement
      equalRanksRead(particleContainer, domain, domainDecomp);
  } else {
      global_log->error() << "[Adios2Reader]: Unkown Mode '" << mode << "'" << std::endl;
  }

    engine->Close();
    global_log->info() << "    [Adios2Reader]: finish." << std::endl;

    return particle_count;
};

void Adios2Reader::rootOnlyRead(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {

	Timer inputTimer;
	inputTimer.start();

    std::vector<double> rx ,ry, rz, vx, vy, vz, qw, qx, qy, qz, Lx, Ly, Lz;
	std::vector<uint64_t> mol_id, comp_id;
    double simtime;
    uint64_t buffer = 1000;
    MPI_Datatype mpi_Particle;
	ParticleData::getMPIType(mpi_Particle);


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
                        global_log->error() << "    [Adios2Reader] Positions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "ry") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, ry, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Positions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "rz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, rz, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Positions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "vx") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, vx, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Velocities should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "vy") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, vy, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Velocities should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "vz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, vz, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Velocities should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qw") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qw, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qx") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qx, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qy") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qy, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "qz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, qz, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Quaternions should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "Lx") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, Lx, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Angular momentum should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "Ly") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, Ly, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Angular momentum should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "Lz") {
                    if (var.second["Type"] == "double") {
                        doTheRead(var.first, Lz, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Angular momentum should be doubles (for now)." << std::endl;
                    }
                }
                if (var.first == "component_id") {
                    if (var.second["Type"] == "uint64_t") {
                        doTheRead(var.first, comp_id, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Component ids should be uint32_t (for now)." << std::endl;
                    }
                }
                if (var.first == "molecule_id") {
                    if (var.second["Type"] == "uint64_t") {
                        doTheRead(var.first, mol_id, buffer, offset);
                    } else {
                        global_log->error() << "    [Adios2Reader] Molecule ids should be uint64_t (for now)." << std::endl;
                    }
                }
            }
            if (var.first == "simulationtime") {
                if (var.second["Type"] == "double") {
                    doTheRead(var.first, simtime);
                } else {
                    global_log->error() << "    [Adios2Reader] Simulation time should be double (for now)." << std::endl;
                }
            }
        }
        engine->PerformGets();

        std::vector<ParticleData> particle_buff(buffer);
        auto& dcomponents = *(_simulation.getEnsemble()->getComponents());
        if(domainDecomp->getRank() == 0) {
            for (int i = 0; i < buffer; i++) {
                    Molecule m1 = Molecule(mol_id[i], &dcomponents[comp_id[i]], rx[i], ry[i], rz[i], vx[i],
                                            vy[i], vz[i], qw[i], qx[i], qy[i], qz[i], Lx[i], Ly[i], Lz[i]);
                    ParticleData::MoleculeToParticleData(particle_buff[i], m1);
            }
        }
        MPI_Bcast(particle_buff.data(), buffer, mpi_Particle, 0, MPI_COMM_WORLD);
        for (int j = 0; j < buffer; j++) {
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

            // Only called inside GrandCanonical
            // TODO
            //global_simulation->getEnsemble()->storeSample(&m, componentid);
        }
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
    global_log->info() << "    [Adios2Reader] simulation time is: " << simtime << std::endl;

};

void Adios2Reader::initAdios2() {
  auto& domainDecomp = _simulation.domainDecomposition();

  try {
    //get adios2 instance
        inst = std::make_shared<adios2::ADIOS>((MPI_Comm) MPI_COMM_WORLD);

        // declare io as output
        io = std::make_shared<adios2::IO>(inst->DeclareIO("Input"));

        // use bp engine
        io->SetEngine("BPFile");

        if (!engine) {
            global_log->info() << "    [Adios2Reader]: Opening File for writing." << fname.c_str() << std::endl;
            engine = std::make_shared<adios2::Engine>(io->Open(fname, adios2::Mode::Read));
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
    global_log->info() << "    [Adios2Reader]: Init complete." << std::endl;
};

#endif // ENABLE_ADIOS2
