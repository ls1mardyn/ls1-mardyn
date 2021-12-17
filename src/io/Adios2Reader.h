#ifndef ADIOS2_READER_H_
#define ADIOS2_READER_H_
/*
 * \file Adios2Reader.h
 *
 * Allows to read ADIOS2 phase space series files as checkpoints.
 *
 */
#ifdef ENABLE_ADIOS2
#include "io/InputBase.h"
#include "molecules/Component.h"
#include "molecules/FullMolecule.h"
#include "molecules/MoleculeForwardDeclaration.h"
#include "utils/Logger.h"

#ifdef MARDYN_AUTOPAS
#include "molecules/AutoPasSimpleMolecule.h"
#else
#include "molecules/FullMolecule.h"
#endif

#include <errno.h>
#include <chrono>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "adios2.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

class Adios2Reader : public InputBase {
public:
	Adios2Reader() = default;
	~Adios2Reader() override = default;

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);

	/** @brief Read in XML configuration for Adios2Reader.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
		<phasespacepoint>
			<file type="adios2">
				<filename>STRING</filename>
				<adios2enginetype>STRING<!-- For possible engines see the ADIOS2 doc (default: BP4)
	 --></adios2enginetype> <adios2step>INTEGER<!-- Step in the ADIOS2 file --></adios2step> <mode>STRING<!-- Possible
	 options are rootOnly (only root node reads the data and broadcasts to all ranks) and parallelRead (each rank reads
	 a portion of the data and distributes it so every rank holds all data) --></mode>
			</file>
		</phasespacepoint>
	   \endcode
	*/
	void readXML(XMLfileUnits& xmlconfig) override;

	void testInit(std::string infile = "mardyn.bp", int step = 0, std::string adios2enginetype = "BP4",
				  std::string mode = "rootOnly");

	void readPhaseSpaceHeader(Domain* domain, double timestep) override;

	std::string getPluginName() { return std::string("Adios2Reader"); }

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) override;

private:
	void initAdios2();
	void rootOnlyRead(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);
	void parallelRead(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);

	// input filename, from XML
	std::string _inputfile;

	std::string _adios2enginetype;
	std::string _mode;
	int _step;
	uint64_t particle_count;
	std::vector<Component> _dcomponents;
	double _simtime = 0.0;
	bool _single_precission = false;
	// main instance
	std::shared_ptr<adios2::ADIOS> mainInstance;
	std::shared_ptr<adios2::Engine> engine;
	std::shared_ptr<adios2::IO> io;

	template <typename T>
	void doTheRead(std::string var_name, std::vector<T>& container) {
		uint64_t num = 1;
		auto advar = io->InquireVariable<T>(var_name);
		advar.SetStepSelection({_step, 1});
		auto info = engine->BlocksInfo(advar, 0);
		auto shape = info[0].Count;
		std::for_each(shape.begin(), shape.end(), [&](decltype(num) n) { num *= n; });
		container.resize(num);
		engine->Get<T>(advar, container);
	}

	/** @brief Read a variable from the provided ADIOS2 file using the provided offset and bufferSize.
	 *   The offset and bufferSize point to a section of the data.
	 *
	 *   @tparam T double or float depending on what is provided in the file
	 *   @param var_name Name of the variable to read
	 *   @param container Target for the actual data
	 *   @param bufferSize Length of the current read
	 *   @param offset offset to the current read
	 *
	 */
	template <typename T>
	void doTheRead(const std::string var_name, std::vector<T>& container, const uint64_t bufferSize, const uint64_t offset) {
		adios2::Dims readsize({bufferSize});
		adios2::Dims offs({offset});

		Log::global_log->debug() << "[Adios2Reader]: Var name " << var_name << std::endl;
		auto advar = io->InquireVariable<T>(var_name);
		advar.SetStepSelection({_step, 1});
		Log::global_log->debug() << "[Adios2Reader]: buffer " << bufferSize << " offset " << offset << std::endl;
		for (auto entry : advar.Shape()) Log::global_log->debug() << "[Adios2Reader]: shape " << entry << std::endl;
		advar.SetSelection(adios2::Box<adios2::Dims>(offs, readsize));
		container.resize(bufferSize);
		engine->Get<T>(advar, container);
	}

	template <typename T>
	void doTheRead(std::string var_name, T& value) {
		auto advar = io->InquireVariable<T>(var_name);
		advar.SetStepSelection({_step, 1});
		engine->Get<T>(advar, value);
	}

	template <typename T>
	Molecule fillMolecule(int i_, std::vector<uint64_t>& mol_id_, std::vector<uint64_t>& comp_id_, T& rx_, T& ry_,
						  T& rz_, T& vx_, T& vy_, T& vz_, T& qw_, T& qx_, T& qy_, T& qz_, T& Lx_, T& Ly_, T& Lz_) {
		Molecule m1;
		if (qw_.empty()) {
			m1 = Molecule(mol_id_[i_], &_dcomponents[comp_id_[i_]], rx_[i_], ry_[i_], rz_[i_], vx_[i_], vy_[i_],
						  vz_[i_], 1, 0, 0, 0, 0, 0, 0);
		} else {
			m1 = Molecule(mol_id_[i_], &_dcomponents[comp_id_[i_]], rx_[i_], ry_[i_], rz_[i_], vx_[i_], vy_[i_],
						  vz_[i_], qw_[i_], qx_[i_], qy_[i_], qz_[i_], Lx_[i_], Ly_[i_], Lz_[i_]);
		}
		return m1;
	}

	template <typename T>
	void performInquire(std::map<std::string, adios2::Params> variables, uint64_t buffer, uint64_t offset,
						std::vector<T>& rx, std::vector<T>& ry, std::vector<T>& rz, std::vector<T>& vx,
						std::vector<T>& vy, std::vector<T>& vz, std::vector<T>& qw, std::vector<T>& qx,
						std::vector<T>& qy, std::vector<T>& qz, std::vector<T>& Lx, std::vector<T>& Ly,
						std::vector<T>& Lz, std::vector<uint64_t>& mol_id, std::vector<uint64_t>& comp_id) {
		for (auto var : variables) {
			if (var.first == "rx") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, rx, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Positions should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "ry") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, ry, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Positions should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "rz") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, rz, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Positions should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "vx") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, vx, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Velocities should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "vy") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, vy, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Velocities should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "vz") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, vz, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Velocities should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "qw") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, qw, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Quaternions should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "qx") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, qx, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Quaternions should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "qy") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, qy, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Quaternions should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "qz") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, qz, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Quaternions should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "Lx") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, Lx, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Angular momentum should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "Ly") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, Ly, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Angular momentum should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "Lz") {
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, Lz, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Angular momentum should be doubles or floats (for now)."
										<< std::endl;
				}
			}
			if (var.first == "component_id") {
				if (var.second["Type"] == "uint64_t") {
					doTheRead(var.first, comp_id, buffer, offset);
				} else {
					global_log->error() << "[Adios2Reader] Component ids should be uint64_t (for now)." << std::endl;
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
				if (var.second["Type"] == "double" || var.second["Type"] == "float") {
					doTheRead(var.first, _simtime);
				} else {
					global_log->error() << "[Adios2Reader] Simulation time should be double or float (for now)."
										<< std::endl;
				}
			}
		}
	}
};
#endif  // ENABLE_ADIOS2
#endif /* ADIOS2_READER_H_*/