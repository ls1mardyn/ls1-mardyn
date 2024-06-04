#ifndef SRC_PLUGIN_ENERGYRAPL_H_
#define SRC_PLUGIN_ENERGYRAPL_H_

#include <chrono>
#include <string>

#include "plugins/PluginBase.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

/**
 * @brief This plugin computes the energy consumption of the simulation using RAPL counters via the sysfs interface.
 * @author Ruben Horn
 *
 * The output of this plugin is either written to the info logger or to a tab separated file (.tsv) with the columns
 * milliseconds, simstep and joules every \c writefrequency simulation steps or once at the end.
 *
 * \b NOTE:
 *  - You must ensure that the files under /sys/class/powercap/intel-rapl/ are readable.
 *  - The total energy is measured over the package and DRAM domains, if present.
 *    (See section 15.10.2 of the [Intel Architectures Software Developer’s Manuals](https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html) Volume 3.
 *    The actual domains that are used are printed to the info log.)
 *    Energy measurements using RAPL are \b not representative of the total energy consumption of the system and comparing different systems is not trivially possible.
 *    Instead, they should only be used to compare the energy consumption of different applications/settings on the same system.
 *    For more practical information on RAPL, see [Khan et al. 2018](https://doi.org/10.1145/3177754). RAPL is an Intel technology and has limited support on AMD (cf. [Schöne et al. 2021](https://arxiv.org/abs/2108.00808)).
 */
class EnergyRAPL : public PluginBase {
private:
	/**
	 * @brief Class to measure the energy consumption using the sysfs interface for the corresponding RAPL counter
	 * domain
	 */
	class RAPLCounter {
	private:
		/**
		 * @brief Path to the file containing the current value for this domain
		 */
		std::string _microJoulesPath;
		/**
		 * @brief Last value measured
		 */
		long long _lastMicroJoules;
		/**
		 * @brief Maximum counter value (used for overflow detection)
		 */
		long long _rangeMicroJoules;
		/**
		 * @brief Aggregated value since initialization or reset
		 */
		long long _microJoules;

	public:
		/**
		 * @brief Creates a new instance of a counter for the energy consumption of a specific domain
		 *
		 * @param domainBasePath The sysfs base path of the corresponding RAPL domain (must contain the files: name,
		 * max_energy_range_uj, energy_uj)
		 */
		RAPLCounter(const std::string &domainBasePath);

		/**
		 * @brief Reset the counter for this domain to 0
		 */
		void reset();

		/**
		 * @brief Takes a new measurement for this domain
		 * @return The aggregated value since initialization or last reset
		 */
		double update();
	};

	/**
	 * @brief Base path for the RAPL sysfs interface
	 */
	const char* const _basePathRAPL = "/sys/class/powercap/intel-rapl/";

	/**
	 * @brief Counter objects for all RAPL domains used to compute the energy consumption
	 */
	std::vector<RAPLCounter> _counters;

	/**
	 * @brief Number of simulation steps between outputs (only output at the end if 0)
	 */
	int _writeFrequency = 0;
	/**
	 * @brief Filename without extension for outputs (output to stdout if empty string)
	 */
	std::string _outputprefix;
	/**
	 * @brief Aggregated energy consumption over all \ref _counters
	 */
	double _joules;
	/**
	 * @brief Simulation step (updated by plugin hook)
	 */
	unsigned long _simstep;
	/**
	 * @brief Start time of the simulation (used to compute timestamp in milliseconds)
	 */
	std::chrono::time_point<std::chrono::steady_clock> _simstart;

#ifdef ENABLE_MPI
	/**
	 * @brief Unique name of the node associated with the MPI rank
	 */
	char _processorName[MPI_MAX_PROCESSOR_NAME];
	/**
	 * @brief The MPI rank
	 */
	int _thisRank;
#endif

	/**
	 * @brief Returns the number of CPU sockets.
	 */
	int getNumberOfPackages();

	/**
	 * @brief Compute and output the energy consumed over all RAPLCounter instances on all nodes
	 * (Only outputs on the root rank using MPI)
	 */
	void outputEnergyJoules();

public:
	/**
	 * @brief Creates an \b uninitialized instance of the plugin
	 */
	EnergyRAPL() = default;

	/**
	 * @brief Cleans up an instance of the plugin
	 */
	~EnergyRAPL() override = default;

	/**
	 * @brief Initializes the plugin
	 *
	 * - Initializes MPI/output related variables of plugin object
	 * - Scans RAPL domains and initializes RAPL counter objects and resets them
	 */
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	/**
	 * @brief Updates RAPL counter objects and triggers output depending on plugin settings
	 */
	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;

	/**
	 * @brief Triggers output depending on plugin settings
	 */
	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	/**
	 * @brief String identifying the plugin \ref EnergyRAPL
	 */
	std::string getPluginName() override { return std::string("EnergyRAPL"); }

	/**
	 * @brief Creates an \b uninitialized instance of the plugin
	 */
	static PluginBase* createInstance() { return new EnergyRAPL(); }

	/** @brief Read in XML configuration for DecompWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <plugin name="EnergyRAPL">
		 <!-- Outputs total energy consumption only once at the end if writefrequency is 0 (default) -->
		 <writefrequency>INTEGER</writefrequency>
		 <!-- Uses info logger instead of a .tsv file if outputprefix is not given -->
		 <outputprefix>STRING</outputprefix>
	   </plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;
};

#endif  // SRC_PLUGIN_ENERGYRAPL_H_
