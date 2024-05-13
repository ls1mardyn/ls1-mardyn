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
 *
 * \b NOTE:
 *  - You must ensure that the files under /sys/class/powercap/intel-rapl/ are readable
 *  - The total energy is measured over the package and DRAM domains, if present.
 */
class EnergyRAPL : public PluginBase {
private:
	/**
	 * Class to measure the energy consumption using the sysfs interface for the corresponding RAPL counter domain
	 */
	class RAPLCounter {
	private:
		std::string _microJoulesPath;
		long long _lastMicroJoules;
		long long _rangeMicroJoules;
		long long _microJoules;

	public:
		/**
		 * Creates a new instance of a counter for the energy consumption of a specific domain
		 *
		 * @param domainBasePath The sysfs base path of the corresponding RAPL domain (must contain the files: name,
		 * max_energy_range_uj, energy_uj)
		 */
		RAPLCounter(const std::string domainBasePath);

		void reset();

		double update();
	};

	const char* const _basePathRAPL = "/sys/class/powercap/intel-rapl/";
	std::vector<RAPLCounter> _counters;
	int _writeFrequency = 0;
	std::string _outputprefix;
	double _joules;
	unsigned long _simstep;
	std::chrono::_V2::steady_clock::time_point _simstart;

#ifdef ENABLE_MPI
	char _processorName[MPI_MAX_PROCESSOR_NAME];
	int _processorNameLength;
	int _thisRank;
#endif

	/**
	 * Returns the number of CPU sockets.
	 */
	int getNumberOfPackages();

	/**
	 * Compute and output the energy consumed over all RAPLCounter instances on all nodes
	 * (Only outputs on the root rank using MPI)
	 */
	void outputEnergyJoules();

public:
	EnergyRAPL() = default;

	~EnergyRAPL() override = default;

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	std::string getPluginName() override { return std::string("EnergyRAPL"); }

	static PluginBase* createInstance() { return new EnergyRAPL(); }

	/** @brief Read in XML configuration for DecompWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <plugin name="EnergyRAPL">
		 <!-- Outputs total energy consumption only once at the end if writefrequency is zero (default) -->
		 <writefrequency>INTEGER</writefrequency>
		 <!-- Uses info logger instead of a .tsv file if outputprefix is not given -->
		 <outputprefix>STRING</outputprefix>
	   </plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;
};

#endif  // SRC_PLUGIN_ENERGYRAPL_H_