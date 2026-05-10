/*
 * PinningInfo.h
 *
 *  Created on: 16 Dec 2025
 *      Author: amartyads
 */

#pragma once

#include <vector>

#include "PluginBase.h"

/**
 * @brief Stores information unique to each openMP thread. This includes thread number, cpu ID of the thread, and NUMA
 * domain of the thread.
 */
struct ThreadwiseInfo {
	int thread;
	int cpuID, numa;
	ThreadwiseInfo(int threadNum = 0, int openMPCPUID = 0, int openMPNUMA = 0)
		: thread(threadNum), cpuID(openMPCPUID), numa(openMPNUMA) {}
};

/**
 * @brief Prints the pinning information of the hardware of current simulation to stdout or to file.
 * @author Amartya Das Sharma
 *
 * The plugin either prints to stdout, or a JSON file. Default behaviour is to print to stdout.
 * Printed information includes rank number, thread number, NUMA domain and name of node on which the rank resides.
 * Thread stats depend on sched.h, which requires glibc and has less functionality on older versions.
 */
class PinningInfo : public PluginBase {
public:
	PinningInfo() = default;
	~PinningInfo() override = default;

	/** @brief Read in XML configuration for PinningInfo.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <plugin name="PinningInfo">
		 <!-- Uses info logger instead of a .json file if filename is not given -->
		 <filename>STRING</filename>
	   </plugin>
	   \endcode
	 * The function also appends ".json" to the filename.
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	/**
	 * @brief Calls the functions to populate hardware data, and write to log or file accordingly
	 *
	 * @param domainDecomp Needed to access the communicator, to pass to other functions
	 */
	void init(ParticleContainer*, DomainDecompBase* domainDecomp, Domain*) override;

	void endStep(ParticleContainer*, DomainDecompBase*, Domain*, unsigned long) override {};
	void finish(ParticleContainer*, DomainDecompBase*, Domain*) override {};

	std::string getPluginName() override { return "PinningInfo"; }

	static PluginBase* createInstance() { return new PinningInfo(); }

private:
	/**
	 * @brief Fills the data members _threadData, _rank, _totalRanks, _nodeName, _maxThreads, _cpuInfo, _cpuArch and
	 * _maxRam. Sets _dataPopulated to true.
	 *
	 * This function is tasked with performing all of the necessary queries to get the pinning information. It uses
	 * MPI functions to fill _rank and _totalRanks, and openMP functions to populate _threadData. Depending on the
	 * availability of functions in sched.h (glibc version dependent) not all information may be available (NUMA domain
	 * and CPU ID). This data is filled with -1 instead. RAM info is populated by calling sysinfo from sys/sysinfo.h,
	 * and CPU info is obtained by parsing /proc/cpuinfo. Also, the value _dataPopulated is set to true for sanity
	 * checks.
	 *
	 * @param domainDecomp Needed to access the communicator to get MPI rank and world size
	 */
	void populateData(DomainDecompBase* domainDecomp);

	/**
	 * @brief Pretty-prints _threadData, _rank, _totalRanks, _nodeName, _maxThreads, _cpuInfo, _cpuArch and _maxRam
	 * using the logger.
	 *
	 * Any data that is -1 is assumed to be incomplete.
	 */
	void printDataToStdout();

	/**
	 * @brief Writes _threadData, _rank, _totalRanks, _nodeName, _maxThreads, _cpuInfo, _cpuArch and _maxRam to a
	 * specified JSON file.
	 *
	 * The nlohmann/json library is used to create and autoindent the JSON before dumping it to file. Each rank
	 * serialises its data and sends it to the root rank, which compiles it into one JSON object. A tab width of 4
	 * spaces is used. Missing values are denoted with -1.
	 * The JSON structure is as follows:
	 * \code{.json}
	   {
		   "total_ranks": _totalRanks,
		   "rank_data": {
			   _rank: {
				   "node_name": _nodeName,
				   "cpu_arch": _cpuArch,
				   "cpu_info": _cpuInfo,
				   "max_avail_ram": _maxRam,
				   "total_threads": _maxThreads,
				   "thread_data": {
				   _threadID: {"cpu_ID": INT, "numa_domain": INT}
			   }
		   }
	   }
	   \endcode
	 * Since JSON is unordered, the elements may be in any order.
	 *
	 * @param domainDecomp Needed to access the communicator to receive the JSON data from participating ranks
	 */
	void writeDataToFile(DomainDecompBase* domainDecomp);

	/**
	 * @brief Stores the filename read in from XML, appended with a ".json".
	 */
	std::string _filename = "";

	/**
	 * @brief Stores the local rank, and number of ranks in the communicator.
	 */
	int _rank, _totalRanks;

	/**
	 * @brief Stores the maximum number of OpenMP threads possible.
	 */
	int _maxThreads;
	/**
	 * @brief Stores name of the current node.
	 */
	std::string _nodeName;

	/**
	 * @brief Stores a string with the max RAM of the current node, in GB and GiB.
	 */
	std::string _maxRam;

	/**
	 * @brief Stores a string containing information needed to identify the CPU.
	 */
	std::string _cpuInfo;

	/**
	 * @brief Stores the CPU architecture.
	 */
	std::string _cpuArch;

	/**
	 * @brief Stores all thread data, in the form of ThreadwiseInfo struct objects.
	 */
	std::vector<ThreadwiseInfo> _threadData;

	/**
	 * @brief Boolean for sanity checks before printing.
	 */
	bool _dataPopulated = false;
};
