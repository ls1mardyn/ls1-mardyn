/*
 * HardwareInfo.h
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
 *
 * The total number of threads isn't unique, but is used to store the output of mardyn_get_num_threads, which may not
 * always match mardyn_get_max_threads, such as in heterogeneous hardware scenarios. totalThreads is printed to stdout.
 * When writing to file, the value from _threadData[0] is used.
 */
struct ThreadwiseInfo {
	int thread;
	int cpuID, numa;
	ThreadwiseInfo(int threadNum = 0, int totalThreadsNum = 1, unsigned int openMPCPUID = 0,
				   unsigned int openMPNUMA = 0)
		: thread(threadNum), cpuID(openMPCPUID), numa(openMPNUMA) {}
};

/**
 * @brief Prints the information of the hardware of current simulation to stdout or to file.
 * @author Amartya Das Sharma
 *
 * The plugin either prints to stdout, or a json file. Default behaviour is to print to stdout.
 * Printed information includes rank number, thread number, NUMA domain and name of node on which the rank resides.
 * Thread stats depend on sched.h, which requires glibc and has less functionality on older versions. In serial mode
 * (non MPI) node name is gotten from unistd.h, which is absent on non POSIX compliant OSs.
 */
class HardwareInfo : public PluginBase {
public:
	HardwareInfo() = default;
	~HardwareInfo() override = default;

	/** @brief Read in XML configuration for HardwareInfo.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <plugin name="HardwareInfo">
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

	std::string getPluginName() override { return "HardwareInfo"; }

	static PluginBase* createInstance() { return new HardwareInfo(); }

private:
	/**
	 * @brief Fills the data members _threadData, _rank, _totalRanks, and _nodeName, and sets _dataPopulated to
	 * true
	 *
	 * This function is tasked with performing all of the necessary quesries to get the hardware information. It uses
	 * MPI functions to fill _rank and _totalRanks, openMP fucntions to populate _threadData. Depending on the
	 * availability of functions in sched.h (glibc version dependent) and on the availability of unistd.h (POSIX OS
	 * required), not all information may be available (NUMA domain and CPU ID). This data is filled with -1 instead. At
	 * the end of this function, an openMP barrier is present. Also, the value _dataPopulated is set to true for sanity
	 * checks.
	 *
	 * @param domainDecomp Needed to access the communicator to get MPI rank and world size
	 */
	void populateData(DomainDecompBase* domainDecomp);

	/**
	 * @brief Pretty-prints _threadData, _rank, _totalRanks, and _nodeName using the logger.
	 *
	 * Any data that is -1 is assumed to be incomplete and skipped.
	 */
	void printDataToStdout();

	/**
	 * @brief Writes _threadData, _rank, _totalRanks, and _nodeName to a specified json file.
	 *
	 * No json libraries are used, so all the data is processed as raw strings and directly written to file. The file is
	 * opened as an MPI file and all ranks write directly to it in specified offsets. Thus, the data will be ordered by
	 * ranks in increasing order.
	 * The root node writes some header information, and closes the ending braces. In case the simulation is compiled in
	 * serial mode, this function also writes the data serially.
	 * The json structure is as follows:
	 * \code{.json}
	   {
		   "total_ranks": _totalRanks,
		   "rank_data": {
			   <data from convertFullDataToJson()>
		   }
	   }
	   \endcode
	 *
	 * @param domainDecomp Needed to access the communicator to set up the MPI file to be written to
	 */
	void writeDataToFile(DomainDecompBase* domainDecomp);

	/**
	 * @brief Wrapper for logic that converts the entire _threadData into a json string.
	 *
	 * This contains all the data that is rank dependent, and assumes that all data that is -1 is absent. Since the data
	 * will be written in order of ranks, the trailing comma for the data of the final rank is removed. Thus, without
	 * the trailing comma, this function returns a fully parseable json, without any closing braces required. Each line
	 * of the json is prepended with two tabs.
	 * The json structure is as follows:
	 * \code{.json}
	 * _rank: {
		   "node_name": _nodeName,
		   "total_threads": omp_get_max_threads(),
		   "thread_data": {
			   _threadID: {"cpu_ID": INT, "numa_domain": INT}
		   }
	   }
	 * \endcode
	 * @return const std::string The rank-dependent hardware information formatted in json.
	 */
	std::string convertFullDataToJson() const;

	/**
	 * @brief Stores the filename read in from XML, appended with a ".json".
	 */
	std::string _filename = "";

	/**
	 * @brief Stores the local rank, and number of ranks in the communicator.
	 */
	int _rank, _totalRanks;

	int _maxThreads;
	/**
	 * @brief Stores name of the current node.
	 */
	std::string _nodeName;

	std::string _maxRam;

	std::string _cpuInfo;

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
