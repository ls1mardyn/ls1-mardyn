#ifndef SRC_IO_LOADBALANCEWRITER_H_
#define SRC_IO_LOADBALANCEWRITER_H_

#include "plugins/PluginBase.h"
#include "utils/Timer.h"

#include <string>
#include <vector>
#include <deque>

/** name of the LoadbalanceWriter's default timer */
#define LB_WRITER_DEFAULT_TIMER_NAME "LoadbalanceWriter_default"

/** @brief The LoadbalanceWriter writes out information about the programs load balance.
 *
 * The LoadbalanceWriter writes out load balance information for MPI processes based on
 * timers. The 'default' timer used measures the time between invocations of this plugin.
 * Additional timers used in the program and registered in the Simulation can be added
 * to compute the load balance of specific code parts.
 *
 * The output file will include for each time step in the simulation one line.
 * In each line the following values are stored for each LB timer:
 * - min, max times over all processes
 * - the factor max/min
 * - the imbalance 1 - average / max. This describes lost performance due to imbalances for the current time step.
 * - an averaged imbalance. Compared to the step-wise imbalance, the time-values are averaged for averageLength time
 * steps. In the first few steps, this averaging is not fully possible, thus only the average of the first steps is
 * taken. The averaged imbalance aims to reduce imbalances caused by noise and performance fluctuations and is thus a
 * more accurate measure for the quality of the partitioning. It is (in average) lower than the step-wise imbalance. In
 * addition, when comparing the step-wise imbalance and the averaged imbalance, information about the fluctuations can
 * be retrieved.
 *
 * @note Warning thresholds (level) (for the factor max / min) for each timer can be set (see documentation of
 * readXML()) which will output a warning to the logfile if the factor max_time / min_time is above the threshold. This
 * warning level should be bigger than 1.!
 *
 * @todo This plugin may be extended to threads
 */
class LoadbalanceWriter : public PluginBase {
public:
	LoadbalanceWriter();
	~LoadbalanceWriter() override = default;

	/** @brief Read in XML configuration for LoadbalanceWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * parameters:
	 * 	name: name of the timer
	 * 	warninglevel: warnings will be printed if "max time / min time > warninglevel"
	 * 	incremental: specifies whether the timer is incremental or not,
	 * 		e.g., a timer just measuring the time for the current time step is not incremental,
	 * 			but one measuring the time since the first time step is incremental.
	 * \code{.xml}
	   <outputplugin name="LoadbalanceWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <averageLength>INTEGER</averageLength>
	     <outputfilename>STRING</outputfilename>
	     <timers> <!-- additional timers -->
	        <timer> <name>LoadbalanceWriter_default</name> <warninglevel>DOUBLE</warninglevel> </timer>
	        <timer> <name>STRING</name> <warninglevel>DOUBLE</warninglevel> <incremental>BOOL</incremental> </timer>
	        <!-- ... -->
	     </timers>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain) override;

	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    ) override;

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override {
		/* nothing to do */
	}

	std::string getPluginName() override {
		return std::string("LoadbalanceWriter");
	}
	static PluginBase* createInstance() { return new LoadbalanceWriter(); }

private:
	void recordTimes(long unsigned int simstep);
	void resetTimes();
	void writeOutputFileHeader();
	void writeLBEntry(size_t id, std::ofstream &outputfile, int numRanks);
	void flush(DomainDecompBase* domainDecomp);
	void displayWarning(unsigned long simstep, const std::string& timername, double f_LB);

	unsigned long _writeFrequency;
	unsigned long _averageLength{10};
	std::string _outputFilename;
	Timer *_defaultTimer;
	std::vector<std::string> _timerNames;
	std::vector<double> _times;

	// Holds multiple time values (not flushed after write) to get averaged load values.
	std::deque<double> _timesForAverageBuffer;
	unsigned long _lastTimesOldValueCount{0};  // How many old values _lastTimes holds from before a flush.
	std::vector<double> _global_sum_average_times;
	std::vector<double> _global_max_average_times;

	std::vector<double> _sum_times;
	std::vector<double> _global_sum_times;
	std::vector<double> _global_times;
	std::vector<unsigned long> _simsteps;
	std::map<std::string, double> _warninglevels;
	std::map<std::string, bool> _incremental;  // describes whether the timer will continuously increase and the difference between two calls should be used as timer value
	std::map<std::string, double> _incremental_previous_times;  // previous times of the incremental timers
	std::vector<double> getAveragedTimes();
};

#endif  // SRC_IO_LOADBALANCEWRITER_H_

