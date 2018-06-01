#ifndef SRC_IO_LOADBALANCEWRITER_H_
#define SRC_IO_LOADBALANCEWRITER_H_

#include <string>

#include "plugins/PluginBase.h"
#include "utils/Timer.h"

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
 *
 * Warning levels for each timer can be set which will output a waring to the logfile.
 *
 * @todo This plugin may be extended to threads
 */
class LoadbalanceWriter : public PluginBase {
public:
	LoadbalanceWriter();
	~LoadbalanceWriter() {}

	/** @brief Read in XML configuration for DecompWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="LoadbalanceWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputfilename>STRING</outputfilename>
	     <timers> <!-- additional timers -->
	        <timer> <name>LoadbalanceWriter_default</name> <warninglevel>DOUBLE</warninglevel> </timer>
	        <timer> <name>STRING</name> <warninglevel>DOUBLE</warninglevel> </timer>
	        <!-- ... -->
	     </timers>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);

	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    );

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
		/* nothing to do */
	}

	std::string getPluginName() {
		return std::string("LoadbalanceWriter");
	}
	static PluginBase* createInstance() { return new LoadbalanceWriter(); }

private:
	void recordTimes(long unsigned int simstep);
	void resetTimes();
	void writeOutputFileHeader();
	void writeLBEntry(size_t id, std::ofstream &outputfile);
	void flush(DomainDecompBase* domainDecomp);
	void displayWarning(unsigned long simstep, std::string timername, double f_LB);

	unsigned long _writeFrequency;
	std::string _outputFilename;
	Timer *_defaultTimer;
	std::vector<std::string> _timerNames;
	std::vector<double> _times;
	std::vector<double> _global_times;
	std::vector<unsigned long> _simsteps;
	std::map<std::string, double> _warninglevels;
};

#endif  // SRC_IO_LOADBALANCEWRITER_H_

