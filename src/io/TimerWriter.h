//
// Created by seckler on 10.10.19.
//
#pragma once

#include <fstream>
#include <vector>

#include "plugins/PluginBase.h"
#include "utils/Timer.h"

/**
 * Output plugin to write the timing info of timers to files.
 * Every rank writes the info to a separate file.
 * The average time spent in one simulation step is printed for the provided timers.
 * Hereby the average is taken for writefrequency steps (printed if (sim step % write frequency == 0), except for the
 * zeroth time step).
 */
class TimerWriter : public PluginBase {
public:
	/** @brief Read in XML configuration for TimerWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * parameters:
	 * 	name: name of the timer
	 * 	incremental: specifies whether the timer is incremental or not, (default is false)
	 * 		e.g., a timer just measuring the time for the current time step is not incremental,
	 * 			but one measuring the time since the first time step is incremental.
	 * \code{.xml}
	   <outputplugin name="TimerWriter">
		 <writefrequency>INTEGER</writefrequency>
		 <outputprefix>STRING</outputprefix><!--default is mardyn-timers. -rankX_TIMESTAMP.dat will be appended. -->
		 <timers> <!-- timers that should be printed -->
			<timer> <name>STRING</name> <incremental>BOOL</incremental> </timer>
			<!-- ... -->
		 </timers>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits &xmlconfig) override;

	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override {
		/* nothing to do */
	}

	std::string getPluginName() override { return std::string("TimerWriter"); }
	static PluginBase *createInstance() { return new TimerWriter(); }

private:
	unsigned long _writeFrequency{0ul};
	std::ofstream _fileStream{};
	std::string _outputPrefix{""};
	std::vector<std::string> _timerNames{};
	std::vector<bool> _incremental{};
	std::vector<double> _times{};
	unsigned long _stepsSinceLastWrite{0ul};
};
