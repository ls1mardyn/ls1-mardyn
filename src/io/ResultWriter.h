#ifndef SRC_IO_RESULTWRITER_H_
#define SRC_IO_RESULTWRITER_H_

#include <fstream>
#include <string>

#include "plugins/PluginBase.h"
#include "utils/Accumulator.h"


/** @brief Writes thermodynamic properties to a file.
 *
 * Writes the current value and the average over a specified number of time
 * steps of values to a file with the file extension '.res'.
 * The following values are written to the file:
 * - Simulation time step
 * - time since the simulation started (dimensionless)
 * - Average potential Energy
 * - Pressure
 * - BetaTrans
 * - BetaRot
 */
class ResultWriter : public PluginBase {
public:
	ResultWriter() : _U_pot_acc(nullptr), _p_acc(nullptr) {}
	~ResultWriter() {
		delete _U_pot_acc;
		delete _p_acc;
	}

	/** @brief Read in XML configuration for ResultWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
		<outputplugin name="ResultWriter">
			<writefrequency>INTEGER</writefrequency>				<!-- Frequency in which the output is written; Default: 1 -->
			<outputprefix>STRING</outputprefix>						<!-- Prefix of the output file; Default: "results" -->
			<accumulation_steps>INTEGER</accumulation_steps>		<!-- Result is accumulated over the specified steps; Default: 1000 -->
			<writeprecision>UINTEGER</writeprecision>				<!-- Precision of output can be set here; Default: 5 -->
		</outputplugin>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);

	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    );

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);
	
	std::string getPluginName() {
		return std::string("ResultWriter");
	}
	static PluginBase* createInstance() { return new ResultWriter(); }

private:
	std::ofstream _resultStream;
	long _writeFrequency;
	int _writePrecision;
	std::string _outputPrefix;
	Accumulator<double> *_U_pot_acc;
	Accumulator<double> *_p_acc;
};

#endif  // SRC_IO_RESULTWRITER_H_
