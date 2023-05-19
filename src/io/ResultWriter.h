#ifndef SRC_IO_RESULTWRITER_H_
#define SRC_IO_RESULTWRITER_H_

#include "plugins/PluginBase.h"
#include "utils/Accumulator.h"

#include <string>


/** @brief Writes thermodynamic properties to a file.
 *
 * Writes the current value and the average over a specified number of time
 * steps of values to a file with the file extension '.res'.
 * The following values are written to the file:
 * - Simulation time step
 * - time since the simulation started (dimensionless)
 * - Average potential energy
 * - Average kinetic energy
 * - Average kinetic translational energy
 * - Average kinetic rotational energy
 * - Average Pressure
 * - Isochoric heat capacity
 * - Global number of particles
 */
class ResultWriter : public PluginBase {
public:
	ResultWriter() = default;
	~ResultWriter() override = default;


	/** @brief Read in XML configuration for ResultWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
		<outputplugin name="ResultWriter">
			<writefrequency>UINTEGER</writefrequency>				<!-- Frequency in which the output is written; Default: 1000 -->
			<outputprefix>STRING</outputprefix>						<!-- Prefix of the output file; Default: "mardyn" -->
			<writeprecision>UINTEGER</writeprecision>				<!-- Precision of output can be set here; Default: 5 -->
		</outputplugin>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
				 Domain *domain, unsigned long simstep) override;

	void finish(ParticleContainer* /* particleContainer */,
				DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}
	
	std::string getPluginName() override {return std::string("ResultWriter");}
	static PluginBase* createInstance() { return new ResultWriter(); }

private:
	uint64_t _writeFrequency{1000UL};
	unsigned int _writePrecision{5};
	unsigned int _writeWidth{20};
	std::string _outputPrefix{"mardyn"};
	uint64_t _numSamples{0UL};
	double _uPot_acc{0.0F};
	double _uKin_acc{0.0F};
	double _uKinTrans_acc{0.0F};
	double _uKinRot_acc{0.0F};
	double _p_acc{0.0F};
};

#endif  // SRC_IO_RESULTWRITER_H_
