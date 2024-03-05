#ifndef SRC_IO_ENERGYLOGWRITER_H_
#define SRC_IO_ENERGYLOGWRITER_H_

#include <string>

#include "plugins/PluginBase.h"


/** @brief This writer creates a global energy log file.
 *
 * The global energy log file will contain
 * - the global number of molecules \f$N\f$
 * - the global potential energy \f$U_\mathrm{pot}\f$
 * - the global kinetic energy \f$U_\mathrm{kin}\f$
 * - the global kinetic motional energy \f$U_\mathrm{kinTrans}\f$
 * - the global kinetic rotational energy \f$U_\mathrm{kinRot}\f$
 * - the global temperature \f$T\f$
 * - the global pressure \f$p\f$
 */
class EnergyLogWriter : public PluginBase {
public:
	EnergyLogWriter() : _outputFilename("global_energy.log"), _writeFrequency(1) {}
	~EnergyLogWriter() {}

	/** @brief Read in XML configuration for EnergyLogWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="EnergyLogWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputfilename>STRING</outputfilename>
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
	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() {
		return std::string("EnergyLogWriter");
	}
	static PluginBase* createInstance() { return new EnergyLogWriter(); }
private:
	std::string _outputFilename;
	unsigned long _writeFrequency;
};

#endif  // SRC_IO_ENERGYLOGWRITER_H_
