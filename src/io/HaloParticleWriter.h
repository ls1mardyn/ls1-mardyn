#ifndef SRC_IO_HALOPARTICLEWRITER_H_
#define SRC_IO_HALOPARTICLEWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

/**
 * Prints the halo particles for each process in a separate files.
 * mainly useful for debugging purposes.
 * uses the afterForces step
 */
class HaloParticleWriter : public PluginBase {
public:

    HaloParticleWriter() = default;
	~HaloParticleWriter() override = default;


	/** @brief Read in XML configuration for HaloParticleWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="HaloParticleWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <incremental>INTEGER</incremental>
	     <appendTimestamp>BOOL</appendTimestamp>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain) override;
	void afterForces(ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, unsigned long simstep) override;

	void endStep(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep) override {
	}

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override;

	std::string getPluginName() override {
		return std::string("HaloParticleWriter");
	}
	static PluginBase* createInstance() { return new HaloParticleWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency{1000};
	bool _incremental{true};
	bool _appendTimestamp{false};
};

#endif  // SRC_IO_HALOPARTICLEWRITER_H_
