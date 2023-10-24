#pragma once

#include <string>

#include "plugins/PluginBase.h"

/**
 * Prints the CommunicationPartners for each rank in a separate file.
 * Useful for debugging.
 */
class CommunicationPartnerWriter : public PluginBase {
public:

	CommunicationPartnerWriter() {}
	~CommunicationPartnerWriter() {}


	/** @brief Read in XML configuration for CommunicationPartnerWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="CommunicationPartnerWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <incremental>INTEGER</incremental>
	     <appendTimestamp>INTEGER</appendTimestamp>
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
		return std::string("CommunicationPartnerWriter");
	}
	static PluginBase* createInstance() { return new CommunicationPartnerWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool	_incremental;
	bool	_appendTimestamp;
};

