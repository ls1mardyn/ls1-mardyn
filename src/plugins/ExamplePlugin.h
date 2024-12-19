/*
 * ExamplePlugin.h
 *
 *  Created on: 4 Jun 2018
 *      Author: tchipevn
 */

#ifndef SRC_PLUGINS_EXAMPLEPLUGIN_H_
#define SRC_PLUGINS_EXAMPLEPLUGIN_H_

#include "PluginBase.h"

#include <string>

/**
 * The purpose of this class is to show the basic usage of Plugins.
 * It just outputs a string that the user specified in the XML at the
 * plugin position that the user specified.
 */
class ExamplePlugin: public PluginBase {
private:
	enum class WhereToDisplay {
		ALL=0,
		BEFORE_EVENT_NEW_TIMESTEP = 1,
		BEFORE_FORCES = 2,
		AFTER_FORCES = 3,
		END_STEP = 4,
		AT_INIT = 5,
		AT_FINISH = 6,
	};

public:
	ExamplePlugin() {}
	virtual ~ExamplePlugin() {}

	void init(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

    void readXML(XMLfileUnits& xmlconfig);

	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	);

    void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep);

    void finish(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp, Domain* domain);

    std::string getPluginName() {
    	return std::string("ExamplePlugin");
    }

	static PluginBase* createInstance() { return new ExamplePlugin(); }

private:
    std::string _message;
	unsigned long _writeFrequency;
	WhereToDisplay _displaySelector;
};

#endif /* SRC_PLUGINS_EXAMPLEPLUGIN_H_ */
