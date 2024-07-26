/*
 * ExamplePlugin.cpp
 *
 *  Created on: 4 Jun 2018
 *      Author: tchipevn
 */

#include "ExamplePlugin.h"

#include "utils/mardyn_assert.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"

#include <string>


void ExamplePlugin::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	_message = "Your code would be called here.";
	xmlconfig.getNodeValue("message", _message);
	Log::global_log->info() << "Output prefix: " << _message << std::endl;

	_displaySelector = WhereToDisplay::ALL;
	std::string displaySelectorString;
	xmlconfig.getNodeValue("where_to_display", displaySelectorString);
	Log::global_log->info() << "where to display: " << displaySelectorString
			<< std::endl;

	const char* str = displaySelectorString.c_str();
	if (std::strcmp(str, "all") == 0) {
		_displaySelector = WhereToDisplay::ALL;
		Log::global_log->info() << "Displaying at all plugin positions."
				<< std::endl;
	} else if (std::strcmp(str, "beforeEventNewTimestep") == 0) {
		_displaySelector = WhereToDisplay::BEFORE_EVENT_NEW_TIMESTEP;
		Log::global_log->info() << "Displaying at beforeEventNewTimestep."
				<< std::endl;
	} else if (std::strcmp(str, "beforeForces") == 0) {
		_displaySelector = WhereToDisplay::BEFORE_FORCES;
		Log::global_log->info() << "Displaying at beforeForces." << std::endl;
	} else if (std::strcmp(str, "afterForces") == 0) {
		_displaySelector = WhereToDisplay::AFTER_FORCES;
		Log::global_log->info() << "Displaying at afterForces." << std::endl;
	} else if (std::strcmp(str, "endStep") == 0) {
		_displaySelector = WhereToDisplay::END_STEP;
		Log::global_log->info() << "Displaying at endStep." << std::endl;
	} else if (std::strcmp(str, "init") == 0) {
		_displaySelector = WhereToDisplay::AT_INIT;
		Log::global_log->info() << "Displaying at init." << std::endl;
	} else if (std::strcmp(str, "finish") == 0) {
		_displaySelector = WhereToDisplay::AT_FINISH;
		Log::global_log->info() << "Displaying at finish." << std::endl;
	} else {
		Log::global_log->error()
				<< "Unknown option specified to ExamplePlugin::where_to_display."
				<< std::endl;
		Log::global_log->error()
				<< "Valid options are: all, beforeEventNewTimestep, beforeForces, afterForces, endStep, init, finish."
				<< std::endl;
		mardyn_exit(11);
	}
}

void ExamplePlugin::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	if (_displaySelector == WhereToDisplay::ALL
			or _displaySelector == WhereToDisplay::AT_INIT) {
		Log::global_log->info() << "ExamplePlugin::init                           " << _message << std::endl;
	}
}

void ExamplePlugin::beforeEventNewTimestep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector
							== WhereToDisplay::BEFORE_EVENT_NEW_TIMESTEP)) {
		Log::global_log->info() << "ExamplePlugin::beforeEventNewTimestep         " << _message << std::endl;
	}
}

void ExamplePlugin::beforeForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector == WhereToDisplay::BEFORE_FORCES)) {
		Log::global_log->info() << "ExamplePlugin::beforeForces                   " << _message << std::endl;
	}
}

void ExamplePlugin::afterForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector == WhereToDisplay::AFTER_FORCES)) {
		Log::global_log->info() << "ExamplePlugin::afterForces                    " << _message << std::endl;
	}
}

void ExamplePlugin::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector == WhereToDisplay::END_STEP)) {
		Log::global_log->info() << "ExamplePlugin::endStep                        " << _message << std::endl;
	}
}

void ExamplePlugin::finish(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	if (_displaySelector == WhereToDisplay::ALL
			or _displaySelector == WhereToDisplay::AT_FINISH) {
		Log::global_log->info() << "ExamplePlugin::finish                         " << _message << std::endl;
	}
}
