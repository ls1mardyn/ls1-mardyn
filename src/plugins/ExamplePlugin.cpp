/*
 * ExamplePlugin.cpp
 *
 *  Created on: 4 Jun 2018
 *      Author: tchipevn
 */

#include "ExamplePlugin.h"
#include "Simulation.h"

#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"


void ExamplePlugin::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	_message = "Your code would be called here.";
	xmlconfig.getNodeValue("message", _message);
	global_log->info() << "Output prefix: " << _message << std::endl;

	_displaySelector = WhereToDisplay::ALL;
	std::string displaySelectorString;
	xmlconfig.getNodeValue("where_to_display", displaySelectorString);
	global_log->info() << "where to display: " << displaySelectorString
			<< std::endl;

	const char* str = displaySelectorString.c_str();
	if (strcmp(str, "all") == 0) {
		_displaySelector = WhereToDisplay::ALL;
		global_log->info() << "Displaying at all plugin positions."
				<< std::endl;
	} else if (strcmp(str, "beforeEventNewTimestep") == 0) {
		_displaySelector = WhereToDisplay::BEFORE_EVENT_NEW_TIMESTEP;
		global_log->info() << "Displaying at beforeEventNewTimestep."
				<< std::endl;
	} else if (strcmp(str, "beforeForces") == 0) {
		_displaySelector = WhereToDisplay::BEFORE_FORCES;
		global_log->info() << "Displaying at beforeForces." << std::endl;
	} else if (strcmp(str, "afterForces") == 0) {
		_displaySelector = WhereToDisplay::AFTER_FORCES;
		global_log->info() << "Displaying at afterForces." << std::endl;
	} else if (strcmp(str, "endStep") == 0) {
		_displaySelector = WhereToDisplay::END_STEP;
		global_log->info() << "Displaying at endStep." << std::endl;
	} else if (strcmp(str, "init") == 0) {
		_displaySelector = WhereToDisplay::AT_INIT;
		global_log->info() << "Displaying at init." << std::endl;
	} else if (strcmp(str, "finish") == 0) {
		_displaySelector = WhereToDisplay::AT_FINISH;
		global_log->info() << "Displaying at finish." << std::endl;
	} else {
		global_log->error()
				<< "Unknown option specified to ExamplePlugin::where_to_display."
				<< std::endl;
		global_log->error()
				<< "Valid options are: all, beforeEventNewTimestep, beforeForces, afterForces, endStep, init, finish."
				<< std::endl;
		Simulation::exit(11);
	}
}

void ExamplePlugin::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	if (_displaySelector == WhereToDisplay::ALL
			or _displaySelector == WhereToDisplay::AT_INIT) {
		global_log->info() << "ExamplePlugin::init                           " << _message << std::endl;
	}
}

void ExamplePlugin::beforeEventNewTimestep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector
							== WhereToDisplay::BEFORE_EVENT_NEW_TIMESTEP)) {
		global_log->info() << "ExamplePlugin::beforeEventNewTimestep         " << _message << std::endl;
	}
}

void ExamplePlugin::beforeForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector == WhereToDisplay::BEFORE_FORCES)) {
		global_log->info() << "ExamplePlugin::beforeForces                   " << _message << std::endl;
	}
}

void ExamplePlugin::afterForces(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector == WhereToDisplay::AFTER_FORCES)) {
		global_log->info() << "ExamplePlugin::afterForces                    " << _message << std::endl;
	}
}

void ExamplePlugin::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	if (simstep % _writeFrequency == 0
			and (_displaySelector == WhereToDisplay::ALL
					or _displaySelector == WhereToDisplay::END_STEP)) {
		global_log->info() << "ExamplePlugin::endStep                        " << _message << std::endl;
	}
}

void ExamplePlugin::finish(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	if (_displaySelector == WhereToDisplay::ALL
			or _displaySelector == WhereToDisplay::AT_FINISH) {
		global_log->info() << "ExamplePlugin::finish                         " << _message << std::endl;
	}
}
