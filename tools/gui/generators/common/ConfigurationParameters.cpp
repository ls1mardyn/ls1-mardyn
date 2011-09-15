/*
 * ConfigurationParameters.cpp
 *
 * @Date: 03.08.2011
 * @Author: eckhardw
 */

#include "ConfigurationParameters.h"
#include "Parameters/ParameterWithDoubleValue.h"
#include "Parameters/ParameterWithChoice.h"

#include <vector>

ConfigurationParameters::ConfigurationParameters()
: ParameterCollection("ConfigurationParameters", "Mardyn/ls1 Configuration", "Mardyn/ls1 Configuration",
			Parameter::BUTTON, false, false) {

	std::vector<string> choices;
	choices.push_back("xml");
	choices.push_back("legacy ls1");
	addParameter(new ParameterWithChoice("ConfigurationParameters.fileformat", "Configuration File Format",
			"The format of the configuration file (xml not yet supported)", Parameter::LIST, false, choices, choices[1]));

	addParameter(new ParameterWithDoubleValue("ConfigurationParameters.timesteplength", "Length of timestep", "Length of timestep",
			Parameter::LINE_EDIT, false, 0.01));

	addParameter(new ParameterWithDoubleValue("ConfigurationParameters.cutoffradius", "Cutoff radius", "Cutoff Radius, i.e. cell length",
				Parameter::LINE_EDIT, false, 5.0));
}

ConfigurationParameters::~ConfigurationParameters() {
	// TODO Auto-generated destructor stub
}
