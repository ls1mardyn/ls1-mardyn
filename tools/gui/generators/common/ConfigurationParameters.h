/*
 * ConfigurationParameters.h
 *
 * @Date: 03.08.2011
 * @Author: eckhardw
 */

#ifndef CONFIGURATIONPARAMETERS_H_
#define CONFIGURATIONPARAMETERS_H_

#include "Parameters/ParameterCollection.h"

/**
 * This class represents a collection of parameters corresponding to the
 * parameters in the config file in Mardyn/ls1.
 */
class ConfigurationParameters: public ParameterCollection {

public:

	ConfigurationParameters();

	virtual ~ConfigurationParameters();
};

#endif /* CONFIGURATIONPARAMETERS_H_ */
