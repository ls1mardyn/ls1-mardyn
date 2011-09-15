/*
 * MardynConfigurationParameters.h
 *
 * @Date: 03.08.2011
 * @Author: eckhardw
 */

#ifndef CONFIGURATIONPARAMETERS_H_
#define CONFIGURATIONPARAMETERS_H_

#include "Parameters/ParameterCollection.h"

class MardynConfiguration;

/**
 * This class represents a collection of parameters corresponding to the
 * parameters in the config file in Mardyn/ls1.
 */
class MardynConfigurationParameters: public ParameterCollection {

public:

	MardynConfigurationParameters(const MardynConfiguration& config);

	virtual ~MardynConfigurationParameters();

	/**
	 * set the value of parameter (which must belong logically to this collection)
	 */
	static void setParameterValue(MardynConfiguration& config, const Parameter* parameter, const std::string valueName);

};

#endif /* CONFIGURATIONPARAMETERS_H_ */
