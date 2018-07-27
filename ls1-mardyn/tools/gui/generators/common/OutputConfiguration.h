/*
 * OutputConfiguration.h
 *
 * @Date: 19.09.2011
 * @Author: eckhardw
 */

#ifndef OUTPUTCONFIGURATION_H_
#define OUTPUTCONFIGURATION_H_

#include <string>

/**
 * contains the configuration for one Output Plugin of Mardyn
 */
class OutputConfiguration {

	std::string _name;
	std::string _outputPrefix;
	unsigned _outputFrequency;

public:

	OutputConfiguration(std::string name);

	virtual ~OutputConfiguration();

    std::string getName() const;

    unsigned int getOutputFrequency() const;

    std::string getOutputPrefix() const;

    void setName(std::string _name);

    void setOutputFrequency(unsigned int _outputFrequency);

    void setOutputPrefix(std::string _outputPrefix);

};

#endif /* OUTPUTCONFIGURATION_H_ */
