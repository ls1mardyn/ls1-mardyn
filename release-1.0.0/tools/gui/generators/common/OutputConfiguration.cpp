/*
 * OutputConfiguration.cpp
 *
 * @Date: 19.09.2011
 * @Author: eckhardw
 */

#include "OutputConfiguration.h"

OutputConfiguration::OutputConfiguration(std::string name)
: _name(name), _outputPrefix("default"), _outputFrequency(1) {
}

std::string OutputConfiguration::getName() const
{
    return _name;
}

unsigned int OutputConfiguration::getOutputFrequency() const
{
    return _outputFrequency;
}

std::string OutputConfiguration::getOutputPrefix() const
{
    return _outputPrefix;
}

void OutputConfiguration::setName(std::string _name)
{
    this->_name = _name;
}

void OutputConfiguration::setOutputFrequency(unsigned int _outputFrequency)
{
    this->_outputFrequency = _outputFrequency;
}

void OutputConfiguration::setOutputPrefix(std::string _outputPrefix)
{
    this->_outputPrefix = _outputPrefix;
}

OutputConfiguration::~OutputConfiguration() {
}
