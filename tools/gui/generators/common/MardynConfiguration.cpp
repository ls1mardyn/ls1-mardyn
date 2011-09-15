/*
 * MardynConfiguration.cpp
 *
 * @Date: 11.09.2011
 * @Author: eckhardw
 */

#include "MardynConfiguration.h"

MardynConfiguration::MardynConfiguration()
: _cutoffRadius(3.0), _timestepLength(0.1), _outputFormat(LEGACY)
{


}

double MardynConfiguration::getCutoffRadius() const
{
    return _cutoffRadius;
}

MardynConfiguration::OutputFormat MardynConfiguration::getOutputFormat() const
{
    return _outputFormat;
}

std::string MardynConfiguration::getScenarioName() const
{
    return _scenarioName;
}

double MardynConfiguration::getTimestepLength() const
{
    return _timestepLength;
}

void MardynConfiguration::setCutoffRadius(double cutoffRadius)
{
    this->_cutoffRadius = cutoffRadius;
}

void MardynConfiguration::setOutputFormat(OutputFormat outputFormat)
{
    this->_outputFormat = outputFormat;
}

void MardynConfiguration::setScenarioName(std::string scenarioName)
{
    this->_scenarioName = scenarioName;
}

void MardynConfiguration::setTimestepLength(double timestepLength)
{
    this->_timestepLength = timestepLength;
}

MardynConfiguration::~MardynConfiguration() {
}
