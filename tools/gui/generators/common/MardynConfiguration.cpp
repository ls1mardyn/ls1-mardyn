/*
 * MardynConfiguration.cpp
 *
 * @Date: 11.09.2011
 * @Author: eckhardw
 */

#include "MardynConfiguration.h"

MardynConfiguration::MardynConfiguration()
: _cutoffRadius(1), _LJcutoffRadius(1), _timestepLength(1), _outputFormat(LEGACY),
  _performPrincipalAxisTransformation(false)
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

bool MardynConfiguration::performPrincipalAxisTransformation() const
{
    return _performPrincipalAxisTransformation;
}

double MardynConfiguration::getLJCutoffRadius() const
{
    return _LJcutoffRadius;
}

void MardynConfiguration::setLJCutoffRadius(double _LJcutoffRadius)
{
    this->_LJcutoffRadius = _LJcutoffRadius;
}

void MardynConfiguration::setPerformPrincipalAxisTransformation(bool performPrincipalAxisTransformation)
{
    _performPrincipalAxisTransformation = performPrincipalAxisTransformation;
}


