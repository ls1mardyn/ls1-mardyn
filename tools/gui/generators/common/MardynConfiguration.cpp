/*
 * MardynConfiguration.cpp
 *
 * @Date: 11.09.2011
 * @Author: eckhardw
 */

#include "MardynConfiguration.h"
#include <iostream>
#include <cstdlib>


const std::string MardynConfiguration::OutputFormat_XML("XML");
const std::string MardynConfiguration::OutputFormat_LEGACY("ls1 LEGACY");

const std::string MardynConfiguration::ParallelisationType_NONE("None");
const std::string MardynConfiguration::ParallelisationType_DOMAINDECOMPOSITION("DomainDecomposition");
const std::string MardynConfiguration::ParallelisationType_KDDECOMPOSITION("KDDecomposition");

const std::string MardynConfiguration::ContainerType_LINKEDCELLS("LinkedCells");
const std::string MardynConfiguration::ContainerType_ADAPTIVELINKEDCELLS("AdaptiveLinkedCells");


MardynConfiguration::MardynConfiguration() :
	_cutoffRadius(1), _LJcutoffRadius(1), _timestepLength(1), _outputFormat(LEGACY),
	        _performPrincipalAxisTransformation(false),
	    	_NVE(false), _parallelisationType(KDDECOMPOSITION), _containerType(LINKEDCELLS),
	    	_hasResultWriter(false), _resultWriterConfig("ResultWriter"),
	    	_hasStatisticsWriter(false), _statisticsWriterConfig("StatisticsWriter"),
	    	_hasVTKMoleculeWriter(false), _vtkMoleculeWriterConfig("VTKMoleculeWriter"),
	    	_hasVTKGridWriter(false), _vtkGridWriterConfig("VTKGridWriter")
{
}

double MardynConfiguration::getCutoffRadius() const {
	return _cutoffRadius;
}

MardynConfiguration::OutputFormat MardynConfiguration::getOutputFormat() const {
	return _outputFormat;
}

std::string MardynConfiguration::getOutputFormatString() const {
	if (_outputFormat == XML) {
		return OutputFormat_XML;
	} else {
		return OutputFormat_LEGACY;
	}
}

std::string MardynConfiguration::getScenarioName() const {
	return _scenarioName;
}

double MardynConfiguration::getTimestepLength() const {
	return _timestepLength;
}

void MardynConfiguration::setCutoffRadius(double cutoffRadius) {
	this->_cutoffRadius = cutoffRadius;
}

void MardynConfiguration::setOutputFormat(std::string outputFormat) {
	if (outputFormat == OutputFormat_LEGACY) {
		_outputFormat = LEGACY;
	} else if (outputFormat == OutputFormat_XML) {
		_outputFormat = XML;
	} else {
		std::cout << "MardynConfiguration::setOutputFormat() Error: Invalid outputFormat: " << outputFormat << std::endl;
		exit(-1);
	}
}

void MardynConfiguration::setScenarioName(std::string scenarioName) {
	this->_scenarioName = scenarioName;
}

void MardynConfiguration::setTimestepLength(double timestepLength) {
	this->_timestepLength = timestepLength;
}

MardynConfiguration::~MardynConfiguration() {
}

bool MardynConfiguration::performPrincipalAxisTransformation() const {
	return _performPrincipalAxisTransformation;
}

double MardynConfiguration::getLJCutoffRadius() const {
	return _LJcutoffRadius;
}

bool MardynConfiguration::isNVE() const {
	return _NVE;
}

std::string MardynConfiguration::getContainerTypeString() const {
	if (_containerType == ADAPTIVELINKEDCELLS) {
		return ContainerType_ADAPTIVELINKEDCELLS;
	} else {
		return ContainerType_LINKEDCELLS;
	}
}

bool MardynConfiguration::hasResultWriter() const {
	return _hasResultWriter;
}

bool MardynConfiguration::hasStatisticsWriter() const {
	return _hasStatisticsWriter;
}

bool MardynConfiguration::hasVTKGridWriter() const {
	return _hasVTKGridWriter;
}

bool MardynConfiguration::hasVTKMoleculeWriter() const {
	return _hasVTKMoleculeWriter;
}

std::string MardynConfiguration::getParallelisationTypeString() const {
	if (_parallelisationType == NONE) {
		return ParallelisationType_NONE;
	} else if (_parallelisationType == KDDECOMPOSITION) {
		return ParallelisationType_KDDECOMPOSITION;
	} else {
		return ParallelisationType_DOMAINDECOMPOSITION;
	}
}

OutputConfiguration& MardynConfiguration::getResultWriterConfig() {
	return _resultWriterConfig;
}

OutputConfiguration& MardynConfiguration::getStatisticsWriterConfig() {
	return _statisticsWriterConfig;
}

OutputConfiguration& MardynConfiguration::getVtkGridWriterConfig() {
	return _vtkGridWriterConfig;
}

OutputConfiguration& MardynConfiguration::getVtkMoleculeWriterConfig() {
	return _vtkMoleculeWriterConfig;
}

const OutputConfiguration& MardynConfiguration::getResultWriterConfig() const {
	return _resultWriterConfig;
}

const OutputConfiguration& MardynConfiguration::getStatisticsWriterConfig() const {
	return _statisticsWriterConfig;
}

const OutputConfiguration& MardynConfiguration::getVtkGridWriterConfig() const {
	return _vtkGridWriterConfig;
}

const OutputConfiguration& MardynConfiguration::getVtkMoleculeWriterConfig() const {
	return _vtkMoleculeWriterConfig;
}

void MardynConfiguration::setContainerType(std::string containerType) {
	if (containerType == ContainerType_LINKEDCELLS) {
		_containerType = LINKEDCELLS;
	} else if (containerType == ContainerType_ADAPTIVELINKEDCELLS) {
		_containerType = ADAPTIVELINKEDCELLS;
	} else {
		std::cout << "MardynConfiguration::setContainerType() Error: Invalid containerType: " << containerType << std::endl;
		exit(-1);
	}
}

void MardynConfiguration::setHasResultWriter(bool _hasResultWriter) {
	this->_hasResultWriter = _hasResultWriter;
}

void MardynConfiguration::setHasStatisticsWriter(bool _hasStatisticsWriter) {
	this->_hasStatisticsWriter = _hasStatisticsWriter;
}

void MardynConfiguration::setHasVTKGridWriter(bool _hasVTKGridWriter) {
	this->_hasVTKGridWriter = _hasVTKGridWriter;
}

void MardynConfiguration::setHasVTKMoleculeWriter(bool _hasVTKMoleculeWriter) {
	this->_hasVTKMoleculeWriter = _hasVTKMoleculeWriter;
}

void MardynConfiguration::setParallelisationType(std::string parallelisationType) {
	if (parallelisationType == ParallelisationType_NONE) {
		_parallelisationType = NONE;
	} else if (parallelisationType == ParallelisationType_KDDECOMPOSITION) {
		_parallelisationType = KDDECOMPOSITION;
	} else if (parallelisationType == ParallelisationType_DOMAINDECOMPOSITION) {
		_parallelisationType = DOMAINDECOMPOSITION;
	}
	else {
		std::cout << "MardynConfiguration::setParallelisationType() Error: Invalid parallelisationType: " << parallelisationType << std::endl;
		exit(-1);
	}
}

void MardynConfiguration::setNVE(bool _NVE) {
	this->_NVE = _NVE;
}

void MardynConfiguration::setLJCutoffRadius(double _LJcutoffRadius) {
	this->_LJcutoffRadius = _LJcutoffRadius;
}

void MardynConfiguration::setPerformPrincipalAxisTransformation(bool performPrincipalAxisTransformation) {
	_performPrincipalAxisTransformation = performPrincipalAxisTransformation;
}


