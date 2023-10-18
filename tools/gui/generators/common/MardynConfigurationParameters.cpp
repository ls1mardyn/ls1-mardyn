/*
 * MardynConfigurationParameters.cpp
 *
 * @Date: 03.08.2011
 * @Author: eckhardw
 */

#include "MardynConfigurationParameters.h"
#include "MardynConfiguration.h"
#include "Parameters/ParameterWithDoubleValue.h"
#include "Parameters/ParameterWithChoice.h"
#include "Parameters/ParameterWithStringValue.h"
#include "Parameters/ParameterWithBool.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Tokenize.h"
#include "../MDGenerator.h"

#include <vector>
#include <iostream>


MardynConfigurationParameters::MardynConfigurationParameters(const MardynConfiguration& other)
: ParameterCollection("ConfigurationParameters", "Mardyn/ls1 Configuration", "Mardyn/ls1 Configuration",
		Parameter::BUTTON, false, false) {

	addParameter(new ParameterWithStringValue("ConfigurationParameters.ScenarioName", "Scenario name", "Name of the scenario (determines the base name of the output files)",
			Parameter::LINE_EDIT, false, other.getScenarioName()));

	std::vector<std::string> choices;
	choices.push_back(MardynConfiguration::OutputFormat_XML);
	choices.push_back(MardynConfiguration::OutputFormat_LEGACY);
	addParameter(new ParameterWithChoice("ConfigurationParameters.fileformat", "Configuration File Format",
			"The format of the configuration file (xml not yet supported)", Parameter::COMBOBOX, false, choices,
			other.getOutputFormatString()));

	addParameter(new ParameterWithBool("ConfigurationParameters.NVE", "NVE Ensemble", "Perform simulation in the NVE Ensemble.",
				Parameter::CHECKBOX, false, other.isNVE()));

	addParameter(new ParameterWithDoubleValue("ConfigurationParameters.timesteplength", "Length of timestep [fs]", "Length of timestep",
			Parameter::LINE_EDIT, false, other.getTimestepLength() / MDGenerator::fs_2_mardyn));

	addParameter(new ParameterWithDoubleValue("ConfigurationParameters.cutoffradius", "Cutoff radius [Angstrom]", "Cutoff Radius, i.e. cell length",
				Parameter::LINE_EDIT, false, other.getCutoffRadius() / MDGenerator::angstroem_2_atomicUnitLength));

	addParameter(new ParameterWithDoubleValue("ConfigurationParameters.LJcutoffradius", "LJ Cutoff radius [Angstrom]", "Lennard-Jones Cutoff Radius, i.e. cell length",
				Parameter::LINE_EDIT, false, other.getLJCutoffRadius() / MDGenerator::angstroem_2_atomicUnitLength));

	addParameter(new ParameterWithBool("ConfigurationParameters.principalAxisTrafo", "Perform principal axis transformation", "Perform a principal axis transformation for each component\nwhen writing the Mardyn phasespace file.",
				Parameter::CHECKBOX, false, other.performPrincipalAxisTransformation()));

	std::vector<std::string> parChoices;
	parChoices.push_back(MardynConfiguration::ParallelisationType_DOMAINDECOMPOSITION);
	parChoices.push_back(MardynConfiguration::ParallelisationType_KDDECOMPOSITION);
	parChoices.push_back(MardynConfiguration::ParallelisationType_NONE);
	addParameter(new ParameterWithChoice("ConfigurationParameters.ParallelisationType", "Parallelisation Type",
				"- DomainDecomposition: standard parallelisation method\n- KDDecomposition: parallelisation with dynamic load balancing (recommended)\n- NONE: no parallelisation",
				Parameter::COMBOBOX, false, parChoices, other.getParallelisationTypeString()));

	std::vector<std::string> containerChoices;
	containerChoices.push_back(MardynConfiguration::ContainerType_LINKEDCELLS);
	containerChoices.push_back(MardynConfiguration::ContainerType_ADAPTIVELINKEDCELLS);
	addParameter(new ParameterWithChoice("ConfigurationParameters.ContainerType", "Container Type",
					"molecule container type",
					Parameter::COMBOBOX, false, containerChoices, other.getContainerTypeString()));

	addParameter(new ParameterWithBool("ConfigurationParameters.hasResultWriter", "Create ResultWriter output", "Write the basic results (temperature, pressure) to a file.",
				Parameter::CHECKBOX, true, other.hasResultWriter()));
	if (other.hasResultWriter()) {
		addOutputConfigurationParameters(other.getResultWriterConfig(), "ResultWriter");
	}
	addParameter(new ParameterWithBool("ConfigurationParameters.hasStatisticsWriter", "Create StatisticsWriter output", "Write statistics (number of molecules per cell) to a file.",
			Parameter::CHECKBOX, true, other.hasStatisticsWriter()));
	if (other.hasStatisticsWriter()) {
		addOutputConfigurationParameters(other.getStatisticsWriterConfig(), "StatisticsWriter");
	}
	addParameter(new ParameterWithBool("ConfigurationParameters.hasVTKMoleculeWriter", "Create VTK output", "Write molecule data to a file for visualization VTK/Paraview.",
			Parameter::CHECKBOX, true, other.hasVTKMoleculeWriter()));
	if (other.hasVTKMoleculeWriter()) {
		addOutputConfigurationParameters(other.getVtkMoleculeWriterConfig(), "VTKWriter");
	}
	addParameter(new ParameterWithBool("ConfigurationParameters.hasVTKGridWriter", "Create VTK output (Celldata)", "Write cell data to a file for visualization VTK/Paraview.",
			Parameter::CHECKBOX, true, other.hasVTKGridWriter()));
	if (other.hasVTKGridWriter()) {
		addOutputConfigurationParameters(other.getVtkGridWriterConfig(), "VTKGridWriter");
	}
}


MardynConfigurationParameters::~MardynConfigurationParameters() {
}


void MardynConfigurationParameters::setParameterValue(MardynConfiguration& config, const Parameter* parameter, const std::string valueName) {
	if (valueName == "ScenarioName") {
		std::string scenarioName = dynamic_cast<const ParameterWithStringValue*>(parameter)->getStringValue();
		config.setScenarioName(scenarioName);
	} else if (valueName == "fileformat") {
		std::string choice = dynamic_cast<const ParameterWithChoice*>(parameter)->getStringValue();
		config.setOutputFormat(choice);
	} else if (valueName == "timesteplength") {
		double timestepLength = dynamic_cast<const ParameterWithDoubleValue*>(parameter)->getValue();
		config.setTimestepLength(timestepLength * MDGenerator::fs_2_mardyn);
	} else if (valueName == "cutoffradius") {
		double cutoffRadius = dynamic_cast<const ParameterWithDoubleValue*>(parameter)->getValue();
		config.setCutoffRadius(cutoffRadius * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "LJcutoffradius") {
		double cutoffRadius = dynamic_cast<const ParameterWithDoubleValue*>(parameter)->getValue();
		config.setLJCutoffRadius(cutoffRadius * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "principalAxisTrafo") {
		bool performTrafo = dynamic_cast<const ParameterWithBool*>(parameter)->getValue();
		config.setPerformPrincipalAxisTransformation(performTrafo);
	} else if (valueName == "NVE") {
		bool NVE = dynamic_cast<const ParameterWithBool*>(parameter)->getValue();
		config.setNVE(NVE);
	} else if (valueName == "ParallelisationType") {
		std::string choice = dynamic_cast<const ParameterWithChoice*>(parameter)->getStringValue();
		config.setParallelisationType(choice);
	} else if (valueName == "ContainerType") {
		std::string choice = dynamic_cast<const ParameterWithChoice*>(parameter)->getStringValue();
		config.setContainerType(choice);
	} else if (valueName == "hasStatisticsWriter") {
		bool hasStatisticsWriter = dynamic_cast<const ParameterWithBool*>(parameter)->getValue();
		config.setHasStatisticsWriter(hasStatisticsWriter);
	} else if (valueName == "hasVTKMoleculeWriter") {
		bool hasVTKMoleculeWriter = dynamic_cast<const ParameterWithBool*>(parameter)->getValue();
		config.setHasVTKMoleculeWriter(hasVTKMoleculeWriter);
	} else if (valueName == "hasVTKGridWriter") {
		bool hasVTKGridWriter = dynamic_cast<const ParameterWithBool*>(parameter)->getValue();
		config.setHasVTKGridWriter(hasVTKGridWriter);
	} else if (valueName == "hasResultWriter") {
		bool hasResultWriter = dynamic_cast<const ParameterWithBool*>(parameter)->getValue();
		config.setHasResultWriter(hasResultWriter);
	} else {
		const std::string firstPart = firstSubString(".", valueName);
		const std::string secondPart = remainingSubString(".", valueName);

		if (firstPart == "ResultWriter") {
			setOutputConfigurationParameter(config.getResultWriterConfig(), parameter, secondPart);
		} else if (firstPart == "StatisticsWriter") {
			setOutputConfigurationParameter(config.getStatisticsWriterConfig(), parameter, secondPart);
		} else if (firstPart == "VTKWriter") {
			setOutputConfigurationParameter(config.getVtkMoleculeWriterConfig(), parameter, secondPart);
		} else if (firstPart == "VTKGridWriter") {
			setOutputConfigurationParameter(config.getVtkGridWriterConfig(), parameter, secondPart);
		} else {
			std::cout << "Invalid Parameter in ConfigurationParameters::setParameterValue(): " << std::endl;
			parameter->print();
			exit(-1);
		}
	}
}


void MardynConfigurationParameters::addOutputConfigurationParameters(const OutputConfiguration& config, const std::string& basename) {
	addParameter(new ParameterWithStringValue("ConfigurationParameters." + basename + ".Prefix", config.getName() + " Output Prefix", "The file name will be prefixed by this value.",
			Parameter::LINE_EDIT, false, config.getOutputPrefix()));
	addParameter(new ParameterWithIntValue("ConfigurationParameters." + basename + ".Frequency", config.getName() + " Output Frequency", "The output file will be written every n-th timestep.",
			Parameter::LINE_EDIT, false, config.getOutputFrequency()));
}

void MardynConfigurationParameters::setOutputConfigurationParameter(OutputConfiguration& config, const Parameter* parameter, const std::string& valueName) {
	if (valueName == "Prefix") {
		std::string prefix = dynamic_cast<const ParameterWithStringValue*>(parameter)->getStringValue();
		config.setOutputPrefix(prefix);
	} else if (valueName == "Frequency") {
		int frequency = dynamic_cast<const ParameterWithIntValue*>(parameter)->getValue();
		config.setOutputFrequency(frequency);
	} else {
		std::cout << "Invalid Parameter in ConfigurationParameters::setParameterValue(): " << std::endl;
		parameter->print();
		exit(-1);
	}
}
