/*
 * MardynConfiguration.h
 *
 * @Date: 11.09.2011
 * @Author: eckhardw
 */

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "OutputConfiguration.h"
#include <string>

/**
 * This class represents the configuration of a scenario for Mardyn.
 */
class MardynConfiguration {

public:

	static const std::string OutputFormat_XML;
	static const std::string OutputFormat_LEGACY;

	enum OutputFormat {
		XML, LEGACY
	};

	static const std::string ParallelisationType_NONE;
	static const std::string ParallelisationType_DOMAINDECOMPOSITION;
	static const std::string ParallelisationType_KDDECOMPOSITION;

	enum ParallelisationType {
		NONE, DOMAINDECOMPOSITION, KDDECOMPOSITION
	};

	static const std::string ContainerType_LINKEDCELLS;
	static const std::string ContainerType_ADAPTIVELINKEDCELLS;

	enum ContainerType {
		LINKEDCELLS, ADAPTIVELINKEDCELLS
	};

private:
	double _cutoffRadius;
	double _LJcutoffRadius;
	double _timestepLength;

	OutputFormat _outputFormat;
	std::string _scenarioName;

	bool _performPrincipalAxisTransformation;
	bool _NVE;

	ParallelisationType _parallelisationType;
	ContainerType _containerType;

	bool _hasResultWriter;
	OutputConfiguration _resultWriterConfig;
	bool _hasStatisticsWriter;
	OutputConfiguration _statisticsWriterConfig;
	bool _hasVTKMoleculeWriter;
	OutputConfiguration _vtkMoleculeWriterConfig;
	bool _hasVTKGridWriter;
	OutputConfiguration _vtkGridWriterConfig;

public:
	MardynConfiguration();

	virtual ~MardynConfiguration();

	double getCutoffRadius() const;

	void setCutoffRadius(double cutoffRadius);

	OutputFormat getOutputFormat() const;
	std::string getOutputFormatString() const;
	void setOutputFormat(std::string outputFormat);

	std::string getScenarioName() const;

	void setScenarioName(std::string scenarioName);

	double getTimestepLength() const;

	void setTimestepLength(double timestepLength);

	bool performPrincipalAxisTransformation() const;

	void setPerformPrincipalAxisTransformation(bool _performPrincipalAxisTransformation);

	double getLJCutoffRadius() const;

	void setLJCutoffRadius(double _LJcutoffRadius);

	bool isNVE() const;

    void setNVE(bool _NVE);

    void setContainerType(std::string _containerType);
    std::string getContainerTypeString() const;

    void setParallelisationType(std::string _parallelisationType);
    std::string getParallelisationTypeString() const;

    bool hasResultWriter() const;
    void setHasResultWriter(bool _hasResultWriter);
    OutputConfiguration& getResultWriterConfig();
    const OutputConfiguration& getResultWriterConfig() const;

    bool hasStatisticsWriter() const;
    void setHasStatisticsWriter(bool _hasStatisticsWriter);
    OutputConfiguration& getStatisticsWriterConfig();
    const OutputConfiguration& getStatisticsWriterConfig() const;

    bool hasVTKGridWriter() const;
    void setHasVTKGridWriter(bool _hasVTKGridWriter);
    OutputConfiguration& getVtkGridWriterConfig();
    const OutputConfiguration& getVtkGridWriterConfig() const;

    bool hasVTKMoleculeWriter() const;
    void setHasVTKMoleculeWriter(bool _hasVTKMoleculeWriter);
    OutputConfiguration& getVtkMoleculeWriterConfig();
    const OutputConfiguration& getVtkMoleculeWriterConfig() const;

};

#endif /* CONFIGURATION_H_ */
