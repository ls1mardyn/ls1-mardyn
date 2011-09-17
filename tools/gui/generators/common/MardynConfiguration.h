/*
 * MardynConfiguration.h
 *
 * @Date: 11.09.2011
 * @Author: eckhardw
 */

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include <string>

/**
 * This class represents the configuration of a scenario for Mardyn.
 */
class MardynConfiguration {

public:

	enum OutputFormat {
		XML, LEGACY
	};

private:
	double _cutoffRadius;
	double _LJcutoffRadius;
	double _timestepLength;

	OutputFormat _outputFormat;
	std::string _scenarioName;

	bool _performPrincipalAxisTransformation;

public:
	MardynConfiguration();

	virtual ~MardynConfiguration();

	double getCutoffRadius() const;

	OutputFormat getOutputFormat() const;

	std::string getScenarioName() const;

	double getTimestepLength() const;

	void setCutoffRadius(double cutoffRadius);

	void setOutputFormat(OutputFormat outputFormat);

	void setScenarioName(std::string scenarioName);

	void setTimestepLength(double timestepLength);

	bool performPrincipalAxisTransformation() const;

	void setPerformPrincipalAxisTransformation(bool _performPrincipalAxisTransformation);

	double getLJCutoffRadius() const;

	void setLJCutoffRadius(double _LJcutoffRadius);
};

#endif /* CONFIGURATION_H_ */
