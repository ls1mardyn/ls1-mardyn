/*
 * CubicGridGenerator.cpp
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt
 */

#include "CubicGridGenerator.h"
#include "Parameters/ParameterWithIntValue.h"
#include "common/MardynConfigurationParameters.h"
#include "Tokenize.h"
#include <cstring>


extern "C" {

	Generator* create_generator() {
		return new CubicGridGenerator();
	}

	void destruct_generator(Generator* generator) {
		delete generator;
	}
}


CubicGridGenerator::CubicGridGenerator() :
	MDGenerator("CubicGridGenerator"), numMoleculesX(4), numMoleculesY(4),
			numMoleculesZ(3), _temperature(1.5), origin(0.5, 0.5, 0.5) {

	dx = 1.0;
	dy = 1.0;
	dz = 1.0;

	_components.resize(1);
	_components[0].addLJcenter(0, 0, 0, 1.0, 1.0, 1.0, 5.0, false);
}

vector<ParameterCollection*> CubicGridGenerator::getParameters() {
	vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("EqvGridParameters", "Parameters of EqvGridGenerator",
			"Parameters of EqvGridGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithIntValue("numMoleculesX", "Number of Molecules X",
					"Number of Molecules in X direction", Parameter::LINE_EDIT,
					false, numMoleculesX));
	tab->addParameter(
			new ParameterWithIntValue("numMoleculesY", "Number of Molecules Y",
					"Number of Molecules in Y direction", Parameter::LINE_EDIT,
					false, numMoleculesY));
	tab->addParameter(
			new ParameterWithIntValue("numMoleculesZ", "Number of Molecules Z",
					"Number of Molecules in Z direction", Parameter::LINE_EDIT,
					false, numMoleculesZ));
	tab->addParameter(
			new ParameterWithDoubleValue("temperature", "Temperature",
					"Temperature in the domain", Parameter::LINE_EDIT,
					false, _temperature ));
	tab->addParameter(
			new ComponentParameters("component1", "component1",
					"Set up the parameters of component 1", _components[0]));
	return parameters;
}


void CubicGridGenerator::setParameter(Parameter* p) {
	string id = p->getNameId();
	if (id == "numMoleculesX") {
		numMoleculesX = static_cast<ParameterWithIntValue*> (p)->getValue();
	} else if (id == "numMoleculesY") {
		numMoleculesY = static_cast<ParameterWithIntValue*> (p)->getValue();
	} else if (id == "numMoleculesZ") {
		numMoleculesZ = static_cast<ParameterWithIntValue*> (p)->getValue();
	} else if (id.find("component1") != std::string::npos) {
		std::string part = remainingSubString(".", id);
		ComponentParameters::setParameterValue(_components[0], p, part);
	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);
	}
}


void CubicGridGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
	simBoxLength[0] = dx * numMoleculesX;
	simBoxLength[1] = dy * numMoleculesY;
	simBoxLength[2] = dz * numMoleculesZ;

	domain->setCurrentTime(0);

	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_temperature);
	domain->setGlobalLength(0, simBoxLength[0]);
	domain->setGlobalLength(1, simBoxLength[1]);
	domain->setGlobalLength(2, simBoxLength[2]);

	for (unsigned int i = 0; i < _components.size(); i++) {
		domain->addComponent(_components[i]);
	}
	domain->setepsilonRF(1e+10);
}


unsigned long CubicGridGenerator::readPhaseSpace(ParticleContainer* particleContainer,
		std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {

	int id = 1;
	for (int i = 0; i < numMoleculesX; i++) {
		for (int j = 0; j < numMoleculesY; j++) {
			for (int k = 0; k < numMoleculesZ; k++) {
				Molecule m(id, 0, origin.x + i * dx,
						origin.y + j * dy, origin.z + k * dz, 0.1 * id,
						-0.2 * id, 0.3 * id, 0, 0, 0, 0, 0, 0, 0, &_components);
				particleContainer->addParticle(m);
				id++;
			}
		}
	}
	_logger->debug() << "EqvGridGenerator: numParticles = "<< particleContainer->getNumberOfParticles() << endl;
}


bool CubicGridGenerator::validateParameters() {
	bool valid = true;

	if (_configuration.getScenarioName() == "") {
		valid = false;
		_logger->error() << "ScenarioName not set!" << endl;
	}

	if (_configuration.getOutputFormat() == MardynConfiguration::XML) {
		valid = false;
		_logger->error() << "OutputFormat XML not yet supported!" << endl;
	}
	return valid;
}

