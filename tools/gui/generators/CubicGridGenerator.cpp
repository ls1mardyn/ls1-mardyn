/*
 * CubicGridGenerator.cpp
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt, eckhardw
 */

#include "CubicGridGenerator.h"
#include "Parameters/ParameterWithIntValue.h"
#include "common/MardynConfigurationParameters.h"
#include "common/PrincipalAxisTransform.h"
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
	MDGenerator("CubicGridGenerator"), _numMolecules(4), _molarDensity(0.6),
	_temperature(1.5) {

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
			new ParameterWithDoubleValue("molarDensity", "Molar density [mol/l]",
					"molar density in mol/l", Parameter::LINE_EDIT,
					false, _molarDensity));
	tab->addParameter(
			new ParameterWithIntValue("numMolecules", "Number of Molecules",
					"Total number of Molecules", Parameter::LINE_EDIT,
					false, _numMolecules));
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
	if (id == "numMolecules") {
		_numMolecules = static_cast<ParameterWithIntValue*> (p)->getValue();
	} else if (id == "molarDensity") {
		_molarDensity = static_cast<ParameterWithDoubleValue*> (p)->getValue();
	} else if (id == "temperature") {
		_temperature = static_cast<ParameterWithDoubleValue*> (p)->getValue();
	} else if (id.find("component1") != std::string::npos) {
		std::string part = remainingSubString(".", id);
		ComponentParameters::setParameterValue(_components[0], p, part);
	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);
	}
}


void CubicGridGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
	// 1 mol/l = 0.6022 / nm^3 = 0.0006022 / Ang^3 = 0.089236726516 / a0^3
	double parts_per_a0 = _molarDensity * MDGenerator::molPerL_2_mardyn;
	double volume = _numMolecules / parts_per_a0;
	_simBoxLength[0] = pow(volume, 1./3.);
	_simBoxLength[1] = _simBoxLength[0];
	_simBoxLength[2] = _simBoxLength[0];

	domain->setCurrentTime(0);

	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_temperature);
	domain->setGlobalLength(0, _simBoxLength[0]);
	domain->setGlobalLength(1, _simBoxLength[1]);
	domain->setGlobalLength(2, _simBoxLength[2]);

	for (unsigned int i = 0; i < _components.size(); i++) {
		Component component = _components[i];
		if (_configuration.performPrincipalAxisTransformation()) {
			principalAxisTransform(component);
		}
		domain->addComponent(component);
	}
	domain->setepsilonRF(1e+10);
}


unsigned long CubicGridGenerator::readPhaseSpace(ParticleContainer* particleContainer,
		std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {

	int numMoleculesPerDimension = pow(_numMolecules / 2, 1./3.);

	int id = 1;
	double spacing = _simBoxLength[0] / numMoleculesPerDimension;
	double origin = spacing / 4.; // origin of the first DrawableMolecule

//	_logger->info() << "Spacing=" << spacing << " numMolPerDim=" << numMoleculesPerDimension <<
//			" origin=" << origin << endl;

	for (int i = 0; i < numMoleculesPerDimension; i++) {
		for (int j = 0; j < numMoleculesPerDimension; j++) {
			for (int k = 0; k < numMoleculesPerDimension; k++) {
				vector<double> velocity = getRandomVelocity(_temperature);
				Molecule m(id, 0, origin +
						i * spacing, origin + j * spacing, origin + k * spacing, // position
						velocity[0], -velocity[1], velocity[2], // velocity
						0, 0, 0, 0, 0, 0, 0, &_components);
				particleContainer->addParticle(m);
				id++;
			}
		}
	}

	origin = spacing / 4. * 3.; // origin of the first DrawableMolecule

	for (int i = 0; i < numMoleculesPerDimension; i++) {
		for (int j = 0; j < numMoleculesPerDimension; j++) {
			for (int k = 0; k < numMoleculesPerDimension; k++) {
				vector<double> velocity = getRandomVelocity(_temperature);
				Molecule m(id, 0, origin +
						i * spacing, origin + j * spacing, origin + k * spacing, // position
						velocity[0], -velocity[1], velocity[2], // velocity
						0, 0, 0, 0, 0, 0, 0, &_components);
				particleContainer->addParticle(m);
				id++;
			}
		}
	}
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

