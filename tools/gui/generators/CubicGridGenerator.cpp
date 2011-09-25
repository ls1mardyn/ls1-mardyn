/*
 * CubicGridGenerator.cpp
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt, eckhardw
 */

#include "CubicGridGenerator.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Parameters/ParameterWithBool.h"
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
	_temperature(300. / 315774.5), _binaryMixture(false) {

	_components.resize(1);
	_components[0].addLJcenter(0, 0, 0, 1.0, 1.0, 1.0, 0.0, false);
	calculateSimulationBoxLength();
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
			new ParameterWithDoubleValue("temperature", "Temperature [K]",
					"Temperature in the domain in Kelvin", Parameter::LINE_EDIT,
					false, _temperature / MDGenerator::kelvin_2_mardyn ));
	tab->addParameter(
			new ComponentParameters("component1", "component1",
					"Set up the parameters of component 1", _components[0]));
	tab->addParameter(
			new ParameterWithBool("binaryMixture", "Binary Mixture",
					"Check this option to simulate a binary mixture.\n(A second component will be added.)",
					Parameter::CHECKBOX, true, _binaryMixture));
	if (_binaryMixture) {
		tab->addParameter(
				new ComponentParameters("component2", "component2",
						"Set up the parameters of component 2", _components[1]));
	}
	return parameters;
}


void CubicGridGenerator::setParameter(Parameter* p) {
	string id = p->getNameId();
	if (id == "numMolecules") {
		_numMolecules = static_cast<ParameterWithIntValue*> (p)->getValue();
		calculateSimulationBoxLength();
	} else if (id == "molarDensity") {
		_molarDensity = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		calculateSimulationBoxLength();
	} else if (id == "temperature") {
		_temperature = static_cast<ParameterWithDoubleValue*> (p)->getValue() * MDGenerator::kelvin_2_mardyn;
	} else if (id.find("component1") != std::string::npos) {
		std::string part = remainingSubString(".", id);
		ComponentParameters::setParameterValue(_components[0], p, part);
	} else if (id == "binaryMixture") {
		_binaryMixture = static_cast<ParameterWithBool*>(p)->getValue();
		if (_binaryMixture && _components.size() == 1) {
			_components.resize(2);
			_components[1].addLJcenter(0, 0, 0, 1.0, 1.0, 1.0, 5.0, false);
		} else if (!_binaryMixture && _components.size() == 2) {
			_components.resize(1);
		}
	} else if (id.find("component2") != std::string::npos) {
		std::string part = remainingSubString(".", id);
		ComponentParameters::setParameterValue(_components[1], p, part);
	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);
	}
}


void CubicGridGenerator::calculateSimulationBoxLength() {
	// 1 mol/l = 0.6022 / nm^3 = 0.0006022 / Ang^3 = 0.089236726516 / a0^3
	double parts_per_a0 = _molarDensity * MDGenerator::molPerL_2_mardyn;
	double volume = _numMolecules / parts_per_a0;
	_simBoxLength[0] = pow(volume, 1./3.);
	_simBoxLength[1] = _simBoxLength[0];
	_simBoxLength[2] = _simBoxLength[0];
}


void CubicGridGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
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
//
// create a body centered cubic layout, by creating by placing the molecules on the
// vertices of a regular grid, then shifting that grid by spacing/2 in all dimensions.

	int numMoleculesPerDimension = pow(_numMolecules / 2, 1./3.);
	_components[0].updateMassInertia();
	if (_binaryMixture) {
		_components[1].updateMassInertia();
	}

	int id = 1;
	double spacing = _simBoxLength[0] / numMoleculesPerDimension;
	double origin = spacing / 4.; // origin of the first DrawableMolecule

	for (int i = 0; i < numMoleculesPerDimension; i++) {
		for (int j = 0; j < numMoleculesPerDimension; j++) {
			for (int k = 0; k < numMoleculesPerDimension; k++) {
				vector<double> velocity = getRandomVelocity(_temperature);

				int componentType = 0;
				if (_binaryMixture) {
					componentType = randdouble(0, 1.999999);
				}

				double I[3] = {0.,0.,0.};
				I[0] = _components[componentType].I11();
				I[1] = _components[componentType].I22();
				I[2] = _components[componentType].I33();
				/*****  Copied from animake - initialize anular velocity *****/
				double w[3];
				for(int d=0; d < 3; d++) {
					w[d] = (I[d] == 0)? 0.0: ((randdouble(0,1) > 0.5)? 1: -1) *
							sqrt(2.0* randdouble(0,1)* _temperature / I[d]);
					w[d] = w[d] * MDGenerator::fs_2_mardyn;
				}
				/************************** End Copy **************************/

				Molecule m(id, componentType, origin +
						i * spacing, origin + j * spacing, origin + k * spacing, // position
						velocity[0], -velocity[1], velocity[2], // velocity
						1, 0, 0, 0, w[0], w[1], w[2], &_components);
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
				int componentType = 0;
				if (_binaryMixture) {
					componentType = randdouble(0, 1.999999);
				}

				double I[3] = {0.,0.,0.};
				I[0] = _components[componentType].I11();
				I[1] = _components[componentType].I22();
				I[2] = _components[componentType].I33();
				/*****  Copied from animake - initialize anular velocity *****/
				double w[3];
				for(int d=0; d < 3; d++) {
					w[d] = (I[d] == 0)? 0.0: ((randdouble(0,1) > 0.5)? 1: -1) *
							sqrt(2.0* randdouble(0,1)* _temperature / I[d]);
					w[d] = w[d] * MDGenerator::fs_2_mardyn;
				}
				/************************** End Copy **************************/

				Molecule m(id, componentType, origin +
						i * spacing, origin + j * spacing, origin + k * spacing, // position
						velocity[0], -velocity[1], velocity[2], // velocity
						1, 0, 0, 0, w[0], w[1], w[2], &_components);
				particleContainer->addParticle(m);
				id++;
			}
		}
	}
	return id;
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

	if (_simBoxLength[0] < 2. * _configuration.getCutoffRadius()) {
		valid = false;
		_logger->error() << "Cutoff radius is too big (there would be only 1 cell in the domain!)" << endl;
		_logger->error() << "Cutoff radius=" << _configuration.getCutoffRadius()
							<< " domain size=" << _simBoxLength[0] << endl;
	}
	return valid;
}

