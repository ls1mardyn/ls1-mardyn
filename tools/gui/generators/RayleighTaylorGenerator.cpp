/*
 * RayleighTaylorGenerator.cpp
 *
 *  Created on: June, 2012
 *      Author: nagashim
 */

#include "RayleighTaylorGenerator.h"
#include "common/MardynConfigurationParameters.h"
#include "common/PrincipalAxisTransform.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Parameters/ParameterWithDoubleValue.h"
#include "Tokenize.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/Timer.h"

#include <cmath>
#include <climits>
#include <iostream>

#ifndef MARDYN
extern "C" {

Generator* create_generator() {
	return new RayleighTaylorGenerator();
}

void destruct_generator(Generator* generator) {
	delete generator;
}
}
#endif

RayleighTaylorGenerator::RayleighTaylorGenerator() :
					MDGenerator("RayleighTaylorGenerator") {
	numOfMolecules = 14000;
	L1 = 144;
	L2 = 60;
	L3 = 60;
	epsilon_A = 1;
	epsilon_B = 1;
	sigma_A = 1;
	sigma_B = 1;
	q_A = q_B = 0.5;
	m_A = m_B = 23;
	N = 14000;
	r_cut = 6;
	delta_t = 0.001;
	T = 0.1;
	G = 0.175;
	h = 1;
	p_max = 4;
	skal = 992.573;


	_temperature = 300. / 315774.5;
	_components.resize(1);
	_components[0].addLJcenter(0, 0, 0, 1.0, 1.0, 1.0, 0.0, false);
}

RayleighTaylorGenerator::~RayleighTaylorGenerator() {
}

void RayleighTaylorGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
	domain->setCurrentTime(0);

	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_temperature);

	domain->setGlobalLength(0, L1);
	domain->setGlobalLength(1, L2);
	domain->setGlobalLength(2, L3);

	/*
	domain->setGlobalLength(0, simBoxLength[0]);
	domain->setGlobalLength(1, simBoxLength[1]);
	domain->setGlobalLength(2, simBoxLength[2]);
	*/

	_logger->debug() << "RayleighTaylorGenerator: set global length=[" << L1
	                 << "," << L2 << "," <<  L3 << "]" << endl;

	/*
	_logger->debug() << "RayleighTaylorGenerator: set global length=[" << simBoxLength[0]
	                 << "," << simBoxLength[1] << "," <<  simBoxLength[2] << "]" << endl;
	                 */

	for (unsigned int i = 0; i < _components.size(); i++) {
		Component component = _components[i];
		if (_configuration.performPrincipalAxisTransformation()) {
			principalAxisTransform(component);
		}
		domain->addComponent(component);
	}
	domain->setepsilonRF(1e+10);
}

unsigned long RayleighTaylorGenerator::readPhaseSpace(
		ParticleContainer* particleContainer,
		std::list<ChemicalPotential>* lmu, Domain* domain,
		DomainDecompBase* domainDecomp) {

	Timer inputTimer;
	inputTimer.start();
	_logger->info() << "Reading phase space file (RayleighTaylorGenerator)." << endl;

	srand(1);
	vector<double> bBoxMin;
	vector<double> bBoxMax;

	bBoxMin.resize(3);
	bBoxMax.resize(3);
	for (int i = 0; i < 3; i++) {
		bBoxMin[i] = domainDecomp->getBoundingBoxMin(i, domain);
		bBoxMax[i] = domainDecomp->getBoundingBoxMax(i, domain);
	}

	_logger->info()
							<< "OneCLJGenerator  generating cluster distribution. " << " T "
							<< _temperature << " #molecules " << numOfMolecules << " rho_gas "
							<< " rho_fluid " << endl;

	vector<Component>& dcomponents = domain->getComponents();
	vector<unsigned long> partsPerComp;
	partsPerComp.resize(1);
	particleContainer->update();
	particleContainer->deleteOuterParticles();
	domain->setglobalNumMolecules(
			domainDecomp->countMolecules(particleContainer, partsPerComp));

	for (unsigned int i = 0; i < partsPerComp.size(); i++) {
		dcomponents[i].setNumMolecules(partsPerComp[i]);
		domain->setglobalRotDOF(
				partsPerComp[i]
				             * dcomponents[i].getRotationalDegreesOfFreedom());
	}

	domain->setglobalRho(
			domain->getglobalNumMolecules() / (L1 * L2 * L3));
/*
	domain->setglobalRho(
			domain->getglobalNumMolecules() / (simBoxLength[0]
			          * simBoxLength[1] * simBoxLength[2]));
*/
	removeMomentum(particleContainer, _components);

	inputTimer.stop();
	_logger->info() << "Initial IO took:                 " << inputTimer.get_etime() << " sec" << endl;
	return 0;
}

vector<ParameterCollection*> RayleighTaylorGenerator::getParameters() {
	vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("RayleighTaylorGenerator", "RayleighTaylorGenerator",
			"Parameters of RayleighTaylorGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithDoubleValue("L1", "L1",
					"L1", Parameter::LINE_EDIT, false, L1));

	tab->addParameter(
			new ParameterWithDoubleValue("L2", "L2",
					"L2", Parameter::LINE_EDIT, false, L2));

	tab->addParameter(
			new ParameterWithDoubleValue("L3", "L3",
					"L3", Parameter::LINE_EDIT, false, L3));

	tab->addParameter(
			new ParameterWithDoubleValue("epsilon_A", "epsilon_A",
					"epsilon_A", Parameter::LINE_EDIT, false, epsilon_A));

	tab->addParameter(
			new ParameterWithDoubleValue("epsilon_B", "epsilon_B",
					"epsilon_B", Parameter::LINE_EDIT, false, epsilon_B));

	tab->addParameter(
			new ParameterWithDoubleValue("sigma_A", "sigma_A",
					"sigma_A", Parameter::LINE_EDIT, false, sigma_A));

	tab->addParameter(
			new ParameterWithDoubleValue("sigma_B", "sigma_B",
					"sigma_B", Parameter::LINE_EDIT, false, sigma_B));

	tab->addParameter(
			new ParameterWithDoubleValue("q_A", "q_A",
					"q_A", Parameter::LINE_EDIT, false, q_A));
	tab->addParameter(
			new ParameterWithDoubleValue("q_B", "q_B",
					"q_B", Parameter::LINE_EDIT, false, q_B));

	tab->addParameter(
			new ParameterWithDoubleValue("m_A", "m_A",
					"m_A", Parameter::LINE_EDIT, false, m_A));

	tab->addParameter(
			new ParameterWithDoubleValue("m_B", "m_B",
					"m_B", Parameter::LINE_EDIT, false, m_B));

	tab->addParameter(
			new ParameterWithIntValue("N", "N",
					"N", Parameter::LINE_EDIT, false, N));

	/*
	tab->addParameter(
			new ParameterWithIntValue("numOfMolecules", "numOfMolecules",
					"Number of Molecules", Parameter::LINE_EDIT, false, numOfMolecules));*/

	tab->addParameter(
			new ParameterWithDoubleValue("r_cut", "r_cut",
					"r_cut", Parameter::LINE_EDIT, false, r_cut));

	tab->addParameter(
			new ParameterWithDoubleValue("delta_t", "delta_t",
					"delta_t", Parameter::LINE_EDIT, false, delta_t));

	tab->addParameter(
			new ParameterWithDoubleValue("T", "T",
					"T", Parameter::LINE_EDIT, false, T));

	/*
	tab->addParameter(
			new ParameterWithDoubleValue("temperature", "Temperature [K]",
					"Temperature in the domain in Kelvin", Parameter::LINE_EDIT,
					false, _temperature / MDGenerator::kelvin_2_mardyn ));
	 */

	tab->addParameter(
			new ParameterWithDoubleValue("G", "G",
					"G", Parameter::LINE_EDIT, false, G));

	tab->addParameter(
			new ParameterWithDoubleValue("h", "h",
					"h", Parameter::LINE_EDIT, false, h));

	tab->addParameter(
			new ParameterWithDoubleValue("p_max", "p_max",
					"p_max", Parameter::LINE_EDIT, false, p_max));

	tab->addParameter(
			new ParameterWithDoubleValue("skal", "skal",
					"skal", Parameter::LINE_EDIT, false, skal));

	tab->addParameter(
			new ComponentParameters("component1", "component1", "Set up the parameters of component 1",
					_components[0]));

	return parameters;
}

void RayleighTaylorGenerator::setParameter(Parameter* p) {

	string id = p->getNameId();

	if (id == "L1") {
		L1 = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: L1: " << L1
				<< endl;

	} else if (id == "L2") {
		L2 = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: L2: " << L2
				<< endl;

	} else if (id == "L3") {
		L2 = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: L3: " << L3
				<< endl;

	} else if (id == "epsilon_A") {
		epsilon_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: epsilon_A: " << epsilon_A
				<< endl;

	} else if (id == "epsilon_B") {
		epsilon_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: epsilon_B: " << epsilon_B
				<< endl;

	} else if (id == "sigma_A") {
		sigma_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: sigma_A: " << sigma_A
				<< endl;

	} else if (id == "sigma_B") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: sigma_B: " << sigma_B
				<< endl;

	} else if (id == "m_A") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: m_A: " << m_A
				<< endl;

	} else if (id == "m_B") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: m_B: " << m_B
				<< endl;

	} else if (id == "N") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: N: " << N
				<< endl;
/*
	} else if (id == "numOfMolecules") {
		numOfMolecules = static_cast<ParameterWithIntValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: numOfMolecules: " << numOfMolecules
				<< endl;
*/

	} else if (id == "r_cut") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: r_cut: " << r_cut
				<< endl;

	} else if (id == "delta_t") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: delta_t: " << delta_t
				<< endl;

	} else if (id == "numSphereSizes") {
		numSphereSizes = static_cast<ParameterWithIntValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: numSphereSizes: " << numSphereSizes
				<< endl;

	} else if (id == "T") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: T: " << T
				<< endl;
/*
	} else if (id == "temperature") {
		_temperature = static_cast<ParameterWithDoubleValue*> (p)->getValue() * MDGenerator::kelvin_2_mardyn;
*/

	} else if (id == "G") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: G: " << G
				<< endl;

	} else if (id == "h") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: h: " << h
				<< endl;

	} else if (id == "p_max") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet:p_max: " << p_max
				<< endl;

	} else if (id == "skal") {
		sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet:skal: " << skal
				<< endl;

	}  else if (id.find("component1") != std::string::npos) {
		std::string part = id.substr(11);
		ComponentParameters::setParameterValue(_components[0], p, part);

	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);

	} else {
		std::cout << "UNKOWN Parameter: id = " << p->getNameId() << " value= " << p->getStringValue() << endl;
		exit(-1);
	}
}

bool RayleighTaylorGenerator::validateParameters() {
	bool valid = true;

	if (_configuration.getScenarioName() == "") {
		valid = false;
		_logger->error() << "ScenarioName not set!" << endl;
	}

	if (_configuration.getOutputFormat() == MardynConfiguration::XML) {
		valid = false;
		_logger->error() << "OutputFormat XML not yet supported!" << endl;
	}

	double L[3];
	L[0] = L1;
	L[1] = L2;
	L[2] = L3;

	for (int i = 0; i < 3; i++) {
		//cout << "\nsimBox is " << simBoxLength[i]<< flush;
		if (L[i] < 2.0 * _configuration.getCutoffRadius()) {
			valid = false;
			_logger->error() << "Cutoff radius is too big (there would be only 1 cell in the domain!)" << endl;
			_logger->error() << "Cutoff radius=" << _configuration.getCutoffRadius()
									<< " domain size=" << L[i]/*simBoxLength[i]*/ << endl;
		}
	}

	return valid;
}
