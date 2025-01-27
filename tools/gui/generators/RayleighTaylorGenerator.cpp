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
		MDGenerator("RayleighTaylorGenerator"), _components(*(global_simulation->getEnsemble()->getComponents())) {

	_L1 		= 144.;// * sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength;
	_L2 		= 60.;//  * sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength;
	_L3 		= 60.;//  * sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength;
	_n_1 		= 14;
	_n_2 		= 5;
	_n_3 		= 5;
	_epsilon_A 	= 1.;// * epsilon_tilde * MDGenerator::kelvin_2_mardyn;
	_epsilon_B 	= 1.;// * epsilon_tilde * MDGenerator::kelvin_2_mardyn;
	_sigma_A 	= 1.;// * sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength;
	_sigma_B 	= 1.;// * sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength;
	_q_A 		=  0.5;// * MDGenerator::unitCharge_2_mardyn;
	_q_B 		= -0.5;// * MDGenerator::unitCharge_2_mardyn;
	_m_A = _m_B = 23.;// * m_tilde * MDGenerator::unitMass_2_mardyn;
	_T 			= 0.1;// * epsilon_tilde * MDGenerator::kelvin_2_mardyn / MDGenerator::boltzmann_constant_kB;
	numSphereSizes = 1.;
	_components.resize(2);
	_components[0].addCharge(0.,0.,0.,0.,_q_A);
	_components[1].addCharge(0.,0.,0.,0.,_q_B);
	_components[0].addLJcenter(0, 0, 0,_m_A, _epsilon_A, _sigma_A, 0., false);
	_components[1].addLJcenter(0, 0, 0,_m_A, _epsilon_B, _sigma_B, 0., false);
}

RayleighTaylorGenerator::~RayleighTaylorGenerator() {
}

void RayleighTaylorGenerator::readPhaseSpaceHeader(Domain* domain, double /*timestep*/) {
	//domain->setCurrentTime(0);

	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_T);
	domain->setGlobalLength(0, _L1);
	domain->setGlobalLength(1, _L2);
	domain->setGlobalLength(2, _L3);
	if (_configuration.performPrincipalAxisTransformation()) {
		for (unsigned int i = 0; i < _components.size(); i++) {
			principalAxisTransform(_components[i]);
		}
	}
	domain->setepsilonRF(1e+10);

	global_simulation->getEnsemble()->setComponentLookUpIDs();
}

unsigned long RayleighTaylorGenerator::readPhaseSpace(ParticleContainer* particleContainer,
		Domain* domain, DomainDecompBase* domainDecomp) {

	global_simulation->timers()->start("REYLEIGH_TAYLOR_GENERATOR_INPUT");
	_logger->info() << "Reading phase space file (RayleighTaylorGenerator)." << std::endl;

	_components[0].updateMassInertia();
	_components[1].updateMassInertia();
	// RayleighTaylor scenario is always a binary mixture of positive and negative particles.
	//	if (_binaryMixture) {
	//		_components[1].updateMassInertia();
	//	}
	//	std::cout << "from RayleighTaylor::readPhaseSpace " << std::flush;
	int id = 1;
	int numOfAddedMolecule = 0;

	double lowerBound_z = 1./4.* _L3;
	double upperBound_z = 3./4.* _L3;
	int componentType;
	double space_x = _L1 / (_n_1 + 1);
	double space_y = _L2 / (_n_2 + 1);
	double space_z = _L3 / (_n_3 + 1);

	for(int idx = 1; idx <= _n_1; idx++){
		for(int idy = 1; idy <= _n_2; idy++){
			for(int idz = 1; idz <=_n_3; idz++){

				double x = space_x * idx;
				double y = space_y * idy;
				double z = space_z * idz;

				if((lowerBound_z <= z)&&(z < upperBound_z)){
					componentType = 1;
				}else{
					componentType = 0;
				}

				if (domainDecomp->procOwnsPos(x,y,z, domain)) {
					addMolecule(x, y, z, id, componentType,particleContainer);
					numOfAddedMolecule ++;
				}

				id++;
			}
		}
	}

	removeMomentum(particleContainer, _components);

	unsigned long int globalNumMolecules = particleContainer->getNumberOfParticles();
	std::cout << "numberOfParticles is " << particleContainer->getNumberOfParticles()<< "\n"<<std::flush;
#ifdef ENABLE_PERSISTENT
	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), globalNumMolecules);
	collComm.persistent();
	collComm.get(globalNumMolecules);
#else
	domainDecomp->collCommInit(1);

	domainDecomp->collCommAppendUnsLong(globalNumMolecules);
	domainDecomp->collCommAllreduceSum();
	globalNumMolecules = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
#endif

	domain->setglobalNumMolecules(globalNumMolecules);
	global_simulation->timers()->stop("REYLEIGH_TAYLOR_GENERATOR_INPUT");
	global_simulation->timers()->setOutputString("REYLEIGH_TAYLOR_GENERATOR_INPUT", "Initial IO took:                 ");
	_logger->info() << "Initial IO took:                 " << global_simulation->timers()->getTime("REYLEIGH_TAYLOR_GENERATOR_INPUT") << " sec" << std::endl;
	return id;
}

std::vector<ParameterCollection*> RayleighTaylorGenerator::getParameters() {
	std::vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("RayleighTaylorGenerator", "RayleighTaylorGenerator",
			"Parameters of RayleighTaylorGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithDoubleValue("SimulationBoxSize(X)", "SimulationBoxSize(X)",
					"SimulationBoxSize(X)", Parameter::LINE_EDIT,
					false, _L1));

	tab->addParameter(
			new ParameterWithDoubleValue("SimulationBoxSize(Y)", "SimulationBoxSize(Y)",
					"SimulationBoxSize(Y)", Parameter::LINE_EDIT,
					false, _L2));

	tab->addParameter(
			new ParameterWithDoubleValue("SimulationBoxSize(Z)", "SimulationBoxSize(Z)",
					"SimulationBoxSize(Z)", Parameter::LINE_EDIT,
					false, _L3));

	tab->addParameter(
			new ParameterWithIntValue("NumOfParticlesAlongX", "NumOfParticlesAlongX",
					"NumOfParticlesAlongX", Parameter::LINE_EDIT, false, _n_1));

	tab->addParameter(
			new ParameterWithIntValue("NumOfParticlesAlongY", "NumOfParticlesAlongY",
					"NumOfParticlesAlongY", Parameter::LINE_EDIT, false, _n_2));

	tab->addParameter(
			new ParameterWithIntValue("NumOfParticlesAlongZ", "NumOfParticlesAlongZ",
					"NumOfParticlesAlongZ", Parameter::LINE_EDIT, false, _n_3));

	tab->addParameter(
			new ParameterWithDoubleValue("Temperature", "Temperature",
					"Temperature", Parameter::LINE_EDIT, false,
					_T));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentA.epsilon", "ComponentA.epsilon",
					"ComponentA.epsilon", Parameter::LINE_EDIT, false,
					_epsilon_A));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentB.epsilon", "ComponentB.epsilon",
					"ComponentB.epsilon", Parameter::LINE_EDIT, false,
					_epsilon_B));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentA.sigma", "ComponentA.sigma",
					"ComponentA.sigma", Parameter::LINE_EDIT, false,
					_sigma_A));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentB.sigma", "ComponentB.sigma",
					"ComponentB.sigma", Parameter::LINE_EDIT, false,
					_sigma_B));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentA.charge", "ComponentA.charge",
					"ComponentA.charge", Parameter::LINE_EDIT, false,
					_q_A));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentB.charge", "ComponentB.charge",
					"ComponentB.charge", Parameter::LINE_EDIT, false,
					_q_B));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentA.mass", "ComponentA.mass",
					"ComponentA.mass", Parameter::LINE_EDIT, false,
					_m_A));

	tab->addParameter(
			new ParameterWithDoubleValue("ComponentB.mass", "ComponentB.mass",
					"ComponentB.mass", Parameter::LINE_EDIT, false,
					_m_B));

	return parameters;
}

//TODO:reference values!
void RayleighTaylorGenerator::setParameter(Parameter* p) {

	std::string id = p->getNameId();

	if (id == "SimulationBoxSize(X)") {
		_L1 = static_cast<ParameterWithDoubleValue*> (p)->
				getValue() * (sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength);
		std::cout << "OneCenterLJRayleighTaylor: SimulationBoxSize(X): "
				<< _L1
				<< std::endl;

	} else if (id == "SimulationBoxSize(Y)") {
		_L2 = static_cast<ParameterWithDoubleValue*> (p)->
				getValue() * (sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength);
		std::cout << "OneCenterLJRayleighTaylor: SimulationBoxSize(Y): "
				<< _L2
				<< std::endl;

	} else if (id == "SimulationBoxSize(Z)") {
		_L3 = static_cast<ParameterWithDoubleValue*> (p)->
				getValue() * (sigma_tilde * MDGenerator::angstroem_2_atomicUnitLength);
		std::cout << "OneCenterLJRayleighTaylor: SimulationBoxSize(Z): "
				<< _L3
				<< std::endl;

	} else if (id == "NumOfParticlesAlongX") {
		_n_1 = static_cast<ParameterWithIntValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: NumOfParticlesAlongX: " << _n_1
				<< std::endl;

	} else if (id == "NumOfParticlesAlongY") {
		_n_2 = static_cast<ParameterWithIntValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: NumOfParticlesAlongY: " << _n_2
				<< std::endl;

	} else if (id == "NumOfParticlesAlongZ") {
		_n_3 = static_cast<ParameterWithIntValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: NumOfParticlesAlonZ: " << _n_3
				<< std::endl;

	} else if (id == "numSphereSizes") {
		numSphereSizes = static_cast<ParameterWithIntValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: numSphereSizes: " << numSphereSizes
				<< std::endl;

	} else if (id == "Temperature") {
		_T = static_cast<ParameterWithDoubleValue*> (p)->getValue()* (epsilon_tilde * MDGenerator::kelvin_2_mardyn / MDGenerator::boltzmann_constant_kB);
		std::cout << "OneCenterLJRayleighTaylor: Temperature: " << _T
				<< std::endl;

	} else if (id == "ComponentA.charge") {
		_q_A = static_cast<ParameterWithDoubleValue*> (p)->getValue() * MDGenerator::unitCharge_2_mardyn;
		std::cout << "OneCenterLJRayleighTaylor: ComponentB.charge: " << _q_A
				<< std::endl;
		_components[0].deleteCharge();
		_components[0].addCharge(0.,0.,0.,0.,_q_A);


	} else if (id == "ComponentB.charge") {
		_q_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: ComponentB.charge: " << _q_B
				<< std::endl;
		_components[1].deleteCharge();
		_components[1].addCharge(0.,0.,0.,0.,_q_B);

	} else if (id == "ComponentA.epsilon") {
		_epsilon_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: ComponentA.epsilon: " << _epsilon_A
				<< std::endl;
		_components[0].deleteLJCenter();
		_components[0].addLJcenter(0, 0, 0,_m_A, _epsilon_A, _sigma_A, 0., false);

	} else if (id == "ComponentB.epsilon") {
		_epsilon_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: ComponentB.epsilon: " << _epsilon_B
				<< std::endl;
		_components[1].deleteLJCenter();
		_components[1].addLJcenter(0, 0, 0,_m_B, _epsilon_B, _sigma_B, 0., false);

	} else if (id == "ComponentA.sigma") {
		_sigma_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: ComponentA.sigma: " << _sigma_A
				<< std::endl;
		_components[0].deleteLJCenter();
		_components[0].addLJcenter(0, 0, 0,_m_A, _epsilon_A, _sigma_A, 0., false);

	} else if (id == "ComponentB.sigma") {
		_sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: ComponentB.sigma: " << _sigma_B
				<< std::endl;
		_components[1].deleteLJCenter();
		_components[1].addLJcenter(0, 0, 0,_m_B, _epsilon_B, _sigma_B, 0., false);

	} else if (id == "ComponentA.mass") {
		_m_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		std::cout << "OneCenterLJRayleighTaylor: ComponentA.mass: " << _m_A
				<< std::endl;
		_components[0].deleteLJCenter();
		_components[0].addLJcenter(0, 0, 0,_m_A, _epsilon_A, _sigma_A, 0., false);

	} else if (id == "ComponentB.mass") {
		_m_B = static_cast<ParameterWithDoubleValue*> (p)->getValue() * (m_tilde * MDGenerator::unitMass_2_mardyn);
		std::cout << "OneCenterLJRayleighTaylor: ComponentB.mass: " << _m_B
				<< std::endl;
		_components[1].deleteLJCenter();
		_components[1].addLJcenter(0, 0, 0,_m_B, _epsilon_B, _sigma_B, 0., false);

	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);

	} else {
		std::cout << "UNKOWN Parameter: id = " << p->getNameId() << " value= " << p->getStringValue() << std::endl;
		exit(-1);
	}
}

void RayleighTaylorGenerator::addMolecule(
		double x, double y, double z, unsigned long id, int componentType, ParticleContainer* particleContainer) {
	std::vector<double> velocity = getRandomVelocity(_T);

	//double orientation[4] = {1, 0, 0, 0}; // default: in the xy plane
	// rotate by 30° along the vector (1/1/0), i.e. the angle bisector of x and y axis
	// o = cos 30° + (1 1 0) * sin 15°
	double orientation[4];
	getOrientation(15, 10, orientation);

	double I[3] = {0.,0.,0.};
	I[0] = _components[componentType].I11();
	I[1] = _components[componentType].I22();
	I[2] = _components[componentType].I33();
	// I don't understand this one, I simply hope it initialized angular velocity correctly.
	//  Copied from animake - initialize anular velocity
	double w[3];
	for(int d=0; d < 3; d++) {
		w[d] = (I[d] == 0)? 0.0: ((randdouble(0,1) > 0.5)? 1: -1) *
				sqrt(2.0* randdouble(0,1)* _T / I[d]);
		//w[d] = w[d] * MDGenerator::fs_2_mardyn;
	}
	//End Copy

	Molecule m(id, &(_components[componentType]), x, y, z, // position
			velocity[0], -velocity[1], velocity[2], // velocity
			orientation[0], orientation[1], orientation[2], orientation[3],
			w[0], w[1], w[2]);
	if (particleContainer->isInBoundingBox(m.r_arr().data())) {
		bool inChecked = true;
		particleContainer->addParticle(m, inChecked);
	}
}

bool RayleighTaylorGenerator::validateParameters() {
	bool valid = true;

	if (_configuration.getScenarioName() == "") {
		valid = false;
		_logger->error() << "ScenarioName not set!" << std::endl;
	}

	if (_configuration.getOutputFormat() == MardynConfiguration::XML) {
		valid = false;
		_logger->error() << "OutputFormat XML not yet supported!" << std::endl;
	}

	double L[3];
	L[0] = _L1;
	L[1] = _L2;
	L[2] = _L3;

	for (int i = 0; i < 3; i++) {
		if (L[i] < 2.0 * _configuration.getCutoffRadius()) {
			valid = false;
			_logger->error()
					<< "Cutoff radius is too big (there would be only 1 cell in the domain!)" << std::endl;
			_logger->error() << "Cutoff radius="
					<< _configuration.getCutoffRadius() << " domain size=" << L[i] << std::endl;
		}
	}

	return valid;
}

