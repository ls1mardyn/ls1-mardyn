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
	_L1 = 144*MDGenerator::angstroem_2_atomicUnitLength;
	_L2 = 60 * MDGenerator::angstroem_2_atomicUnitLength;
	_L3 = 60 * MDGenerator::angstroem_2_atomicUnitLength;
	_epsilon_A = 1;
	_epsilon_B = 1;
	_sigma_A = 1;
	_sigma_B = 1;
	_q_A = _q_B = 0.5;
	_m_A = _m_B = 23 * MDGenerator::unitMass_2_mardyn;
	_N = 140;
	_r_cut = 6;
	_delta_t = 0.001;
	_T = 300 * MDGenerator::kelvin_2_mardyn;
	_G = 0.175;//TODO
	_h = 1;//TODO
	_p_max = 4;//TODO
	_skal = 992.573;//TODO

	//	_temperature = 300. / 315774.5;
	_components.resize(2);
	_components[0].addCharge(0.,0.,0.,_m_A,_q_A);
	_components[1].addCharge(0.,0.,0.,_m_B,_q_B);
	_components[0].addLJcenter(0, 0, 0,_m_A, _epsilon_A, _sigma_A, _r_cut, false);
	_components[1].addLJcenter(0, 0, 0,_m_A, _epsilon_B, _sigma_B, _r_cut, false);
}


RayleighTaylorGenerator::~RayleighTaylorGenerator() {
}

bool RayleighTaylorGenerator::getRandomPosition(
		double boxSize_x,double boxSize_y,double boxSize_z,
		double &x,double &y,double &z,int componentType){

	x = boxSize_x * ((double) rand() / (double) RAND_MAX);
	y = boxSize_y * ((double) rand() / (double) RAND_MAX);
	z = boxSize_z * ((double) rand() / (double) RAND_MAX);

	double lowerBound_z = (2-sqrt(2))/4.*boxSize_z;
	double upperBound_z = (2+sqrt(2))/4.*boxSize_z;
	double lowerBound_y = (2-sqrt(2))/4.*boxSize_y;
	double upperBound_y = (2+sqrt(2))/4.*boxSize_y;

	if(componentType==0){
		if((lowerBound_z <= z)&&(lowerBound_y<= y)&&
				(z < upperBound_z)&&(y < upperBound_y)){
			return true;
		}else if(getRandomPosition(boxSize_x,boxSize_y,boxSize_z,x,y,z,componentType)){
			return true;
		}
	}else if(componentType==1){
		if(!((lowerBound_z <= z)&&(lowerBound_y<= y)&&
				(z < upperBound_z)&&(y < upperBound_y))){
			return true;
		}else if(getRandomPosition(boxSize_x,boxSize_y,boxSize_z,x,y,z,componentType)){
			return true;
		}
	}
	else{// If componentType is not either 0 or 1
		std::cout<<"Invalid component type is detected in getRandomPosition!" << std::flush;
		return false;
	}

	return true;
}

void RayleighTaylorGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
	domain->setCurrentTime(0);

	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_T);
	domain->setGlobalLength(0, _L1);
	domain->setGlobalLength(1, _L2);
	domain->setGlobalLength(2, _L3);

	for (unsigned int i = 0; i < _components.size(); i++) {
		Component component = _components[i];
		if (_configuration.performPrincipalAxisTransformation()) {
			principalAxisTransform(component);
		}
		domain->addComponent(component);
	}
	domain->setepsilonRF(1e+10);
}


unsigned long RayleighTaylorGenerator::readPhaseSpace(ParticleContainer* particleContainer,
		std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {

	Timer inputTimer;
	inputTimer.start();
	_logger->info() << "Reading phase space file (RayleighTaylorGenerator)." << endl;

	_components[0].updateMassInertia();
	_components[1].updateMassInertia();
	// RayleighTaylor scenario is always a binary mixture of positive and negative particles.
	//	if (_binaryMixture) {
	//		_components[1].updateMassInertia();
	//	}
	//	std::cout << "from RayleighTaylor::readPhaseSpace " << std::flush;
	int id = 1;
	int numOfAddedMolecule = 0;

	srand((unsigned) time(NULL));
	for(int idx = 0; idx < _N; idx++){

		int componentType = 0;
		componentType = randdouble(0, 1.999999);

		double position_x;
		double position_y;
		double position_z;

		double &x = position_x;
		double &y = position_y;
		double &z = position_z;

		getRandomPosition(_L1,_L2,_L3,x,y,z,componentType);

		if (domainDecomp->procOwnsPos(x,y,z, domain)) {
			addMolecule(x, y, z, id, componentType,particleContainer);
			numOfAddedMolecule ++;
		}
		//std::cout << "x is" << x << std::flush;
		id++;
	}

	//std::cout << "numOfAddedMolecule" << numOfAddedMolecule << std::flush;

	removeMomentum(particleContainer, _components);

	unsigned long int globalNumMolecules = particleContainer->getNumberOfParticles();
	std::cout << "numberOfParticles is " << particleContainer->getNumberOfParticles()<< "\n"<<std::flush;
	domainDecomp->collCommInit(1);

	domainDecomp->collCommAppendUnsLong(globalNumMolecules);
	domainDecomp->collCommAllreduceSum();
	globalNumMolecules = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();

	domain->setglobalNumMolecules(globalNumMolecules);
	inputTimer.stop();
	_logger->info() << "Initial IO took:                 " << inputTimer.get_etime() << " sec" << endl;
	return id;
}

vector<ParameterCollection*> RayleighTaylorGenerator::getParameters() {
	vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("RayleighTaylorGenerator", "RayleighTaylorGenerator",
			"Parameters of RayleighTaylorGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithDoubleValue("SimulationBoxSize(X)", "SimulationBoxSize(X)",
					"SimulationBoxSize(X)", Parameter::LINE_EDIT,
					false, _L1 / MDGenerator::angstroem_2_atomicUnitLength));

	tab->addParameter(
			new ParameterWithDoubleValue("SimulationBoxSize(y)", "SimulationBoxSize(y)",
					"SimulationBoxSize(y)", Parameter::LINE_EDIT,
					false, _L2 / MDGenerator::angstroem_2_atomicUnitLength));

	tab->addParameter(
			new ParameterWithDoubleValue("SimulationBoxSize(Z)", "SimulationBoxSize(Z)",
					"SimulationBoxSize(Z)", Parameter::LINE_EDIT,
					false, _L3 / MDGenerator::angstroem_2_atomicUnitLength));

	tab->addParameter(
			new ParameterWithIntValue("NumOfParticles", "NumOfParticles",
					"NumOfParticles", Parameter::LINE_EDIT, false, _N));

	tab->addParameter(
			new ParameterWithDoubleValue("Temperature[K]", "Temperature[K]",
					"Temperature[K]", Parameter::LINE_EDIT, false, _T / MDGenerator::kelvin_2_mardyn));

	tab->addParameter(
			new ParameterWithDoubleValue("epsilon_A", "epsilon_A",
					"epsilon_A", Parameter::LINE_EDIT, false, _epsilon_A));

	tab->addParameter(
			new ParameterWithDoubleValue("epsilon_B", "epsilon_B",
					"epsilon_B", Parameter::LINE_EDIT, false, _epsilon_B));

	tab->addParameter(
			new ParameterWithDoubleValue("sigma_A", "sigma_A",
					"sigma_A", Parameter::LINE_EDIT, false, _sigma_A));

	tab->addParameter(
			new ParameterWithDoubleValue("sigma_B", "sigma_B",
					"sigma_B", Parameter::LINE_EDIT, false, _sigma_B));

	tab->addParameter(
			new ParameterWithDoubleValue("q_A", "q_A",
					"q_A", Parameter::LINE_EDIT, false, _q_A));
	tab->addParameter(
			new ParameterWithDoubleValue("q_B", "q_B",
					"q_B", Parameter::LINE_EDIT, false, _q_B));

	tab->addParameter(
			new ParameterWithDoubleValue("m_A", "m_A",
					"m_A", Parameter::LINE_EDIT, false, _m_A));

	tab->addParameter(
			new ParameterWithDoubleValue("m_B", "m_B",
					"m_B", Parameter::LINE_EDIT, false, _m_B));

	tab->addParameter(
			new ParameterWithDoubleValue("r_cut", "r_cut",
					"r_cut", Parameter::LINE_EDIT, false, _r_cut));

	tab->addParameter(
			new ParameterWithDoubleValue("delta_t", "delta_t",
					"delta_t", Parameter::LINE_EDIT, false, _delta_t));

	tab->addParameter(
			new ParameterWithDoubleValue("G", "G",
					"G", Parameter::LINE_EDIT, false, _G));

	tab->addParameter(
			new ParameterWithDoubleValue("h", "h",
					"h", Parameter::LINE_EDIT, false, _h));

	tab->addParameter(
			new ParameterWithDoubleValue("p_max", "p_max",
					"p_max", Parameter::LINE_EDIT, false, _p_max));

	tab->addParameter(
			new ParameterWithDoubleValue("skal", "skal",
					"skal", Parameter::LINE_EDIT, false, _skal));

	tab->addParameter(
			new ComponentParameters("component1", "component1", "Set up the parameters of component 1",
					_components[0]));

	tab->addParameter(
			new ComponentParameters("component2", "component2", "Set up the parameters of component 2",
					_components[1]));

	return parameters;
}

void RayleighTaylorGenerator::setParameter(Parameter* p) {

	string id = p->getNameId();

	if (id == "SimulationBoxSize(X)") {
		_L1 = static_cast<ParameterWithDoubleValue*> (p)->
				getValue() * MDGenerator::angstroem_2_atomicUnitLength;
		cout << "OneCenterLJRayleighTaylor: SimulationBoxSize(X): "
				<< _L1 / MDGenerator::angstroem_2_atomicUnitLength
				<< endl;

	} else if (id == "SimulationBoxSize(y)") {
		_L2 = static_cast<ParameterWithDoubleValue*> (p)->
				getValue() * MDGenerator::angstroem_2_atomicUnitLength;
		cout << "OneCenterLJRayleighTaylor: SimulationBoxSize(y): "
				<< _L2 / MDGenerator::angstroem_2_atomicUnitLength
				<< endl;

	} else if (id == "SimulationBoxSize(Z)") {
		_L3 = static_cast<ParameterWithDoubleValue*> (p)->
				getValue() * MDGenerator::angstroem_2_atomicUnitLength;
		cout << "OneCenterLJRayleighTaylor: SimulationBoxSize(Z): "
				<< _L3 / MDGenerator::angstroem_2_atomicUnitLength
				<< endl;

	} else if (id == "epsilon_A") {
		_epsilon_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: epsilon_A: " << _epsilon_A
				<< endl;

	} else if (id == "epsilon_B") {
		_epsilon_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: epsilon_B: " << _epsilon_B
				<< endl;

	} else if (id == "sigma_A") {
		_sigma_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: sigma_A: " << _sigma_A
				<< endl;

	} else if (id == "sigma_B") {
		_sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: sigma_B: " << _sigma_B
				<< endl;

	} else if (id == "m_A") {
		_m_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: m_A: " << _m_A
				<< endl;

	} else if (id == "m_B") {
		_m_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: m_B: " << _m_B
				<< endl;

	} else if (id == "NumOfParticles") {
		_N = static_cast<ParameterWithIntValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: NumOfParticles: " << _N
				<< endl;

	} else if (id == "r_cut") {
		_r_cut = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: r_cut: " << _r_cut
				<< endl;

	} else if (id == "delta_t") {
		_delta_t = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: delta_t: " << _delta_t
				<< endl;

	} else if (id == "numSphereSizes") {
		numSphereSizes = static_cast<ParameterWithIntValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: numSphereSizes: " << numSphereSizes
				<< endl;

	} else if (id == "Temperature[K]") {
		_T = static_cast<ParameterWithDoubleValue*> (p)->getValue()* MDGenerator::kelvin_2_mardyn;
		cout << "OneCenterLJRayleighTaylor: Temperature[K]: " << _T
				<< endl;

	} else if (id == "q_A") {
		_q_A = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: Temperature[K]: " << _q_A
				<< endl;

	} else if (id == "q_B") {
		_q_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: Temperature[K]: " << _q_B
				<< endl;

	} else if (id == "G") {
		_G = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: G: " << _G
				<< endl;

	} else if (id == "h") {
		_sigma_B = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor: h: " << _h
				<< endl;

	} else if (id == "p_max") {
		_p_max = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor:p_max: " << _p_max
				<< endl;

	} else if (id == "skal") {
		_skal = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJRayleighTaylor:skal: " << _skal
				<< endl;

	}  else if (id.find("component1") != std::string::npos) {
		std::string part = id.substr(11);
		ComponentParameters::setParameterValue(_components[0], p, part);

	}  else if (id.find("component2") != std::string::npos) {
		std::string part = id.substr(11);
		ComponentParameters::setParameterValue(_components[1], p, part);

	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);

	} else {
		std::cout << "UNKOWN Parameter: id = " << p->getNameId() << " value= " << p->getStringValue() << endl;
		exit(-1);
	}
}

void RayleighTaylorGenerator::addMolecule(
		double x, double y, double z, unsigned long id, int componentType, ParticleContainer* particleContainer) {
	vector<double> velocity = getRandomVelocity(_T);

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
		w[d] = w[d] * MDGenerator::fs_2_mardyn;
	}
	//End Copy

	Molecule m(id, componentType, x, y, z, // position
			velocity[0], -velocity[1], velocity[2], // velocity
			orientation[0], orientation[1], orientation[2], orientation[3],
			w[0], w[1], w[2], &_components);
	particleContainer->addParticle(m);
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
	L[0] = _L1;
	L[1] = _L2;
	L[2] = _L3;

	for (int i = 0; i < 3; i++) {
		//cout << "\nsimBox is " << simBoxLength[i]<< flush;
		if (L[i] < 2.0 * _configuration.getCutoffRadius()) {
			valid = false;
			_logger->error() << "Cutoff radius is too big (there would be only 1 cell in the domain!)" << endl;
			_logger->error() << "Cutoff radius=" << _configuration.getCutoffRadius()
																	<< " domain size=" << L[i] << endl;
		}
	}

	return valid;
}

