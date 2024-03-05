/*
 * MS2RSTGenerator.cpp
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt, eckhardw
 */

#include "MS2RSTGenerator.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Parameters/ParameterWithStringValue.h"
#include "Parameters/ParameterWithLongIntValue.h"
#include "Parameters/ParameterWithBool.h"
#include "common/MardynConfigurationParameters.h"
#include "common/PrincipalAxisTransform.h"
#include "molecules/Molecule.h"
#include "common/MS2RestartReader.h"
#include "Tokenize.h"
#include <cstring>

#ifndef MARDYN
extern "C" {

	Generator* create_generator() {
		return new MS2RSTGenerator();
	}

	void destruct_generator(Generator* generator) {
		delete generator;
	}
}
#endif

MS2RSTGenerator::MS2RSTGenerator() :
	MDGenerator("MS2RSTGenerator"), _molarDensity(0),
	_temperature(0), _simBoxLength(0), _ms2_to_angstroem(0),
	_filePath(""), _hasRotationalDOF(false), _numMolecules(0), _components(*(global_simulation->getEnsemble()->getComponents()))
{
	_components.resize(1);
	_components[0].addLJcenter(0, 0, 0, 1.0, 1.0, 1.0, 0.0, false);
	_components[0].setID(0);
}

std::vector<ParameterCollection*> MS2RSTGenerator::getParameters() {
	std::vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("MS2RSTParameters", "Parameters of MS2RSTGenerator",
			"Parameters of MS2RSTGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithDoubleValue("ms2ToAngstroem", "Factor MS2 unit length to Angstroem",
					"The conversion factor for unit length in the ms2 scenario to Angstroem", Parameter::LINE_EDIT,
						false, _ms2_to_angstroem));
	tab->addParameter(
			new ParameterWithStringValue("path", "Path to ms2 restart file",
					"Path to ms2 restart file", Parameter::LINE_EDIT,
					false, _filePath ));
	tab->addParameter(
			new ParameterWithDoubleValue("molarDensity", "Molar density [mol/l]",
					"molar density in mol/l in the scenario", Parameter::LINE_EDIT,
					false, _molarDensity));
	tab->addParameter(
			new ParameterWithLongIntValue("numMolecules", "Number of Molecules",
					"Number of molecules in the scenario that has to be read", Parameter::LINE_EDIT,
					true, _numMolecules));
	tab->addParameter(
			new ParameterWithDoubleValue("temperature", "Temperature [K]",
					"Temperature in the domain in Kelvin", Parameter::LINE_EDIT,
					false, _temperature / MDGenerator::kelvin_2_mardyn ));
	tab->addParameter(
			new ComponentParameters("component1", "component1",
					"Set up the parameters of component 1", _components[0]));
	tab->addParameter(
			new ParameterWithBool("rotDOF", "Rotational DOF",
					"Check this option if orientation and angular momentum have to be read in per molecule",
					Parameter::CHECKBOX, false, _hasRotationalDOF));

	return parameters;
}


void MS2RSTGenerator::setParameter(Parameter* p) {
	std::string id = p->getNameId();
	if (id == "numMolecules") {
		_numMolecules = static_cast<ParameterWithLongIntValue*> (p)->getValue();
		calculateSimulationBoxLength();
	} else if (id == "path") {
		_filePath = static_cast<ParameterWithStringValue*> (p)->getStringValue();
	} else if (id == "ms2ToAngstroem") {
		_ms2_to_angstroem = static_cast<ParameterWithDoubleValue*> (p)->getValue();
	} else if (id == "molarDensity") {
		_molarDensity = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		calculateSimulationBoxLength();
	} else if (id == "temperature") {
		_temperature = static_cast<ParameterWithDoubleValue*> (p)->getValue() * MDGenerator::kelvin_2_mardyn;
	} else if (id.find("component1") != std::string::npos) {
		std::string part = remainingSubString(".", id);
		ComponentParameters::setParameterValue(_components[0], p, part);
	} else if (id == "rotDOF") {
		_hasRotationalDOF = static_cast<ParameterWithBool*>(p)->getValue();
	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);
	}
}


void MS2RSTGenerator::calculateSimulationBoxLength() {
	// 1 mol/l = 0.6022 / nm^3 = 0.0006022 / Ang^3 = 0.089236726516 / a0^3
	double parts_per_a0 = _molarDensity * MDGenerator::molPerL_2_mardyn;
	double volume = _numMolecules / parts_per_a0;
	_simBoxLength = pow(volume, 1. / 3.);
}


void MS2RSTGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
	_logger->info() << "Reading PhaseSpaceHeader from MS2RSTGenerator..." << std::endl;

	//domain->setCurrentTime(0);
	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_temperature);
	domain->setGlobalLength(0, _simBoxLength);
	domain->setGlobalLength(1, _simBoxLength);
	domain->setGlobalLength(2, _simBoxLength);

	for (unsigned int i = 0; i < _components.size(); i++) {
		_components[i].updateMassInertia();
		if (_configuration.performPrincipalAxisTransformation()) {
			principalAxisTransform(_components[i]);
		}
	}
	domain->setepsilonRF(1e+10);
	_logger->info() << "Reading PhaseSpaceHeader from MS2RSTGenerator done." << std::endl;

    /* silence compiler warnings */
    (void) timestep;
}


unsigned long MS2RSTGenerator::readPhaseSpace(ParticleContainer* particleContainer,
		Domain* domain, DomainDecompBase* domainDecomp) {

	global_simulation->timers()->start("MS2RST_GENERATOR_INPUT");
	_logger->info() << "Reading phase space file (MS2RSTGenerator)." << std::endl;

	std::vector<bool> rotationDOF(1);
	rotationDOF[0] = _hasRotationalDOF;

	MS2RestartReader::MoleculeData* ms2mols = MS2RestartReader::readMS2RestartFile(
			_filePath, 1, _numMolecules, rotationDOF);

	for (unsigned int i = 0; i < _numMolecules; i++) {
		for (int j = 0; j < 3; j++) {
			ms2mols[i].x[j] = ms2mols[i].x[j] * _ms2_to_angstroem * angstroem_2_atomicUnitLength;
			if (ms2mols[i].x[j] < 0 || ms2mols[i].x[j] > _simBoxLength) {
				std::cout << "Error: molecule out of box: " << std::endl;
				ms2mols[i].print(cout);
			}
		}
		addMolecule(ms2mols[i], particleContainer);
	}

	// todo: temperature!

	delete ms2mols;

	removeMomentum(particleContainer, _components);

	thermostat(particleContainer);

	domain->evaluateRho(particleContainer->getNumberOfParticles(), domainDecomp);
	_logger->info() << "Calculated Rho=" << domain->getglobalRho() << std::endl;
	global_simulation->timers()->start("MS2RST_GENERATOR_INPUT");
	global_simulation->timers()->setOutputString("MS2RST_GENERATOR_INPUT", "Initial IO took:                 ");
	_logger->info() << "Initial IO took:                 " << global_simulation->timers()->getTime("MS2RST_GENERATOR_INPUT") << " sec" << std::endl;
	return _numMolecules;
}

void MS2RSTGenerator::addMolecule(MS2RestartReader::MoleculeData& ms2mol, ParticleContainer* particleContainer) {
	Molecule m(ms2mol.id, &_components[ms2mol.cid], ms2mol.x[0], ms2mol.x[1], ms2mol.x[2], // position
			ms2mol.v[0], ms2mol.v[1], ms2mol.v[2], // velocity
			ms2mol.q[0], ms2mol.q[1], ms2mol.q[2], ms2mol.q[3],
			ms2mol.d[0], ms2mol.d[1], ms2mol.d[2] );
	if (particleContainer->isInBoundingBox(m.r_arr().data())) {
		bool inChecked = true;
		particleContainer->addParticle(m, inChecked);
	}
}


void MS2RSTGenerator::thermostat(ParticleContainer* container) {
	double summv2 = 0;
	double sumIw2 = 0;
	unsigned long numberOfRotationalDOF = 0;
	for (auto tM = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
		tM->calculate_mv2_Iw2(summv2, sumIw2);
		numberOfRotationalDOF += _components[ tM->componentid() ].getRotationalDegreesOfFreedom();
	}

	double currentTemperature = (summv2 + sumIw2) / (double)(3 * _numMolecules + numberOfRotationalDOF);

	double betaTrans = pow( 3.0 * _numMolecules * _temperature / summv2, 0.5);
	double betaRot = 1.0;
	if (_hasRotationalDOF) {
		betaRot = pow( numberOfRotationalDOF * _temperature / sumIw2, 0.5);
	}

	for (auto tM = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
		tM->scale_v(betaTrans);
		tM->scale_D(betaRot);
	}

	_logger->info() << "Current Temperature: " << currentTemperature << ", target Temperature: " << _temperature << std::endl;
	_logger->info() << " bTrans=" << betaTrans << ", bRot=" << betaRot << " #RotDOF=" << numberOfRotationalDOF << std::endl;

}


bool MS2RSTGenerator::validateParameters() {
	bool valid = true;

	if (_configuration.getScenarioName() == "") {
		valid = false;
		_logger->error() << "ScenarioName not set!" << std::endl;
	}

	if (_configuration.getOutputFormat() == MardynConfiguration::XML) {
		valid = false;
		_logger->error() << "OutputFormat XML not yet supported!" << std::endl;
	}

	if (_simBoxLength < 2. * _configuration.getCutoffRadius()) {
		valid = false;
		_logger->error() << "Cutoff radius is too big (there would be only 1 cell in the domain!)" << std::endl;
		_logger->error() << "Cutoff radius=" << _configuration.getCutoffRadius()
							<< " domain size=" << _simBoxLength << std::endl;
	}
	return valid;
}


