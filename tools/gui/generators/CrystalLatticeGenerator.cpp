/*
 * CrystalLatticeGenerator.cpp
 *
 *  Created on: Oct 13, 2012
 *      Author: Wolfgang Eckhardt
 */

#include "CrystalLatticeGenerator.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Parameters/ParameterWithLongIntValue.h"
#include "Parameters/ParameterWithBool.h"
#include "common/MardynConfigurationParameters.h"
#include "common/PrincipalAxisTransform.h"
#include "molecules/Molecule.h"
#include "Tokenize.h"
#include <cstring>

#ifndef MARDYN
extern "C" {

	Generator* create_generator() {
		return new CrystalLatticeGenerator();
	}

	void destruct_generator(Generator* generator) {
		delete generator;
	}
}
#endif

CrystalLatticeGenerator::CrystalLatticeGenerator() :
	MDGenerator("CrystalLatticeGenerator"), _numMoleculesPerDim(4), _h(3.0), _charge(1.0), _components(*(global_simulation->getEnsemble()->getComponents())) {
	_components.resize(2);
	_components[0].addCharge(0, 0, 0, 1.0, _charge);
	_components[1].addCharge(0, 0, 0, 1.0, - _charge);
	_configuration.setNVE(true);
	calculateSimulationBoxLength();
}

std::vector<ParameterCollection*> CrystalLatticeGenerator::getParameters() {
	std::vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("CrystalLatticeParameters", "Parameters of CrystalLatticeGenerator",
			"Parameters of CrystalLatticeGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithDoubleValue("h", "h [A]",
					"The grid spacing h in Angstroem", Parameter::LINE_EDIT,
					false, _h / MDGenerator::angstroem_2_atomicUnitLength));
	tab->addParameter(
			new ParameterWithLongIntValue("numMolecules", "Number of Molecules per dimension",
					"Number of Molecules per dimension (must be an even number!)", Parameter::LINE_EDIT,
					false, _numMoleculesPerDim));
	tab->addParameter(
			new ParameterWithDoubleValue("charge", "Charge [e]",
					"Charge of the ions in elementary charges [e]", Parameter::LINE_EDIT,
					false, _charge ));

	return parameters;
}


void CrystalLatticeGenerator::setParameter(Parameter* p) {
	std::string id = p->getNameId();
	if (id == "numMolecules") {
		_numMoleculesPerDim = static_cast<ParameterWithLongIntValue*> (p)->getValue();
		calculateSimulationBoxLength();
	} else if (id == "h") {
		_h = static_cast<ParameterWithDoubleValue*> (p)->getValue() * MDGenerator::angstroem_2_atomicUnitLength;
		calculateSimulationBoxLength();
	} else if (id == "charge") {
		_charge = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		_components[0].charge(0).setQ(_charge);
		_components[1].charge(0).setQ(-_charge);
	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);
	}
}


void CrystalLatticeGenerator::calculateSimulationBoxLength() {
	// 1 mol/l = 0.6022 / nm^3 = 0.0006022 / Ang^3 = 0.089236726516 / a0^3
	double length = _numMoleculesPerDim * _h;
	_simBoxLength = length;
}


void CrystalLatticeGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
	_logger->info() << "Reading PhaseSpaceHeader from CubicGridGenerator..." << std::endl;
	//domain->setCurrentTime(0);

	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(0);
	domain->setGlobalLength(0, _simBoxLength);
	domain->setGlobalLength(1, _simBoxLength);
	domain->setGlobalLength(2, _simBoxLength);

	domain->setepsilonRF(1e+10);
	_logger->info() << "Reading PhaseSpaceHeader from CubicGridGenerator done." << std::endl;

    /* silence compiler warnings */
    (void) timestep;

	global_simulation->getEnsemble()->setComponentLookUpIDs();
}


unsigned long CrystalLatticeGenerator::readPhaseSpace(ParticleContainer* particleContainer,
		Domain* domain, DomainDecompBase* domainDecomp) {

	global_simulation->timers()->start("CRYSTAL_LATTICE_GENERATOR_INPUT");
	_logger->info() << "Reading phase space file (CubicGridGenerator)." << std::endl;

	unsigned long int id = 1;
	double spacing = _simBoxLength / _numMoleculesPerDim;
	double origin = spacing / 2.; // origin of the first DrawableMolecule

	// only for console output
	double percentage = 1.0 / (_numMoleculesPerDim) * 100.0;
	int percentageRead = 0;

	for (unsigned i = 0; i < _numMoleculesPerDim; i+=2) {
		for (unsigned j = 0; j < _numMoleculesPerDim; j+=2) {
			for (unsigned k = 0; k < _numMoleculesPerDim; k+=2) {

				double x = origin + i * spacing;
				double y = origin + j * spacing;
				double z = origin + k * spacing;
				if (domainDecomp->procOwnsPos(x,y,z, domain)) {
					addMolecule(x, y, z, id, 0, particleContainer);
				}
				if (domainDecomp->procOwnsPos(x + spacing ,y,z, domain)) {
					addMolecule(x + spacing, y, z, id, 1, particleContainer);
				}
				if (domainDecomp->procOwnsPos(x,y + spacing,z, domain)) {
					addMolecule(x, y + spacing, z, id, 1, particleContainer);
				}
				if (domainDecomp->procOwnsPos(x + spacing ,y + spacing,z, domain)) {
					addMolecule(x + spacing, y + spacing, z, id, 0, particleContainer);
				}

				if (domainDecomp->procOwnsPos(x,y,z + spacing, domain)) {
					addMolecule(x, y, z + spacing, id, 1, particleContainer);
				}
				if (domainDecomp->procOwnsPos(x + spacing ,y,z + spacing, domain)) {
					addMolecule(x + spacing, y, z + spacing, id, 0, particleContainer);
				}
				if (domainDecomp->procOwnsPos(x,y + spacing,z + spacing, domain)) {
					addMolecule(x, y + spacing, z + spacing, id, 0, particleContainer);
				}
				if (domainDecomp->procOwnsPos(x + spacing ,y + spacing,z + spacing, domain)) {
					addMolecule(x + spacing, y + spacing, z + spacing, id, 1, particleContainer);
				}
				// increment id in any case, because this particle will probably
				// be added by some other process
				id+=4;
			}
		}

		percentageRead = i * percentage;
		_logger->info() << "Finished reading molecules: " << (percentageRead) << "%\r" << std::flush;
	}

	domain->evaluateRho(particleContainer->getNumberOfParticles(), domainDecomp);
	_logger->info() << "Calculated Rho=" << domain->getglobalRho() << std::endl;
	global_simulation->timers()->stop("CRYSTAL_LATTICE_GENERATOR_INPUT");
	global_simulation->timers()->setOutputString("CRYSTAL_LATTICE_GENERATOR_INPUT", "Initial IO took:                 ");
	_logger->info() << "Initial IO took:                 " << global_simulation->timers()->getTime("CRYSTAL_LATTICE_GENERATOR_INPUT") << " sec" << std::endl;
	return id;
}

void CrystalLatticeGenerator::addMolecule(double x, double y, double z, unsigned long id, unsigned cid, ParticleContainer* particleContainer) {
	Molecule m(id, &_components[cid], x, y, z, // position
			0.0, 0.0, 0.0, // velocity
			1.0, 0.0, 0.0, 0.0, // orientation
			0.0, 0.0, 0.0);
	if (particleContainer->isInBoundingBox(m.r_arr().data())) {
		bool inChecked = true;
		particleContainer->addParticle(m, inChecked);
	}
}


bool CrystalLatticeGenerator::validateParameters() {
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

	if (_numMoleculesPerDim % 2) {
		valid = false;
		_logger->error() << "Number of molecules per dimension must be even!" << std::endl;
	}
	return valid;
}


