/*
 * AqueousNaClGenerator.cpp
 *
 *  Created on: Oct 14, 2012
 *      Author: eckhardw
 */

#include "AqueousNaClGenerator.h"
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
		return new AqueousNaClGenerator();
	}

	void destruct_generator(Generator* generator) {
		delete generator;
	}
}
#endif


AqueousNaClGenerator::AqueousNaClGenerator() :
		MDGenerator("AqueousNaClGenerator"), _numMolecules(4), _temperature(298. / 315774.5), _molarDensity(57.556), _components(
				*(global_simulation->getEnsemble()->getComponents())) {

	_components.resize(3);
	// set up water
	_components[0].addLJcenter(0, 0.123891518, 0, 0.016, 0.000246810271, 5.95953717, 0.0, false);
	_components[0].addCharge(0, -0.159567514, 0, 0, -1.04);
	_components[0].addCharge(1.43042933, -0.983266012, 0, 0.001008, 0.52);
	_components[0].addCharge(-1.43042933, -0.983266012, 0, 0.001008, 0.52);
	_components[0].setID(0);

	// Na+
	_components[1].addLJcenter(0, 0, 0, 0.0229897, 0.062 * MDGenerator::kelvin_2_mardyn,
			2.58 * MDGenerator::angstroem_2_atomicUnitLength, 0.0, false);
	_components[1].addCharge(0, 0, 0, 0, 1.0);
	_components[1].setID(1);
	// Cl-
	_components[2].addLJcenter(0, 0, 0, 0.0354532, 0.446 * MDGenerator::kelvin_2_mardyn,
				4.45 * MDGenerator::angstroem_2_atomicUnitLength, 0.0, false);
	_components[2].addCharge(0, 0, 0, 0, -1.0);
	_components[2].setID(2);

	calculateSimulationBoxLength();
}

std::vector<ParameterCollection*> AqueousNaClGenerator::getParameters() {
	std::vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("AqueousNaClParameters", "Parameters of EqvGridGenerator",
			"Parameters of AqueousNaClGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithLongIntValue("numMolecules", "Number of Molecules",
					"Total number of Molecules", Parameter::LINE_EDIT,
					true, _numMolecules));
	tab->addParameter(
			new ParameterWithDoubleValue("temperature", "Temperature [K]",
					"Temperature in the domain in Kelvin", Parameter::LINE_EDIT,
					false, _temperature / MDGenerator::kelvin_2_mardyn ));

	return parameters;
}


void AqueousNaClGenerator::setParameter(Parameter* p) {
	std::string id = p->getNameId();
	if (id == "numMolecules") {
		_numMolecules = static_cast<ParameterWithLongIntValue*> (p)->getValue();
		calculateSimulationBoxLength();
	} else if (id == "temperature") {
		_temperature = static_cast<ParameterWithDoubleValue*> (p)->getValue() * MDGenerator::kelvin_2_mardyn;
	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);
	}
}


void AqueousNaClGenerator::calculateSimulationBoxLength() {
	// taken from the Tironi paper
	double volumePerMolecule = 66400 * pow(MDGenerator::angstroem_2_atomicUnitLength, 3) / 2207;
	int numMoleculesPerDimension = pow((double) _numMolecules, 1./3.);
	_numMolecules = pow((double) numMoleculesPerDimension, 3);
	double volume = _numMolecules * volumePerMolecule;
	_simBoxLength = pow(volume, 1./3.);
}


void AqueousNaClGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
	_logger->info() << "Reading PhaseSpaceHeader from AqueousNaClGenerator..." << std::endl;
	//domain->setCurrentTime(0);

	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_temperature);
	domain->setGlobalLength(0, _simBoxLength);
	domain->setGlobalLength(1, _simBoxLength);
	domain->setGlobalLength(2, _simBoxLength);

	for (unsigned int i = 0; i < _components.size(); i++) {
		principalAxisTransform(_components[i]);
	}
	domain->setepsilonRF(1e+10);
	_logger->info() << "Reading PhaseSpaceHeader from AqueousNaClGenerator done." << std::endl;

    /* silence compiler warnings */
    (void) timestep;

	global_simulation->getEnsemble()->setComponentLookUpIDs();
}


unsigned long AqueousNaClGenerator::readPhaseSpace(ParticleContainer* particleContainer,
		Domain* domain, DomainDecompBase* domainDecomp) {

	global_simulation->timers()->start("AQUEOUS_NA_CL_GENERATOR_INPUT");
	_logger->info() << "Reading phase space file (AqueousNaClGenerator)." << std::endl;

	int numMoleculesPerDimension = round(pow((double) _numMolecules, 1./3.));
	int numIons = round(0.01812415 * (double) _numMolecules);
	_logger->info() << "Generating " << numIons << " of Na+ and Cl-" << std::endl;
	Ion* ions = new Ion[2 * numIons];
	for (int i = 0; i < numIons; i++) {
		ions[i].cid = 1;
		do {
			for (int j = 0; j < 3; j++) {
				ions[i].position[j] = (double) rand() / ((double) RAND_MAX) * (numMoleculesPerDimension-1);
			}
			std::cout << "Testing " << ions[i].position[0] << "," << ions[i].position[1] << "," << ions[i].position[2] << std::endl;
		} while (isNearIon(ions, i));
		std::cout << "Generate Na+ at " << ions[i].position[0] << "," << ions[i].position[1] << "," << ions[i].position[2] << ";" << std::endl;
	}
	for (int i = numIons; i < 2*numIons; i++) {
		ions[i].cid = 2;
		do {
			for (int j = 0; j < 3; j++) {
				ions[i].position[j] = (double) rand() / ((double) RAND_MAX) * (numMoleculesPerDimension-1);
			}
			std::cout << "Testing " << ions[i].position[0] << "," << ions[i].position[1] << "," << ions[i].position[2] << std::endl;
		} while (isNearIon(ions, i));
		std::cout << "Generate Na+ at " << ions[i].position[0] << "," << ions[i].position[1] << "," << ions[i].position[2] << ";" << std::endl;
	}

	unsigned long int id = 1;
	double spacing = _simBoxLength / numMoleculesPerDimension;
	_logger->info() << "SimBoxLength=" << _simBoxLength << " spacing=" << spacing << std::endl;
	double origin = spacing / 2.; // origin of the first DrawableMolecule

	// only for console output
	double percentage = 1.0 / numMoleculesPerDimension * 100.0;
	int percentageRead = 0;

	for (int i = 0; i < numMoleculesPerDimension; i++) {
		for (int j = 0; j < numMoleculesPerDimension; j++) {
			for (int k = 0; k < numMoleculesPerDimension; k++) {

				double x = origin + i * spacing;
				double y = origin + j * spacing;
				double z = origin + k * spacing;

				int cid = getCID(ions, 2 * numIons, i, j, k);
				if (domainDecomp->procOwnsPos(x,y,z, domain)) {
					addMolecule(x, y, z, id, cid, particleContainer);
				}
				// increment id in any case, because this particle will probably
				// be added by some other process
				id++;
			}
		}

		percentageRead = i * percentage;
		_logger->info() << "Finished reading molecules: " << (percentageRead) << "%\r" << std::flush;
	}

	delete[] ions;

	removeMomentum(particleContainer, _components);
	domain->evaluateRho(particleContainer->getNumberOfParticles(), domainDecomp);
	_logger->info() << "Calculated Rho=" << domain->getglobalRho() << std::endl;
	global_simulation->timers()->stop("AQUEOUS_NA_CL_GENERATOR_INPUT");
	global_simulation->timers()->setOutputString("AQUEOUS_NA_CL_GENERATOR_INPUT", "Initial IO took:                 ");
	_logger->info() << "Initial IO took:                 " << global_simulation->timers()->getTime("AQUEOUS_NA_CL_GENERATOR_INPUT") << " sec" << std::endl;
	return id;
}

void AqueousNaClGenerator::addMolecule(double x, double y, double z, unsigned long id, int cid, ParticleContainer* particleContainer) {
	std::cout << "Add molecule at " << x << ", " << y << ", " << z << " with cid=" << cid << std::endl;

	std::vector<double> velocity = getRandomVelocity(_temperature);

	//double orientation[4] = {1, 0, 0, 0}; // default: in the xy plane
	// rotate by 30° along the vector (1/1/0), i.e. the angle bisector of x and y axis
	// o = cos 30° + (1 1 0) * sin 15°
	double orientation[4];
	getOrientation(15, 10, orientation);

	double I[3] = {0.,0.,0.};
	I[0] = _components[cid].I11();
	I[1] = _components[cid].I22();
	I[2] = _components[cid].I33();
	/*****  Copied from animake - initialize anular velocity *****/
	double w[3];
	for(int d=0; d < 3; d++) {
		w[d] = (I[d] == 0)? 0.0: ((randdouble(0,1) > 0.5)? 1: -1) *
				sqrt(2.0* randdouble(0,1)* _temperature / I[d]);
		w[d] = w[d] * MDGenerator::fs_2_mardyn;
	}
	/************************** End Copy **************************/

	Molecule m(id, &(_components[cid]), x, y, z, // position
			velocity[0], -velocity[1], velocity[2], // velocity
			orientation[0], orientation[1], orientation[2], orientation[3],
			w[0], w[1], w[2] );
	if (particleContainer->isInBoundingBox(m.r_arr().data())) {
		bool inChecked = true;
		particleContainer->addParticle(m, inChecked);
	}
}

bool AqueousNaClGenerator::isNearIon(Ion* ions, int ni) const {
	for (int i = 0; i < ni; i++) {
		int* curPos = ions[ni].position;
		int* pos = ions[i].position;
		if (abs(pos[0] - curPos[0]) < 2
				&& abs(pos[1] - curPos[1]) < 2
				&& abs(pos[2] - curPos[2]) < 2) return true;
	}

	return false;
}

int AqueousNaClGenerator::getCID(Ion ions[], int numIons, int x, int y, int z) const {
	for (int i = 0; i < numIons; i++) {
		if (ions[i].position[0] == x && ions[i].position[1] == y && ions[i].position[2] == z) {
			return ions[i].cid;
		}
	}

	return 0;
}

bool AqueousNaClGenerator::validateParameters() {
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
