/*
 * MDGenerator.cpp
 *
 * @Date: 08.06.2011
 * @Author: eckhardw
 */

#include "MDGenerator.h"

#include <memory>

#include "parallel/DomainDecompBase.h"
#include "io/CheckpointWriter.h"
#include "particleContainer/LinkedCells.h"
#include "Domain.h"

#include "common/MardynConfigLegacyWriter.h"

#ifndef MARDYN
#include "common/DrawableMolecule.h"
#include "QObjects/ScenarioGenerator.h"
#include "utils/mardyn_assert.h"
#endif


double const MDGenerator::angstroem_2_atomicUnitLength = 1.889726878;

double const MDGenerator::unitMass_2_mardyn = 0.001; // Mardyn calculates with 1/1000 u as base unit.

double const MDGenerator::debye_2_mardyn = 0.393429724;

double const MDGenerator::buckingham_2_mardyn = 0.743472313;

double const MDGenerator::fs_2_mardyn = 0.030619994;

const double MDGenerator::molPerL_2_mardyn = 8.923894e-5;

const double MDGenerator::kelvin_2_mardyn = 1.0 / 315774.5;

const double MDGenerator::unitCharge_2_mardyn = 1./96485.3;

const double MDGenerator::abogadro_constant = 6.02214078e23;

const double MDGenerator::boltzmann_constant_kB = 1.38065e-23;

MDGenerator::MDGenerator(std::string name) :
Generator(name) {
	// initialize monolithic Mardyn's global_log and silence it...
#ifndef MARDYN
	// if mardyn is not the main program, silence it's logger;
	Log::global_log = std::make_shared<Log::Logger>(Log::Warning, &(std::cout));
	_logger = std::make_shared<Log::Logger>(Log::Debug, &ScenarioGeneratorApplication::getInstance()->getTextMessageStream());
#else
	_logger = std::make_shared<Log::Logger>(Log::Debug, &(std::cout));
#endif
}

void MDGenerator::setLogger(std::shared_ptr<Log::Logger> logger) {
	_logger = logger;
}

void MDGenerator::createSampleObject() const {
#ifndef MARDYN
	std::vector<Component> components;
	components.resize(1);
	components[0].addLJcenter(0,0,0,1,1,1,1,false);
	Molecule m(1, &(components[0]), 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
	ScenarioGeneratorApplication::getInstance()->addObject(
			new DrawableMolecule(m, 1));
#endif
}

void MDGenerator::generatePreview() {
	srand(1);
#ifndef MARDYN
	int rank = 0;
	Domain domain(rank);
	DomainDecompBase domainDecomposition;

	double bBoxMin[3] = { 0,0,0};
	double bBoxMax[3] = { 0,0,0};
	double cutoffRadius = 3.0;

	readPhaseSpaceHeader(&domain, 0);
	bBoxMax[0] = domain.getGlobalLength(0);
	bBoxMax[1] = domain.getGlobalLength(1);
	bBoxMax[2] = domain.getGlobalLength(2);

	_logger->info() << "MDGenerator: bounding box of domain is [" <<
			bBoxMin[0] << "," << bBoxMin[1] << "," << bBoxMin[2] << "] to ["
			<< bBoxMax[0] << "," << bBoxMax[1] << "," << bBoxMax[2] << "]" << std::endl;
	_logger->info() << "MDGenerator: temperature=" << domain.getTargetTemperature(0) << std::endl;

	LinkedCells container(bBoxMin, bBoxMax, cutoffRadius);

	readPhaseSpace(&container, &domain, &domainDecomposition);
	_logger->info() << "MDGenerator: " << container.getNumberOfParticles() << " particles were created." << std::endl;

	auto molecule = container.iterator(ParticleIterator::ALL_CELLS);
	while (molecule.isValid()) {
		ScenarioGeneratorApplication::getInstance()->addObject(new DrawableMolecule(*molecule, global_simulation->getEnsemble()->getComponents()->size()-1));
		++molecule;
	}
#endif
}


void MDGenerator::generateOutput(const std::string& directory) {

	srand(1);
	std::cout << "MDGenerator::generateOutput called!" << std::endl;

	if (_configuration.getOutputFormat() == MardynConfiguration::LEGACY) {
		MardynConfigLegacyWriter::writeConfigFile(directory, _configuration.getScenarioName() + ".cfg", _configuration);
	} else if (_configuration.getOutputFormat() == MardynConfiguration::XML) {
		_logger->error() << "XML Output not yet supported!" << std::endl;
		_logger->error() << "Generating nothing!" << std::endl;
		return;
	} else {
		_logger->error() << "Invalid File format for Output!" << _configuration.getOutputFormat() << std::endl;
	}

	std::cout << "MDGenerator::config file written!" << std::endl;
	int rank = 0;
	Domain domain(rank);
	DomainDecompBase domainDecomposition;

#ifndef MARDYN
	global_simulation = new Simulation();
	global_simulation->setcutoffRadius(3.0);
#endif

	double bBoxMin[3] = { 0,0,0};
	double bBoxMax[3] = { 0,0,0};
	double cutoffRadius = 3.0;

	std::cout << "MDGenerator::generateOutput before read phasespace header!" << std::endl;
	readPhaseSpaceHeader(&domain, 0);
	std::cout << "MDGenerator::generateOutput after read phasespace header!" << std::endl;
	bBoxMax[0] = domain.getGlobalLength(0);
	bBoxMax[1] = domain.getGlobalLength(1);
	bBoxMax[2] = domain.getGlobalLength(2);
	LinkedCells container(bBoxMin, bBoxMax, cutoffRadius);
	std::cout << "MDGenerator::generateOutput before read phasespace!" << std::endl;
	readPhaseSpace(&container, &domain, &domainDecomposition);
	std::cout << "MDGenerator::generateOutput read phasespace done!" << std::endl;
	domain.setglobalNumMolecules(container.getNumberOfParticles());
	std::cout << "NumMolecules in Container: " << container.getNumberOfParticles() << std::endl;

	std::string destination = directory + "/" + _configuration.getScenarioName() + ".inp";
	_logger->info() << "Writing output to: " << destination << std::endl;
	domain.writeCheckpoint(destination, &container, &domainDecomposition, 0.);

#ifndef MARDYN
	delete global_simulation;
#endif
}

std::vector<double> MDGenerator::getRandomVelocity(double temperature) const {
	std::vector<double> v_;
	v_.resize(3);

	// Velocity
	for (int dim = 0; dim < 3; dim++) {
		v_[dim] = randdouble(-0.5, 0.5);
	}
	double dotprod_v = 0;
	for (unsigned int i = 0; i < v_.size(); i++) {
		dotprod_v += v_[i] * v_[i];
	}
	// Velocity Correction
	double vCorr = sqrt(3.0 * temperature / dotprod_v);
	for (unsigned int i = 0; i < v_.size(); i++) {
		v_[i] *= vCorr;
	}

	return v_;
}


void MDGenerator::getOrientation(int base, int delta, double orientation[4]) {
	double offset = randdouble(-delta / 2., delta / 2.) / 180. * M_PI;
	double rad = base / 180. * M_PI;
	double angle = rad + offset;

	double cosinePart = cos(angle);
	double sinePart = sin(angle);

	double length = sqrt(cosinePart * cosinePart + 2 * (sinePart * sinePart));
	orientation[0] = cosinePart / length;
	orientation[1] = sinePart / length;
	orientation[2] = sinePart / length;
	orientation[3] = 0;
}

bool MDGenerator::isInsideDomain(Domain* domain, double position[3]) {
	for (int i = 0; i < 3; i++) {
		if (position[i] < 0 || position[i] > domain->getGlobalLength(i)) {
			return false;
		}
	}

	return true;
}

void MDGenerator::removeMomentum(ParticleContainer* particleContainer, const std::vector<Component>& components) {
	double mass = 0.;
	double mass_sum = 0.;
	double momentum_sum[3] = { 0., 0., 0. };

	auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
	while(molecule.isValid()){
		mass = components[molecule->componentid()].m();
		mass_sum = mass_sum + mass;
		momentum_sum[0] = momentum_sum[0] + mass * molecule->v(0);
		momentum_sum[1] = momentum_sum[1] + mass * molecule->v(1);
		momentum_sum[2] = momentum_sum[2] + mass * molecule->v(2);
		++molecule;
	}

	double momentum_sub0 = momentum_sum[0] / mass_sum;
	double momentum_sub1 = momentum_sum[1] / mass_sum;
	double momentum_sub2 = momentum_sum[2] / mass_sum;

	molecule = particleContainer->iterator(ParticleIterator::ALL_CELLS);
	while (molecule.isValid()) {
		molecule->vsub(momentum_sub0, momentum_sub1, momentum_sub2);
		++molecule;
	}

	//test
	momentum_sum[0] = 0.;
	momentum_sum[1] = 0.;
	momentum_sum[2] = 0.;

	molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
	while (molecule.isValid()) {
		mass = components[molecule->componentid()].m();
		mass_sum = mass_sum + mass;
		momentum_sum[0] = momentum_sum[0] + mass * molecule->v(0);
		momentum_sum[1] = momentum_sum[1] + mass * molecule->v(1);
		momentum_sum[2] = momentum_sum[2] + mass * molecule->v(2);
		++molecule;
	}

	//printf("momentum_sum[0] from removeMomentum is %lf\n", momentum_sum[0]);
	//printf("momentum_sum[1] from removeMomentum is %lf\n", momentum_sum[1]);
	//printf("momentum_sum[2] from removeMomentum is %lf\n", momentum_sum[2]);
}

