/*
 * MDGenerator.cpp
 *
 * @Date: 08.06.2011
 * @Author: eckhardw
 */

#include "MDGenerator.h"

#include "ensemble/GrandCanonical.h"
#include "parallel/DomainDecompDummy.h"
#include "io/CheckpointWriter.h"
#include "ensemble/PressureGradient.h"
#include "particleContainer/LinkedCells.h"
#include "Domain.h"

#include "common/MardynConfigLegacyWriter.h"
#include "common/DrawableMolecule.h"
#include "QObjects/ScenarioGenerator.h"

MDGenerator::MDGenerator(std::string name) :
Generator(name), _deleteLogger(true) {
	// initialize monolithic Mardyn's global_log and silence it...
#ifndef MARDYN
	// if mardyn is not the main program, silence it's logger;
	Log::global_log = new Log::Logger(Log::Warning, &(std::cout));
	_logger = new Log::Logger(Log::Debug, &ScenarioGeneratorApplication::getInstance()->getTextMessageStream());
#else
	_logger = new Log::Logger(Log::Debug, &(std::cout));
#endif
}

MDGenerator::~MDGenerator() {
	if (_deleteLogger) {
		delete _logger;
	}
}

void MDGenerator::setLogger(Log::Logger* logger) {
	if (_logger != NULL && _deleteLogger) {
		delete logger;
	}

	_logger = logger;
	_deleteLogger = false;
}

const Object* MDGenerator::getSampleObject() const {
	return new DrawableMolecule();
}

void MDGenerator::generatePreview() {
	int rank = 0;
	PressureGradient gradient(rank);
	Domain domain(rank, &gradient);
	DomainDecompDummy domainDecomposition;
	list<ChemicalPotential> lmu;

	double bBoxMin[3] = { 0,0,0};
	double bBoxMax[3] = { 20,20, 20};
	double cutoffRadius = 3.0;
	double LJCutoffRadius = 3.0;
	double tersoffCutoffRadius = 3.0;
	double cellsInCutoffRadius = 1;
	ParticlePairsHandler* partPairsHandler = NULL;

	readPhaseSpaceHeader(&domain, 0);
	bBoxMax[0] = domain.getGlobalLength(0);
	bBoxMax[1] = domain.getGlobalLength(1);
	bBoxMax[2] = domain.getGlobalLength(2);

	_logger->info() << "MDGenerator: bounding box of domain is [" <<
			bBoxMin[0] << "," << bBoxMin[1] << "," << bBoxMin[2] << "] to ["
			<< bBoxMax[0] << "," << bBoxMax[1] << "," << bBoxMax[2] << "]" << endl;

	_logger->info() << "MDGenerator: temperature=" << domain.getGlobalCurrentTemperature() << endl;

	LinkedCells container(bBoxMin, bBoxMax, cutoffRadius,
			LJCutoffRadius, tersoffCutoffRadius, cellsInCutoffRadius, partPairsHandler);

	readPhaseSpace(&container, &lmu, &domain, &domainDecomposition);
	_logger->info() << "MDGenerator: " << container.getNumberOfParticles() << " particles were created." << endl;

	Molecule* molecule = container.begin();
	while (molecule != container.end()) {
		ScenarioGeneratorApplication::getInstance()->addObject(new DrawableMolecule(*molecule));
		molecule = container.next();
	}
}


void MDGenerator::generateOutput(const std::string& directory) {

	if (_configuration.getOutputFormat() == MardynConfiguration::LEGACY) {
		MardynConfigLegacyWriter::writeConfigFile(directory, _configuration.getScenarioName() + ".cfg", _configuration);
	} else if (_configuration.getOutputFormat() == MardynConfiguration::XML) {
		_logger->error() << "XML Output not yet supported!" << endl;
		_logger->error() << "Generating nothing!" << endl;
		return;
	} else {
		_logger->error() << "Invalid File format for Output!" << _configuration.getOutputFormat() << endl;
	}

	int rank = 0;
	PressureGradient gradient(rank);
	Domain domain(rank, &gradient);
	DomainDecompDummy domainDecomposition;
	list<ChemicalPotential> lmu;

	double bBoxMin[3] = { 0,0,0};
	double bBoxMax[3] = { 10,10,10};
	double cutoffRadius = 3.0;
	double LJCutoffRadius = 3.0;
	double tersoffCutoffRadius = 3.0;
	double cellsInCutoffRadius = 1;
	ParticlePairsHandler* partPairsHandler = NULL;
	LinkedCells container(bBoxMin, bBoxMax, cutoffRadius,
			LJCutoffRadius, tersoffCutoffRadius, cellsInCutoffRadius, partPairsHandler);

	readPhaseSpaceHeader(&domain, 0);
	readPhaseSpace(&container, &lmu, &domain, &domainDecomposition);

	string destination = directory + "/molecules.inp";
	_logger->info() << "Writing output to: " << destination << endl;
	CheckpointWriter writer(1, destination, 0, false);
	writer.doOutput(&container, &domainDecomposition, &domain, 1, &lmu);
}

std::vector<double> MDGenerator::getRandomVelocity(double temperature) const {
	vector<double> v_;
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

bool MDGenerator::isInsideDomain(Domain* domain, double position[3]) {
	for (int i = 0; i < 3; i++) {
		if (position[i] < 0 || position[i] > domain->getGlobalLength(i)) {
			return false;
		}
	}

	return true;
}
