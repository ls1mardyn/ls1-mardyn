#include "io/MultiObjectGenerator.h"
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <string>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "Domain.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "io/ObjectGenerator.h"
#include "molecules/MoleculeIdPool.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/generator/EqualVelocityAssigner.h"
#include "utils/generator/MaxwellVelocityAssigner.h"
#include "utils/generator/VelocityAssignerBase.h"
#include "utils/xmlfileUnits.h"


MultiObjectGenerator::MultiObjectGenerator::~MultiObjectGenerator() {
	for(auto& generator : _generators) {
		delete generator;
	}
}


void MultiObjectGenerator::readXML(XMLfileUnits& xmlconfig) {

	XMLfile::Query query = xmlconfig.query("objectgenerator");
	Log::global_log->info() << "Number of sub-objectgenerators: " << query.card() << std::endl;
	std::string oldpath = xmlconfig.getcurrentnodepath();
	for(auto generatorIter = query.begin(); generatorIter; ++generatorIter) {
		xmlconfig.changecurrentnode(generatorIter);
		ObjectGenerator* generator = new ObjectGenerator();
		generator->readXML(xmlconfig);
		_generators.push_back(generator);
	}
	xmlconfig.changecurrentnode(oldpath);
}


unsigned long MultiObjectGenerator::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain,
												   DomainDecompBase* domainDecomp) {
	unsigned long numMolecules = 0;
	std::shared_ptr<MoleculeIdPool> moleculeIdPool = std::make_shared<MoleculeIdPool>(
			std::numeric_limits<unsigned long>::max(), domainDecomp->getNumProcs(), domainDecomp->getRank());

	for(auto generator : _generators) {
		generator->setMoleculeIDPool(moleculeIdPool);
		numMolecules += generator->readPhaseSpace(particleContainer, domain, domainDecomp);
	}
	particleContainer->updateMoleculeCaches();
	Log::global_log->info() << "Number of locally inserted molecules: " << numMolecules << std::endl;
	_globalNumMolecules = numMolecules;
#ifdef ENABLE_MPI
	MPI_Allreduce(MPI_IN_PLACE, &_globalNumMolecules, 1, MPI_UNSIGNED_LONG, MPI_SUM, domainDecomp->getCommunicator());
#endif
	Log::global_log->info() << "Number of globally inserted molecules: " << _globalNumMolecules << std::endl;
	//! @todo Get rid of the domain class calls at this place here...
	return _globalNumMolecules;
}
