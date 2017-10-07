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


using Log::global_log;
using namespace std;

MultiObjectGenerator::MultiObjectGenerator::~MultiObjectGenerator() {}


void MultiObjectGenerator::readXML(XMLfileUnits& xmlconfig) {
	if(xmlconfig.changecurrentnode("velocityAssigner")) {
		std::string defaultVelocityAssignerName;
		xmlconfig.getNodeValue("@type", defaultVelocityAssignerName);
		if(defaultVelocityAssignerName == "EqualVelocityDistribution") {
			_defaultVelocityAssigner = std::make_shared<EqualVelocityAssigner>();
		} else if(defaultVelocityAssignerName == "MaxwellVelocityDistribution") {
			_defaultVelocityAssigner = std::make_shared<MaxwellVelocityAssigner>();
		}
		xmlconfig.changecurrentnode("..");
	}

	XMLfile::Query query = xmlconfig.query("objectgenerator");
	global_log->info() << "Number of sub-objectgenerators: " << query.card() << endl;
	string oldpath = xmlconfig.getcurrentnodepath();
	for( auto generatorIter = query.begin(); generatorIter; ++generatorIter ) {
		xmlconfig.changecurrentnode(generatorIter);
		ObjectGenerator *generator = new ObjectGenerator();
		generator->setVelocityAssigner(_defaultVelocityAssigner);
		generator->readXML(xmlconfig);
		_generators.push_back(generator);
	}
	xmlconfig.changecurrentnode(oldpath);
}


long unsigned int MultiObjectGenerator::readPhaseSpace(ParticleContainer* particleContainer, list<ChemicalPotential>* lmu,
		Domain* domain, DomainDecompBase* domainDecomp) {
	unsigned long numMolecules = 0;
	std::shared_ptr<MoleculeIdPool> moleculeIdPool = std::make_shared<MoleculeIdPool>(std::numeric_limits<unsigned long>::max(), domainDecomp->getNumProcs(), domainDecomp->getRank());

	Ensemble* ensemble = _simulation.getEnsemble();
	_defaultVelocityAssigner->setTemperature(ensemble->T());
	for(auto generator : _generators) {
		generator->setMoleculeIDPool(moleculeIdPool);
		numMolecules += generator->readPhaseSpace(particleContainer, lmu, domain, domainDecomp);
	}
	particleContainer->updateMoleculeCaches();
	unsigned long globalNumMolecules = numMolecules;
#ifdef ENABLE_MPI
	MPI_Allreduce(MPI_IN_PLACE, &globalNumMolecules, 1, MPI_UNSIGNED_LONG, MPI_SUM, domainDecomp->getCommunicator());
#endif
	global_log->debug() << "Number of locally inserted molecules: " << numMolecules << endl;
	global_log->info() << "Number of inserted molecules: " << globalNumMolecules<< endl;
	//! @todo Get rid of the domain class calls at this place here...
	domain->setGlobalTemperature(ensemble->T());
	domain->setglobalNumMolecules(globalNumMolecules);
	domain->setglobalRho(numMolecules / ensemble->V() );
	return numMolecules;
}
