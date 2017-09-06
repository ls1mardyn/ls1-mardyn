#include "io/GridGenerator.h"

#include "Domain.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"
#include "utils/Random.h"
#include "utils/xmlfileUnits.h"
#include "molecules/MoleculeIdPool.h"

#include <cmath>
#include <limits>
#include <map>
#include <string>

#include <mpi.h>

using Log::global_log;
using namespace std;


void GridGenerator::readXML(XMLfileUnits& xmlconfig) {
	XMLfile::Query query = xmlconfig.query("subgenerator");
	global_log->info() << "Number of sub-generators: " << query.card() << endl;
	string oldpath = xmlconfig.getcurrentnodepath();
	for( auto generatorIter = query.begin(); generatorIter; ++generatorIter ) {
		xmlconfig.changecurrentnode(generatorIter);
		_generators.push_back(new Generator);
		_generators.back()->readXML(xmlconfig);
	}
	xmlconfig.changecurrentnode(oldpath);
}

long unsigned int GridGenerator::readPhaseSpace(ParticleContainer* particleContainer, list<ChemicalPotential>* lmu,
		Domain* domain, DomainDecompBase* domainDecomp) {
	unsigned long numMolecules = 0;

	Ensemble* ensemble = _simulation.getEnsemble();
	Random rng;
	double bBoxMin[3];
	double bBoxMax[3];
	domainDecomp->getBoundingBoxMinMax(domain, bBoxMin, bBoxMax);
	MoleculeIdPool moleculeIdPool(std::numeric_limits<unsigned long>::max(), domainDecomp->getNumProcs(), domainDecomp->getRank());

	for(auto generator : _generators) {
		Molecule molecule;
		generator->setBoudingBox(bBoxMin, bBoxMax);
		generator->init();
		while(generator->getMolecule(&molecule) > 0) {
			double v_abs = sqrt(/*kB=1*/ ensemble->T() / molecule.component()->m());
			double phi, theta;
			phi = rng.rnd();
			theta = rng.rnd();
			double v[3];
			v[0] = v_abs * sin(phi);
			v[1] = v_abs * cos(phi) * sin(theta);
			v[2] = v_abs * cos(phi) * cos(theta);
			for(int d = 0; d < 3; d++) {
				molecule.setv(d, v[d]);
			}
			Quaternion q(1.0, 0., 0., 0.); /* orientation of molecules has to be set to a value other than 0,0,0,0! */
			molecule.setq(q);
			molecule.setid(moleculeIdPool.getNewMoleculeId());
			bool inserted = particleContainer->addParticle(molecule);
			if(inserted){
				numMolecules++;
			}
		}
	}
	global_log->debug() << "Number of locally inserted molecules: " << numMolecules << endl;
	particleContainer->updateMoleculeCaches();
	unsigned long globalNumMolecules = 0;
	MPI_Allreduce(&numMolecules, &globalNumMolecules, 1, MPI_UNSIGNED_LONG, MPI_SUM, domainDecomp->getCommunicator());
	global_log->info() << "Number of inserted molecules: " << numMolecules << endl;
	//! @todo Get rid of the domain class calls at this place here...
	domain->setGlobalTemperature(ensemble->T());
	domain->setglobalNumMolecules(globalNumMolecules);
	domain->setglobalRho(numMolecules / ensemble->V() );
	//! @todo reduce numMolecules?!
	return numMolecules;
}
