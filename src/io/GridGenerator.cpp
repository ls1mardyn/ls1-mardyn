#include "io/GridGenerator.h"

#include "Domain.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/Random.h"
#include "utils/xmlfileUnits.h"
#include "utils/generator/Objects.h"
#include "utils/generator/ObjectFactory.h"

#include <string>
#include <map>

using Log::global_log;
using namespace std;

void GridGenerator::readXML(XMLfileUnits& xmlconfig) {
	if(xmlconfig.changecurrentnode("lattice")) {
		_lattice.readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("basis")) {
		_basis.readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("latticeOrigin")) {
		xmlconfig.getNodeValueReduced("x", _origin[0]);
		xmlconfig.getNodeValueReduced("y", _origin[1]);
		xmlconfig.getNodeValueReduced("z", _origin[2]);
		global_log->info() << "Origin: " << _origin[0] << ", " << _origin[1] << ", " << _origin[2] << endl;
		xmlconfig.changecurrentnode("..");
	}
	if(xmlconfig.changecurrentnode("object")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		ObjectFactory object_factory;
		global_log->debug() << "Obj name: " << object_type << endl;
		_object = object_factory.create(object_type);
		if(_object == nullptr) {
			global_log->debug() << "Unknown object type: " << object_type << endl;
		}
		global_log->error() << "Created object of type: " << _object->getName() << endl;
		_object->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}
	_generator.init(_lattice, _basis, _origin, _object);
}

long unsigned int GridGenerator::readPhaseSpace(ParticleContainer* particleContainer, list<ChemicalPotential>* lmu,
		Domain* domain, DomainDecompBase* domainDecomp) {
	unsigned long numMolecules = 0;
	Molecule molecule;

	Ensemble* ensemble = _simulation.getEnsemble();
	Random rng;
	
	while(_generator.getMolecule(&molecule) > 0) {
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
		molecule.setid(numMolecules);
		bool inserted = particleContainer->addParticle(molecule);
		if(inserted){
			numMolecules++;
		}
	}
	global_log->info() << "Number of inserted molecules: " << numMolecules << endl;
	particleContainer->updateMoleculeCaches();
	//! @todo Get rid of the domain class calls at this place here...
	domain->setGlobalTemperature(ensemble->T());
	domain->setglobalNumMolecules(numMolecules);
	domain->setglobalRho(numMolecules / ensemble->V() );
	//! @todo reduce numMolecules?!
	return numMolecules;
}
