#include "io/ObjectGenerator.h"

#include <limits>

#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "molecules/MoleculeIdPool.h"
#include "parallel/DomainDecompBase.h"
#include "utils/generator/ObjectFillerBase.h"
#include "utils/generator/ObjectFillerFactory.h"
#include "utils/generator/ObjectFactory.h"
#include "utils/generator/EqualVelocityAssigner.h"
#include "utils/generator/MaxwellVelocityAssigner.h"
#include "utils/generator/VelocityAssignerBase.h"


using std::endl;

void ObjectGenerator::readXML(XMLfileUnits& xmlconfig) {
	if(xmlconfig.changecurrentnode("filler")) {
		std::string fillerType;
		xmlconfig.getNodeValue("@type", fillerType);
		global_log->debug() << "Filler type: " << fillerType << endl;
		ObjectFillerFactory objectFillerFactory;
		_filler = std::shared_ptr<ObjectFillerBase>(objectFillerFactory.create(fillerType));
		if(_filler == nullptr) {
			global_log->error() << "Object Filler could not be created" << endl;
			Simulation::exit(1);
		}
		global_log->debug() << "Using object filler of type: " << _filler->getPluginName() << endl;
		_filler->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}

	if(xmlconfig.changecurrentnode("object")) {
		std::string objectType;
		xmlconfig.getNodeValue("@type", objectType);
		global_log->debug() << "Obj name: " << objectType << endl;
		ObjectFactory objectFactory;
		_object = std::shared_ptr<Object>(objectFactory.create(objectType));
		if(_object == nullptr) {
			global_log->error() << "Unknown object type: " << objectType << endl;
		}
		global_log->debug() << "Created object of type: " << _object->getPluginName() << endl;
		_object->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}

	if(xmlconfig.changecurrentnode("velocityAssigner")) {
		std::string defaultVelocityAssignerName;
		xmlconfig.getNodeValue("@type", defaultVelocityAssignerName);
		if(defaultVelocityAssignerName == "EqualVelocityDistribution") {
			_velocityAssigner = std::make_shared<EqualVelocityAssigner>();
		} else if(defaultVelocityAssignerName == "MaxwellVelocityDistribution") {
			_velocityAssigner = std::make_shared<MaxwellVelocityAssigner>();
		}
		Ensemble* ensemble = _simulation.getEnsemble();
		_velocityAssigner->setTemperature(ensemble->T());
		xmlconfig.changecurrentnode("..");
	}
}

unsigned long ObjectGenerator::readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {
	unsigned long numMolecules = 0;
	if(_moleculeIdPool == nullptr) {
		_moleculeIdPool = std::make_shared<MoleculeIdPool>(std::numeric_limits<unsigned long>::max(), domainDecomp->getNumProcs(), domainDecomp->getRank());
	}

	double bBoxMin[3];
	double bBoxMax[3];
	domainDecomp->getBoundingBoxMinMax(domain, bBoxMin, bBoxMax);
	std::shared_ptr<Object> bBox = std::make_shared<Cuboid>(bBoxMin, bBoxMax);
	std::shared_ptr<Object> boundedObject = std::make_shared<ObjectIntersection>(bBox, _object);
	_filler->setObject(boundedObject.get());
	_filler->init();

	Molecule molecule;
	while(_filler->getMolecule(&molecule) > 0) {
		molecule.setid(_moleculeIdPool->getNewMoleculeId());
		_velocityAssigner->assignVelocity(&molecule);
		Quaternion q(1.0, 0., 0., 0.); /* orientation of molecules has to be set to a value other than 0,0,0,0! */
		molecule.setq(q);
		bool inserted = particleContainer->addParticle(molecule);
		if(inserted){
			numMolecules++;
		}
	}
	return numMolecules;
}
