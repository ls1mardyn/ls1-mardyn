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
		if(!_filler) {
			global_log->error() << "Object filler could not be created" << endl;
			Simulation::exit(1);
		}
		global_log->debug() << "Using object filler of type: " << _filler->getPluginName() << endl;
		_filler->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->error() << "No filler specified." << endl;
		Simulation::exit(1);
	}

	if(xmlconfig.changecurrentnode("object")) {
		std::string objectType;
		xmlconfig.getNodeValue("@type", objectType);
		global_log->debug() << "Obj name: " << objectType << endl;
		ObjectFactory objectFactory;
		_object = std::shared_ptr<Object>(objectFactory.create(objectType));
		if(!_object) {
			global_log->error() << "Unknown object type: " << objectType << endl;
		}
		global_log->debug() << "Created object of type: " << _object->getPluginName() << endl;
		_object->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->error() << "No object specified." << endl;
		Simulation::exit(1);
	}

	if(xmlconfig.changecurrentnode("velocityAssigner")) {
		std::string velocityAssignerName;
		xmlconfig.getNodeValue("@type", velocityAssignerName);
		global_log->info() << "Velocity assigner: " << velocityAssignerName << endl;
		if(velocityAssignerName == "EqualVelocityDistribution") {
			_velocityAssigner = std::make_shared<EqualVelocityAssigner>();
		} else if(velocityAssignerName == "MaxwellVelocityDistribution") {
			_velocityAssigner = std::make_shared<MaxwellVelocityAssigner>();
		} else {
			global_log->error() << "Unknown velocity assigner specified." << endl;
			Simulation::exit(1);
		}
		Ensemble* ensemble = _simulation.getEnsemble();
		global_log->info() << "Setting temperature for velocity assigner to " << ensemble->T() << endl;
		_velocityAssigner->setTemperature(ensemble->T());
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->warning() << "No velocityAssigner specified.  Will not change velocities provided by filler." << endl;
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
		Quaternion q(1.0, 0., 0., 0.); /* orientation of molecules has to be set to a value other than 0,0,0,0! */
		molecule.setq(q);
		if(_velocityAssigner) {
			_velocityAssigner->assignVelocity(&molecule);
		}
		bool inserted = particleContainer->addParticle(molecule);
		if(inserted){
			numMolecules++;
		}
	}
	return numMolecules;
}
