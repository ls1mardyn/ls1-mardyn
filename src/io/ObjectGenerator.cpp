#include "io/ObjectGenerator.h"


#include "Domain.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "molecules/MoleculeIdPool.h"
#include "parallel/DomainDecompBase.h"
#include "utils/generator/GridFiller.h"
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
		_filler = std::make_shared<GridFiller>();
		_filler->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	}

	if(xmlconfig.changecurrentnode("object")) {
		std::string object_type;
		xmlconfig.getNodeValue("@type", object_type);
		ObjectFactory object_factory;
		global_log->debug() << "Obj name: " << object_type << endl;
		_object = std::shared_ptr<Object>(object_factory.create(object_type));
		global_log->info() << "COUNT : " << _object.use_count() << endl;
		if(_object == nullptr) {
			global_log->error() << "Unknown object type: " << object_type << endl;
		}
		global_log->debug() << "Created object of type: " << _object->getName() << endl;
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

	_filler->setObject(_object.get());
	_filler->setBoudingBox(bBoxMin, bBoxMax);
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
