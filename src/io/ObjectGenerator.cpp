#include "io/ObjectGenerator.h"

#include <limits>
#include <chrono>

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

		const int seed = [&]() -> int {
			bool enableRandomSeed = false;
			xmlconfig.getNodeValue("@enableRandomSeed", enableRandomSeed);
			if(enableRandomSeed) {
				/** A random seed for the velocity generator is created.
				 *  The current rank is added to make sure that, if multiple simulations are instantiated across
				 *  multiple ranks (for ex. when coupling to MaMiCo), the seeds will be unique if they are created
				 *  at the same time
				 */
				return std::chrono::duration_cast<std::chrono::seconds>
				(std::chrono::system_clock::now().time_since_epoch()).count() + _simulation.domainDecomposition().getRank();
			
			} else {
				return 0;
			}
		}();
		global_log->info() << "Seed for velocity assigner: " << seed << endl;
		if(velocityAssignerName == "EqualVelocityDistribution") {
			_velocityAssigner = std::make_shared<EqualVelocityAssigner>(0, seed);
		} else if(velocityAssignerName == "MaxwellVelocityDistribution") {
			_velocityAssigner = std::make_shared<MaxwellVelocityAssigner>(0, seed);
		} else {
			global_log->error() << "Unknown velocity assigner specified." << endl;
			Simulation::exit(1);
		}
		Ensemble* ensemble = _simulation.getEnsemble();
		global_log->info() << "Setting temperature for velocity assigner to " << ensemble->T() << endl;
		_velocityAssigner->setTemperature(ensemble->T());
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->warning() << "No velocityAssigner specified.  Will not change velocities provided by filler."
							  << endl;
	}
}


/** Class implementing a cuboidal bounding box.
 *
 * Implements a cuboid with open upper end, i.e. the upper border is excluded
 * from the volume in contrast to the normal cuboid, where it is included.
 */
class BoundingBox : public Cuboid {
public:
	BoundingBox(double bBoxMin[3], double bBoxMax[3]) : Cuboid(bBoxMin, bBoxMax) {}

	bool isInside(double r[3]) {
		return (lowerCorner(0) <= r[0] && r[0] < upperCorner(0))
			   && (lowerCorner(1) <= r[1] && r[1] < upperCorner(1))
			   && (lowerCorner(2) <= r[2] && r[2] < upperCorner(2));
	}

	std::string getName() { return std::string("BoundingBox"); };
};

unsigned long
ObjectGenerator::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {
	unsigned long numMolecules = 0;
	if(_moleculeIdPool == nullptr) {
		_moleculeIdPool = std::make_shared<MoleculeIdPool>(std::numeric_limits<unsigned long>::max(),
														   domainDecomp->getNumProcs(), domainDecomp->getRank());
	}

	double bBoxMin[3];
	double bBoxMax[3];
	domainDecomp->getBoundingBoxMinMax(domain, bBoxMin, bBoxMax);
	std::shared_ptr<Object> bBox = std::make_shared<BoundingBox>(bBoxMin, bBoxMax);
	std::shared_ptr<Object> boundedObject = std::make_shared<ObjectIntersection>(bBox, _object);
	_filler->setObject(boundedObject);
	_filler->init();

	Molecule molecule;
	unsigned long moleculeID = _moleculeIdPool->getNewMoleculeId();
	while(_filler->getMolecule(&molecule) > 0) {
		molecule.setid(moleculeID);
		if(_velocityAssigner) {
			_velocityAssigner->assignVelocity(&molecule);
		}
		// only add particle if it is inside of the own domain!
		if(particleContainer->isInBoundingBox(molecule.r_arr().data())) {
			particleContainer->addParticle(molecule, true, false);
			numMolecules++;
			moleculeID = _moleculeIdPool->getNewMoleculeId();
		}
	}
	return numMolecules;
}
