#ifndef SRC_IO_OBJECTGENERATOR_H_
#define SRC_IO_OBJECTGENERATOR_H_

#include <memory>

#include "io/InputBase.h"

class ObjectFillerBase;

class Object;

class VelocityAssignerBase;

class MoleculeIdPool;

/** @brief The ObjectGenerator sets up a phase space by filling volumetric objects.
 *
 * The idea of the ObjectGenerator is to create a composite 3D volumetric Object and fill this with molecules.
 * The molecule placement into the object is performed by a Filler. The assignment of molecule velocities is
 * performed by a VelocityAssigner. The molecule IDs are provided by a MoleculeIdPool.
 */
class ObjectGenerator : public InputBase {
public:
	ObjectGenerator() : _filler(nullptr), _object(nullptr), _velocityAssigner(nullptr), _moleculeIdPool(nullptr) {};

	/** @brief Read in XML configuration for ObjectGenerator and all its included objects.
	 * 
	 * The velocityAssigner can take two additional parameters. enableRandomSeed adds the option of having random
	 * intial molecule velocities at the beginning of the simulation. If a seed is specified, that value is used
	 * instead. Leaving both blank gives the default behaviour (seed = 0). Both cannot be nonzero simultaneously.
	 * The removeDrift parameter uses the function IOHelpers::removeMomentum() to remove overall drift from the 
	 * phasespace after initialisation, similar to what is done in CubicGridGenerator. This is off by default, as
	 * some experiments may want to load checkpoints with initial drift.
	 *
	 * The following XML object structure is handled by this method:
	 * @note This structure is not fixed yet and may see changes
	 * \code{.xml}
	   <objectgenerator>
	     <filler type="STRING"> <!-- see Filler documentation --> </filler>
	     <object type="STRING"> <!-- see Object documentation --> </object>
	     <velocityAssigner type="STRING" enableRandomSeed="BOOL" seed="LONG"> 
			<!-- see VelocityAssignerBase documentation --> </velocityAssigner>
		 <removeDrift>BOOL</removeDrift>
	   </objectgenerator>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

	void setFiller(std::shared_ptr<ObjectFillerBase> filler) { _filler = filler; }

	void setObject(std::shared_ptr<Object> object) { _object = object; }

	void setVelocityAssigner(std::shared_ptr<VelocityAssignerBase> vAssigner) { _velocityAssigner = vAssigner; }

	void setMoleculeIDPool(std::shared_ptr<MoleculeIdPool> moleculeIdPool) { _moleculeIdPool = moleculeIdPool; }

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {}

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);

private:
	std::shared_ptr<ObjectFillerBase> _filler;
	std::shared_ptr<Object> _object;
	std::shared_ptr<VelocityAssignerBase> _velocityAssigner;
	std::shared_ptr<MoleculeIdPool> _moleculeIdPool;
	bool _removeDrift = false;
};

#endif  // SRC_IO_OBJECTGENERATOR_H_
