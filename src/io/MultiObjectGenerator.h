#ifndef SRC_IO_MULTIOBJECTGENERATOR_H_
#define SRC_IO_MULTIOBJECTGENERATOR_H_

#include <list>
#include <memory>

#include "io/InputBase.h"

class ObjectGenerator;
class VelocityAssignerBase;

class MultiObjectGenerator : public InputBase {
public:
	MultiObjectGenerator() : _globalNumMolecules(0) {}
	~MultiObjectGenerator();

	/** @brief Read in XML configuration for Generator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * @note This structure is not fixed yet and may see changes
	 * \code{.xml}
	   <generator name="MultiObjectGenerator">
	     <objectgenerator> <!-- ... --> </objectgenerator>
	     ...
	     <velocityAssigner> <!-- ... --> </velocityAssigner>
	   </generator >
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {}

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);

private:

	std::list<ObjectGenerator*> _generators;
	unsigned long _globalNumMolecules;
};

#endif  // SRC_IO_MULTIOBJECTGENERATOR_H_
