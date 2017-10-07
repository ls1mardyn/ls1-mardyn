#ifndef SRC_IO_OBJECTGENERATOR_H_
#define SRC_IO_OBJECTGENERATOR_H_

#include <list>

#include "io/InputBase.h"

class GridFiller;
class VelocityAssignerBase;

class MultiObjectGenerator : public InputBase {
public:
	MultiObjectGenerator() : _defaultVelocityAssigner(nullptr) {};
	virtual ~MultiObjectGenerator() {
		delete _defaultVelocityAssigner;
	}

	/** @brief Read in XML configuration for Generator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * @note This structure is not fixed yet and may see changes
	 * \code{.xml}
	   <generator name="MultiObjectGenerator">
	     <objectgenerator><!-- ... --></objectgenerator>
	   </generator >
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {}

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);

private:

	std::list<GridFiller*> _generators;
	VelocityAssignerBase *_defaultVelocityAssigner;
};

#endif  // SRC_IO_OBJECTGENERATOR_H_
