#ifndef SRC_IO_OBJECTGENERATOR_H_
#define SRC_IO_OBJECTGENERATOR_H_

#include <list>

#include "io/InputBase.h"

class GridFiller;
class VelocityAssignerBase;

class ObjectGenerator : public InputBase {
public:
    ObjectGenerator() {};
    virtual ~ObjectGenerator() {}
	virtual void readXML(XMLfileUnits& xmlconfig);

	void setPhaseSpaceFile(std::string /*filename*/) {}

	void setPhaseSpaceHeaderFile(std::string /*filename*/) {}

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {}

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);

private:

	std::list<GridFiller*> _generators;
	VelocityAssignerBase *_velocityAssigner;
};

#endif  // SRC_IO_OBJECTGENERATOR_H_
