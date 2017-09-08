#ifndef SRC_IO_GRIDGENERATOR_H_
#define SRC_IO_GRIDGENERATOR_H_

#include <list>

#include "io/InputBase.h"
#include "utils/generator/Generator.h"

class VelocityAssignerBase;

class GridGenerator : public InputBase {
public:
    GridGenerator() {};
    virtual ~GridGenerator() {}
	virtual void readXML(XMLfileUnits& xmlconfig);

	void setPhaseSpaceFile(std::string /*filename*/) {}

	void setPhaseSpaceHeaderFile(std::string /*filename*/) {}

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {}

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);

private:

	std::list<Generator*> _generators;
	VelocityAssignerBase *_velocityAssigner;
};

#endif  // SRC_IO_GRIDGENERATOR_H_
