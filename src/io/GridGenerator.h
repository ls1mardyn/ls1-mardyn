#ifndef GRIDGENERATOR_H_
#define GRIDGENERATOR_H_

#include "../utils/generator/GeneratorExt.h"
#include "io/InputBase.h"
#include "utils/generator/Basis.h"
#include "utils/generator/Lattice.h"


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
	Basis _basis;
	Lattice _lattice;
	double _origin[3];
	GeneratorExt _generator;
};

#endif /* GRIDGENERATOR_H_ */
