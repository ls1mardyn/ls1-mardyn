#ifndef TCTS_H_
#define TCTS_H_

#include "io/InputBase.h"


class MkTcTSGenerator : public InputBase {

public:
	MkTcTSGenerator(){}
	~MkTcTSGenerator(){}

	void setPhaseSpaceFile(std::string filename){}
	void setPhaseSpaceHeaderFile(std::string filename){}

	void readPhaseSpaceHeader(Domain* domain, double timestep){}
	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);

	void readXML(XMLfileUnits& xmlconfig);

private:
	double heigth1;

	double rho1;
	double rho2;
    void _msimulation(int arg1);
};

#endif /* TCTS_H_ */
