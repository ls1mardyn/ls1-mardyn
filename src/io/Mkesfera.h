#ifndef MKESFERA_H
#define MKESFERA_H

#include "io/InputBase.h"


class MkesferaGenerator : public InputBase {

public:
	MkesferaGenerator(){}
	~MkesferaGenerator(){}

	void setPhaseSpaceFile(std::string filename){}
	void setPhaseSpaceHeaderFile(std::string filename){}

	void readPhaseSpaceHeader(Domain* domain, double timestep){}
	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);

	void readXML(XMLfileUnits& xmlconfig);

private:
	double R_i, R_o;
	double rho_i, rho_o;
	double center[3]; /**< droplet center */
};

#endif /* MKESFERA_H */

