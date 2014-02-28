#ifndef MKESFERA_H
#define MKESFERA_H

#include "io/InputBase.h"

/** @brief Single droplet scenario generator.
 *
 * Creates a droplet with a given center position and radius inside the simulation box.
 * The density of the droplet and its surrounding can be chosen seperately.
 * The droplet and its surrounding consist both out of molecules from the component
 * with ID=1 in the xml input file.
 */
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

