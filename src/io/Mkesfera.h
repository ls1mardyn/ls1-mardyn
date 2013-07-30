#ifndef MKESFERA_H
#define MKESFERA_H

#include "utils/OptionParser.h"
#include <list>

class Domain;
class DomainDecompBase;
class Integrator;
class OutputBase;
class ParticleContainer;
class Simulation;

class Mkesfera {

public:
	Mkesfera(optparse::Values &options);

	void generate(Domain* domain, DomainDecompBase** domainDecomposition, Integrator** integrator, ParticleContainer** moleculeContainer, std::list<OutputBase*> &outputPlugins, Simulation* simulation);

private:
	double cutoff;
	bool do_shift;
	double R_i, R_o;
	double rho_i, rho_o;
	double T;
	double box[3];

};

#endif /* MKESFERA_H */

