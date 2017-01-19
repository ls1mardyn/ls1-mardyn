
#ifndef GRANDCANONICAL_H_
#define GRANDCANONICAL_H_

#include <list>

#include "utils/Random.h"
#include "molecules/Molecule.h"

class DomainDecompBase;
class ParticleContainer;

typedef ParticleContainer TMoleculeContainer;

//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
class ChemicalPotential {
public:
	ChemicalPotential();

	void setMu(int cid, double chempot) { this->mu = chempot; this->componentid = cid; }
	unsigned getInterval() { return this->interval; }
	void setInterval(unsigned delta) { this->interval = delta; }
	void setInstances(unsigned n) { this->instances = n; }
	void setSystem(double x, double y, double z, double m);
	void setGlobalN(unsigned long N) { this->globalN = N; }
	void setNextID(unsigned long id) { this->nextid = id; }
	void setSubdomain(int rank, double x0, double x1, double y0, double y1, double z0, double z1);
	void setIncrement(unsigned idi) { this->id_increment = idi; }

	void prepareTimestep(TMoleculeContainer* cell, DomainDecompBase* comm);  // C must not contain the halo!

	// false if no deletion remains for this subdomain
	bool getDeletion(TMoleculeContainer* cell, double* minco, double* maxco);
	unsigned long getInsertion(double* ins);  // 0 if no insertion remains for this subdomain
	bool decideDeletion(double deltaUTilde);
	bool decideInsertion(double deltaUTilde);

	Molecule loadMolecule();
	void storeMolecule( Molecule& old )
	{
		if(hasSample()) return;
		assert(old.componentid() == componentid);
#ifndef NDEBUG
		old.check(old.id());
#endif
		this->reservoir = new Molecule(old);
	}
	bool hasSample() { return this->reservoir != NULL; }

	void setPlanckConstant(double h_in) { this->h = h_in; }
	void submitTemperature(double T_in);
	void setControlVolume(
			double x0, double y0, double z0, double x1, double y1, double z1
	);

	unsigned long getGlobalN() { return this->globalN; }
	double getGlobalRho() { return (double)(this->globalN) / this->globalV; }

	void outputIX() { std::cout << "  r" << ownrank << "[IX" << rnd.getIX() << "]  "; }

	void assertSynchronization(DomainDecompBase* comm);

	double getMu() { return this->mu; }
	unsigned int getComponentID() { return this->componentid; }
	int rank() { return this->ownrank; }

        void disableWidom() { this->widom = false; }
        void enableWidom() { this->widom = true; }
        bool isWidom() { return this->widom; }

	double getLambda() { return this->lambda; }
	float getDensityCoefficient() { return this->decisive_density; }

private:
	int ownrank;  // only for debugging purposes (indicate rank in console output)

	double h;  // Plancksches Wirkungsquantum
	double T;

	double mu;
	double muTilde;
	unsigned int componentid;
	unsigned interval;  // how often?
	unsigned instances;  // how many trial insertions and deletions?
	Random rnd, rndmomenta;
	double system[3];  // extent of the system
	float minredco[3];  // minimal coordinates of the subdomain reduced w. r. t. the system size
	float maxredco[3];   // maximal coordinates of the subdomain reduced w. r. t. the system size

	unsigned long nextid;  // ID given to the next inserted particle
	unsigned id_increment;
	std::list<unsigned> remainingDeletions;  // position of the vector that should be deleted
	std::list<double> remainingInsertions[3];
	std::list<unsigned long> remainingInsertionIDs;
	std::list<float> remainingDecisions;  // first deletions, then insertions

	unsigned long globalN;
	double globalV;
	double molecularMass;
	double globalReducedVolume;

	bool restrictedControlVolume;
	double control_bottom[3];
	double control_top[3];

	float decisive_density;

	double lambda;

        bool widom;  // Widom method -> determine mu by test insertions which are all rejected

	Molecule* reservoir;
};

#endif /* GRANDCANONICAL_H_ */
