
#pragma once

#include <list>
#include <memory>

#include "utils/Random.h"
#include "molecules/Molecule.h"

class DomainDecompBase;
class ParticleContainer;
class CellProcessor;
class Domain;
class ParticleIterator;

//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
class ChemicalPotential {
public:
	ChemicalPotential();

	void setMu(int cid, double chempot) { _mu = chempot; _componentid = cid; }
	unsigned getInterval() { return _interval; }
	void setInterval(unsigned delta) { _interval = delta; }
	void setInstances(unsigned n) { _instances = n; }
	void setSystem(double x, double y, double z, double m);
	void setGlobalN(unsigned long N) { _globalN = N; }
	void setNextID(unsigned long id) { _nextid = id; }
	void setSubdomain(int rank, double x0, double x1, double y0, double y1, double z0, double z1);
	void setIncrement(unsigned idi) { _id_increment = idi; }

	void prepareTimestep(ParticleContainer* moleculeContainer, DomainDecompBase* comm);  // C must not contain the halo!

	// false if no deletion remains for this subdomain
	ParticleIterator getDeletion(ParticleContainer* moleculeContainer, double* minco, double* maxcoy);
	unsigned long getInsertion(double* ins);  // 0 if no insertion remains for this subdomain
	bool decideDeletion(double deltaUTilde);
	bool decideInsertion(double deltaUTilde);

	Molecule loadMolecule();
	void storeMolecule( Molecule& old )
	{
		if(hasSample()) return;
		mardyn_assert(old.componentid() == _componentid);
#ifndef NDEBUG
		old.check(old.getID());
#endif
		_reservoir = std::make_unique<Molecule>(old);
	}
	bool hasSample() { return _reservoir != nullptr; }

	void setPlanckConstant(double h_in) { _h = h_in; }
	void submitTemperature(double T_in);
	void setControlVolume(
			double x0, double y0, double z0, double x1, double y1, double z1
	);

	unsigned long getGlobalN() { return _globalN; }
	double getGlobalRho() { return (double)(_globalN) / _globalV; }

	void outputIX() { std::cout << "  r" << _ownrank << "[IX" << _rnd.getIX() << "]  "; }

	void assertSynchronization(DomainDecompBase* comm);

	double getMu() { return _mu; }
	unsigned int getComponentID() { return _componentid; }
	int rank() { return _ownrank; }

	void disableWidom() { _widom = false; }
	void enableWidom() { _widom = true; }
	bool isWidom() { return _widom; }

	double getLambda() { return _lambda; }
	float getDensityCoefficient() { return _decisive_density; }

	/* Moved from LinkedCells! */
	int getLocalGrandcanonicalBalance() {
		return _localInsertionsMinusDeletions;
	}
	/* Moved from LinkedCells! */
	void grandcanonicalStep(ParticleContainer * moleculeContainer, double T, Domain* domain, CellProcessor* cellProcessor);
	/* Moved from LinkedCells! */
	int grandcanonicalBalance(DomainDecompBase* comm);


private:
	//! @brief counts all particles inside the bounding box of this container
	unsigned countParticles(ParticleContainer * moleculeContainer, unsigned int cid) const;
	//! @brief counts particles in the intersection of bounding box and control volume
	unsigned countParticles(ParticleContainer * moleculeContainer, unsigned int cid, double * cbottom, double * ctop) const;

	bool moleculeStrictlyNotInBox(const Molecule& m, const double l[3], const double u[3]) const;

	int _ownrank;  // only for debugging purposes (indicate rank in console output)

	double _h;  // Plancksches Wirkungsquantum
	double _T;

	double _mu;
	double _muTilde;
	unsigned int _componentid;
	unsigned _interval;  // how often?
	unsigned _instances;  // how many trial insertions and deletions?
	Random _rnd, _rndmomenta;
	double _system[3];  // extent of the system
	float _minredco[3];  // minimal coordinates of the subdomain reduced w. r. t. the system size
	float _maxredco[3];   // maximal coordinates of the subdomain reduced w. r. t. the system size

	unsigned long _nextid;  // ID given to the next inserted particle
	unsigned _id_increment;
	std::list<unsigned> _remainingDeletions;  // position of the vector that should be deleted
	std::list<double> _remainingInsertions[3];
	std::list<unsigned long> _remainingInsertionIDs;
	std::list<float> _remainingDecisions;  // first deletions, then insertions

	unsigned long _globalN;
	double _globalV;
	double _molecularMass;
	double _globalReducedVolume;

	bool _restrictedControlVolume;
	double _control_bottom[3];
	double _control_top[3];

	float _decisive_density;

	double _lambda;

	bool _widom;  // Widom method -> determine mu (the chemical potential) by test insertions which are all rejected.
	// Using the widom method is just a way to determine the potential. This has nothing to do with actual insertions or
	// deletions!

	std::unique_ptr<Molecule> _reservoir;

	/* Moved from LinkedCells! */
	int _localInsertionsMinusDeletions; //!< balance of the grand canonical ensemble
};
