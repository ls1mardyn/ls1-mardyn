/*
 * RDFDummyDecomposition.h
 *
 * To be used if RDF boundaries are used, makes sure no halo layers are created
 *  Created on: Jul 5, 2012
 *      Author: tijana
 */
#include "parallel/DomainDecompDummy.h"
#include "ensemble/PressureGradient.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/LinkedCells.h"

#ifndef RDFDUMMYDECOMPOSITION_H_
#define RDFDUMMYDECOMPOSITION_H_
class Simulation;
class RDFDummyDecomposition: public DomainDecompDummy {
public:
	RDFDummyDecomposition(ParticlePairsHandler* ph, int boundary, int insertion_type, Simulation* sim);
	virtual ~RDFDummyDecomposition();

	/*
	 * calls exchange moleucles, to fulfill the interface
	 */
	void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer,
			const std::vector<Component>& components, Domain* domain);

	/*
	 * interface function, calls addPeriodicCopies
	 */
	void exchangeMolecules(ParticleContainer* moleculeContainer,
			const std::vector<Component>& components, Domain* domain);

	/*
	 * Insert with usher, currently in the strip of width rc from the boundary
	 */
	void insertUsher(ParticleContainer* moleculeContainer,
			const std::vector<Component>& components, Domain* domain);
private:
	int rdfBoundary; // boundary
	int particle_insertion_type; // 0 - wall, 1 - usher
	static bool have_avg_energy;
	static int num_calles; // making sure usher is not called in the dummy timestep when molecules are read
	static double* unif_rand; // random numbers
	static bool first_unif; // if there already is a random number
	ParticlePairsHandler* _particlePairsHandler; // needed for linked cells
	Simulation* simulation;

	double randdouble(double a, double b) const {
		return a + rand() * (b - a) / (RAND_MAX);
	}

	// average energy per molecule for molecules in the rmin to rmax region
	double getAverageEnergy(LinkedCells* linkedCells, double* rmin,
			double* rmax);

	// random translation velocity sampled according to the temperature
	void generateRandomVelocity(double temperature, double m, double* v);

	// random angular velocity sampled according to the temperature
	void generateRandomAngularVelocity(double temperature, double* w,
			Domain* domain, Molecule* currentMolecule);

	double getGaussianRandomNumber();
	double getUniformRandomNumber();

	/*
	 * Adds periodic copies in all direcontions except rdf boundary direction
	 * @param moleculeContainer particle container
	 * @param rmax box size maximal dimensions
	 * @param rmin box size minimal dimensions
	 * @param phaseSpaceSize rmax-rmin in each dimension
	 * @param halo_L halo size
	 * @param components, needed when creating periodic copies
	 */
	void addPeriodicCopies(ParticleContainer* moleculeContainer, double* rmax,
			double* rmin, double* phaseSpaceSize, double* halo_L, const std::vector<Component>& components);
};

#endif /* RDFDUMMYDECOMPOSITION_H_ */
