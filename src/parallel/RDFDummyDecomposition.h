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
	void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer,
			const std::vector<Component>& components, Domain* domain);
	void exchangeMolecules(ParticleContainer* moleculeContainer,
			const std::vector<Component>& components, Domain* domain);
private:
	int rdfBoundary;
	int particle_insertion_type;
	static bool have_avg_energy;
	static int num_calles;
	static double* unif_rand;
	static bool first_unif;
	ParticlePairsHandler* _particlePairsHandler;
	Simulation* simulation;

	double randdouble(double a, double b) const {
		return a + rand() * (b - a) / (RAND_MAX);
	}

	double getAverageEnergy(LinkedCells* linkedCells, double* rmin,
			double* rmax);
	void generateRandomVelocity(double temperature, double m, double* v);

	void generateRandomAngularVelocity(double temperature, double* w,
			Domain* domain, Molecule* currentMolecule);

	double getGaussianRandomNumber();
	double getUniformRandomNumber();

	void addPeriodicCopies(ParticleContainer* moleculeContainer, double* rmax,
			double* rmin, double* phaseSpaceSize, double* halo_L, const std::vector<Component>& components);
};

#endif /* RDFDUMMYDECOMPOSITION_H_ */
