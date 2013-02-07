/*
 * RDFForceIntegratorExact.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATOREXACT_H_
#define RDFFORCEINTEGRATOREXACT_H_

#include "RDFForceIntegrator.h"
#include "molecules/potforce.h"
#include "RDF.h"

/**
 * Molecule-based integration scheme for the
 * rdf-based boundary force
 * To be used with the RDFDomainDecompDummy class
 */
class RDFForceIntegratorExact: public RDFForceIntegrator {
public:
	/*
	 * constructor
	 * @param moleculeContainer particle container
	 * @param rc cutoff radius
	 * @param d size of the quadrature
	 * @param globalADist a vector containing molecule-molecule rdf vectors per component (normally 1 component so only one vector)
	 * @param globalSiteADist a vector containing site-site rdfs per component (1 component normally), then 0-0, 0-1, 1-0, 1-1 interaction
	 */
	RDFForceIntegratorExact(ParticleContainer* moleculeContainer, double rc, double d,
			std::vector<std::vector<double> >* globalADist, std::vector<
					std::vector<std::vector<double> > >* globalSiteADist,
					double randomizeValue = 0, double randomizePercentage = 0, int randomizeNumSteps = 1);
	virtual ~RDFForceIntegratorExact();

	/*
	 * go through all the molecules and compute the force if close to boundary
	 * ATTENTION: right now it is considered that the open (rdf) boundary is in x direction
	 */
	double traverseMolecules();

	/*
	 * This function processes one molecule
	 * @param currentMolecule molecule to process
	 * @param force output for the force on that molecule
	 * @param add_influence if the influence should be added. please change this to false if you just want to see the quality of the scheme
	 * @param unit_test if the unit testing code is ran, linear rdf is used. do not use this
	 */
	double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true, bool unit_test = false);


	/*
	 * Just for the unit test, making sure there are only few integration points
	 */
	void prepareUnitTest(Molecule* m) {
		_d = _d_level = 10;
		_d_alpha = 180;
		_rho = 1;
		_extension = m->ljcenter_disp(0);

	}

	/*
	 * precomputing the scaling factors
	 * considers the boundary to be in x direction
	 * @param unit_test this adapts the scaling factors for the unit test, don't use this normally
	 */
	double* precomputeScalingFactorsX(bool unit_test = false);
private:
	double  _extension; // distance between center and site
	double*_scaling_factors_x; // scaling factors, x direction
	double _d_alpha; // angle spacing for the scaling factor computation
	double _d_level; // level spacing for the scaling factor computation
	double _rho; // density
	int _n_r, _n_n, _n_levels, _n_alpha; // number of points in normal and radial direction, of levels and angles for the scaling factor computation
	bool called_x; // if computing scaling factors was already called
	int timestep; // counts timesteps
	bool first_unif; // if there already is one random number
	double unif_rand[2]; // stores two uniform random numbers
	double _randomizeValue; // randomize between (-value, value)
	double _randomizePercentage; // randomize between (-pF, pF), where p is percentage, F force
	int _randomizeNumSteps; // randomize every _randomizeNumSteps step, default value every step so 1
	/*
	 * gets gaussian random number
	 * to be used for force randomizations
	 */
	double getGaussianRandomNumber();

	/*
	 * This method performs the integration
	 * ATTENTION: considers rdf boundary in x direction
	 * @param xlim integration limits in x direction
	 * @param ylim integration limits in y direction
	 * @param zlim integration limits in z direction
	 * @param mol molecule it is integrated for
	 * @param site which site is integrated for
	 * @param boundary for each direction, stores -1 or 1
	 * @param add_influence if false, the force is not added to the molecule
	 * @param return_force output containing the force on this molecule
	 * @param unit_test do not use this, for unit testing
	 */
	double integrateRDFSiteCartesian(double xlim[2], double ylim[2],
			double zlim[2], Molecule* mol, int plane, unsigned int site,
			int boundary[3], bool add_influence, double* return_force, bool unit_test = false);


};

#endif /* RDFFORCEINTEGRATOREXACT_H_ */
