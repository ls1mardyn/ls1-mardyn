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

class RDFForceIntegratorExact: public RDFForceIntegrator {
public:
	RDFForceIntegratorExact(ParticleContainer* moleculeContainer, double rc, double d,
			std::vector<std::vector<double> >* globalADist, std::vector<
					std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegratorExact();

	double traverseMolecules();

	void getScalingFactor(double* mol_r, double* site_r,
			double x, double y, double z, int site_i, double* scale);

	double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true, bool unit_test = false);


	void prepareUnitTest(Molecule* m) {
		_d = _d_level = 10;
		_d_alpha = 180;
		_rho = 1;
		_extension = m->ljcenter_disp(0);

	}
	double* precomputeScalingFactorsX(bool unit_test = false);
private:
	double  _extension, *_scaling_factors_x,  _d_alpha, _d_level, _rho, _g_start;
	int _n_r, _n_n, _n_levels, _n_alpha;
	bool called_x;
	int timestep;
	double* rhos;
	bool first_unif;
	double unif_rand[2];
	double getGaussianRandomNumber();

	double integrateRDFSiteCartesian(double xlim[2], double ylim[2],
			double zlim[2], Molecule* mol, int plane, unsigned int site,
			int boundary[3], bool add_influence, double* return_force, bool unit_test = false);

	double checkScalingFactor(int idx_level,
			int idx_n, int idx_r);

};

#endif /* RDFFORCEINTEGRATOREXACT_H_ */
