/*
 * RDFForceIntegratorExtendedSite.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATOREXTENDEDSITE_H_
#define RDFFORCEINTEGRATOREXTENDEDSITE_H_

#include "RDFForceIntegrator.h"


class RDFForceIntegratorExtendedSite: public RDFForceIntegrator {
public:
	RDFForceIntegratorExtendedSite(ParticleContainer* moleculeContainer, double rc, std::vector<std::vector<double> >* globalADist,
			std::vector<std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegratorExtendedSite();

	void traverseMolecules();

private:
	static double _normal_lim[2], _dn, _dr, _extension, *_scaling_factors, _d_alpha, _d_level;
	static int _n_r, _n_n, _n_levels, _n_alpha;
	static bool called;
	void integrateRDFSite(Molecule* currentMolecule, double* normal_dim, int* boundary, int plane, unsigned int site);
	static void precomputeScalingFactors();
};

#endif /* RDFFORCEINTEGRATOREXTENDEDSITE_H_ */
