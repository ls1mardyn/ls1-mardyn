/*
 * RDFForceIntegratorSiteSimpleScale.h
 *
 *  Created on: Sep 6, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATORSITESIMPLESCALE_H_
#define RDFFORCEINTEGRATORSITESIMPLESCALE_H_

#include "RDFForceIntegrator.h"

class RDFForceIntegratorSiteSimpleScale: public RDFForceIntegrator {
public:
	RDFForceIntegratorSiteSimpleScale(ParticleContainer* moleculeContainer,
			double rc, std::vector<std::vector<double> >* globalADist,
			std::vector<std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegratorSiteSimpleScale();
	double traverseMolecules();
	double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true){}
private:
	static double _dn, _dr, _extension, _normal_lim[2];
	void integrateRDFSite(Molecule* currentMolecule, double* normal_dim,
			int* boundary, int plane, unsigned int site);
};

#endif /* RDFFORCEINTEGRATORSITESIMPLESCALE_H_ */
