/*
 * RDFForceIntegratorSite.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATORSITE_H_
#define RDFFORCEINTEGRATORSITE_H_

#include "RDFForceIntegrator.h"
#include "molecules/potforce.h"

class RDFForceIntegratorSite: public RDFForceIntegrator {
public:
	RDFForceIntegratorSite(ParticleContainer* moleculeContainer, double rc, double d, std::vector<std::vector<double> >* globalADist,
			std::vector<std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegratorSite();

	double traverseMolecules();

	double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true, bool unit_test = false);


private:
	double _rho;
	double integrateRDFSite(Molecule* currentMolecule, double* normal_dim, int* boundary, int plane, unsigned int site, double* force, bool add_influence);
	void integrateRDFSiteCartesian(double xlim[2], double ylim[2],
			double zlim[2], Molecule* mol, int plane, unsigned int site,
			int boundary[3]);

};

#endif /* RDFFORCEINTEGRATORSITE_H_ */
