/*
 * RDFForceIntegrator.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#include "ParticleContainer.h"
#include <vector>
#include <cstdlib>
#include "molecules/Molecule.h"

#define PI 3.1415926535

#ifndef RDFFORCEINTEGRATOR_H_
#define RDFFORCEINTEGRATOR_H_

class RDFForceIntegrator {
public:
	RDFForceIntegrator(ParticleContainer* moleculeContainer, double rc, double d, std::vector<std::vector<double> >* globalADist,
			std::vector<std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegrator();

	virtual double traverseMolecules() = 0;

	virtual double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true, bool unit_test = false) = 0;


protected:
	ParticleContainer* _moleculeContainer; // container of the molecules
	double _rc, _d; // cutoff radius and quadrature spacing
	double _high_limit[3], _low_limit[3]; // if the molecule is within these limits, it doesn't need to be integrated for
	double _rmax[3], _rmin[3]; // domain boundaries
	int _numMolecules; // number of molecules
	std::vector<std::vector<double> >* _globalADist; // moleucle-molecule rdf (accumulated)
	std::vector<std::vector<std::vector<double> > >* _globalSiteADist; // site-site rdf (accumulated)

	void initTraversal();
};

#endif /* RDFFORCEINTEGRATOR_H_ */
