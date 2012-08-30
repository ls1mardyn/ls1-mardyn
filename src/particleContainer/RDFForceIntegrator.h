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
	RDFForceIntegrator(ParticleContainer* moleculeContainer, double rc, std::vector<std::vector<double> >* globalADist,
			std::vector<std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegrator();

	virtual void traverseMolecules() = 0;



protected:
	ParticleContainer* _moleculeContainer;
	double _rc;
	double _high_limit[3], _low_limit[3], _rmax[3], _rmin[3];
	int _numMolecules;
	std::vector<std::vector<double> >* _globalADist;
	std::vector<std::vector<std::vector<double> > >* _globalSiteADist;

	void initTraversal();
};

#endif /* RDFFORCEINTEGRATOR_H_ */
