/*
 * RDFForceIntegrator.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#include "RDFForceIntegrator.h"
#include <cmath>


RDFForceIntegrator::RDFForceIntegrator(ParticleContainer* moleculeContainer, double rc, double d, std::vector<std::vector<double> >* globalADist,
		std::vector<std::vector<std::vector<double> > >* globalSiteADist) {

	_moleculeContainer = moleculeContainer;
	_rc = rc;
	_d = d;
	_globalADist = globalADist;
	_globalSiteADist = globalSiteADist;
	initTraversal();
}

RDFForceIntegrator::~RDFForceIntegrator() {
	// TODO Auto-generated destructor stub
}

void RDFForceIntegrator::initTraversal() {
	for (int i = 0; i < 3; i++) {
		_rmin[i] = _moleculeContainer->getBoundingBoxMin(i);
		_low_limit[i] = _rmin[i] + _rc;
		_rmax[i] = _moleculeContainer->getBoundingBoxMax(i);
		_high_limit[i] = _rmax[i] - _rc;
	}
}
