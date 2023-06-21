/*
 * PosNegComp.cpp
 *
 *  Created on: 03.12.2018
 *      Author: mheinen
 */

#include "PosNegComp.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"

// constructor and destructor
PosNegComp::PosNegComp()
{
}

PosNegComp::~PosNegComp()
{
}

void PosNegComp::readXML(XMLfileUnits& xmlconfig)
{
	_cid_ub.pos = 1;
	_cid_ub.neg = 1;
	_limitY.left = 0.;
	_limitY.right = 100.;
	xmlconfig.getNodeValue("cid_ub/pos", _cid_ub.pos);
	xmlconfig.getNodeValue("cid_ub/neg", _cid_ub.neg);
	xmlconfig.getNodeValue("cid_ub/ignore", _cid_ub.ignore);
	xmlconfig.getNodeValue("limit_y/left", _limitY.left);
	xmlconfig.getNodeValue("limit_y/right", _limitY.right);
}

void PosNegComp::init(ParticleContainer *particleContainer,
		  DomainDecompBase *domainDecomp, Domain *domain)
{
}

void PosNegComp::beforeForces(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep
)
{
	this->changeComponents(particleContainer);
}

void PosNegComp::afterForces(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep
)
{
	this->changeComponents(particleContainer);
}

void PosNegComp::changeComponents(ParticleContainer* particleContainer)
{
	Component* comp_pos = global_simulation->getEnsemble()->getComponent(_cid_ub.pos-1);
	Component* comp_neg = global_simulation->getEnsemble()->getComponent(_cid_ub.neg-1);

	double regionLowCorner[3], regionHighCorner[3];
	for (unsigned d = 0; d < 3; d++) {
		regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
		regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
	}
	// ensure that we do not iterate over things outside of the container.
	regionLowCorner[1] = std::max(_limitY.left, regionLowCorner[1]);
	regionHighCorner[1] = std::min(_limitY.right, regionHighCorner[1]);

	for (auto it = particleContainer->regionIterator(regionLowCorner, regionHighCorner, ParticleIterator::ALL_CELLS);
		 it.isValid(); ++it) {
		// check if particle should be ignored because of its cid
		const uint32_t cid_ub = it->componentid()+1;
		if(cid_ub == _cid_ub.ignore)
			continue;

		// check if particle within limitY
		double ry = it->r(1);
		if(ry <= _limitY.left || ry >= _limitY.right)
			continue;

		double vy = it->v(1);
		if(vy > 0.)
			it->setComponent(comp_pos);
		else
			it->setComponent(comp_neg);
	}
}
