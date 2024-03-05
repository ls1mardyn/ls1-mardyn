/*
 * ExtractPhase.cpp
 *
 *  Created on: 31.12.2019
 *      Author: mheinen
 */

#include "ExtractPhase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "utils/CommVar.h"

ExtractPhase::ExtractPhase()
{
	_bDone.beforeForces = false;
	_bDone.afterForces = false;
}

ExtractPhase::~ExtractPhase()
{
}

void ExtractPhase::readXML(XMLfileUnits& xmlconfig)
{
	// Component change
	_compChange.enabled = false;
	xmlconfig.getNodeValue("change/enabled", _compChange.enabled);
	xmlconfig.getNodeValue("change/cid_ub", _compChange.cid_ub);

	// Density target
	_densityTarget.value = 0.7;
	_densityTarget.percent = 1.0;
	_densityTarget.cutoff = 2.5;
	_densityTarget.range.left = 0.;
	_densityTarget.range.right = 100.;
	xmlconfig.getNodeValue("density/value", _densityTarget.value);
	xmlconfig.getNodeValue("density/percent", _densityTarget.percent);
	xmlconfig.getNodeValue("density/cutoff", _densityTarget.cutoff);
	bool bRet = true;
	bRet = bRet && xmlconfig.getNodeValue("density/range/left", _densityTarget.range.left);
	bRet = bRet && xmlconfig.getNodeValue("density/range/right", _densityTarget.range.right);
	_densityTarget.range.enabled = bRet;

	// Interface params
	_interface.type = INTT_VAPOR_LIQUID;
	_interface.range.left = 0.;
	_interface.range.right = 100.;
	std::string str = "unknown";
	xmlconfig.getNodeValue("interface/type", str);
	xmlconfig.getNodeValue("interface/range/left", _interface.range.left);
	xmlconfig.getNodeValue("interface/range/right", _interface.range.right);
	if(str == "1" || str == "vl" || str == "VL")
		_interface.type = INTT_VAPOR_LIQUID;
	else if (str == "2" || str == "lv" || str == "LV")
		_interface.type = INTT_LIQUID_VAPOR;
}

void ExtractPhase::init(ParticleContainer *particleContainer,
		  DomainDecompBase *domainDecomp, Domain *domain)
{
}

void ExtractPhase::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	if(_bDone.beforeForces)
		return;

	double regionLowCorner[3], regionHighCorner[3];
	for (unsigned d = 0; d < 3; d++) {
		regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
		regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
	}
	// ensure that we do not iterate over things outside of the container.
	regionLowCorner[1] = std::max(_densityTarget.range.left, regionLowCorner[1]);
	regionHighCorner[1] = std::min(_densityTarget.range.right, regionHighCorner[1]);

	double A_xz, V;
	Domain* domain = global_simulation->getDomain();
	A_xz = domain->getGlobalLength(0) * domain->getGlobalLength(2);
	V = A_xz * (_densityTarget.range.right - _densityTarget.range.left);
	CommVar<uint64_t> numParticles;
	numParticles.local = 0;

	for (auto it = particleContainer->regionIterator(regionLowCorner, regionHighCorner,
													 ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		 it.isValid(); ++it)
		numParticles.local++;

	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(numParticles.local);
	domainDecomp->collCommAllreduceSum();
	numParticles.global = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
	_densityTarget.value = numParticles.global / V;

	// Perform action only once
	_bDone.beforeForces = true;
}

void ExtractPhase::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	if(_bDone.afterForces)
		return;

	double regionLowCorner[3], regionHighCorner[3];
	for (unsigned d = 0; d < 3; d++) {
		regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
		regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
	}

	// Delete particles of vapor phase outside of interface range
	{
		// ensure that we do not iterate over things outside of the container.
		if(INTT_LIQUID_VAPOR == _interface.type)
			regionLowCorner[1] = std::max(_interface.range.right, regionLowCorner[1]);
		if(INTT_VAPOR_LIQUID == _interface.type)
			regionHighCorner[1] = std::min(_interface.range.left, regionHighCorner[1]);

		for (auto it = particleContainer->regionIterator(regionLowCorner, regionHighCorner,
														 ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			 it.isValid(); ++it) {
			particleContainer->deleteMolecule(it, false);
		}
	}

	// Delete/mark particles of vapor particles inside the interface range
	// ensure that we do not iterate over things outside of the container.
	for (unsigned d = 0; d < 3; d++) {
		regionLowCorner[d] = particleContainer->getBoundingBoxMin(d) - _densityTarget.cutoff;
		regionHighCorner[d] = particleContainer->getBoundingBoxMax(d) + _densityTarget.cutoff;
	}
	regionLowCorner[1] = std::max(_interface.range.left - _densityTarget.cutoff, regionLowCorner[1]);
	regionHighCorner[1] = std::min(_interface.range.right + _densityTarget.cutoff, regionHighCorner[1]);

	double pi = 3.14159265;
	double rad1 = _densityTarget.cutoff;
	double rad2 = rad1*rad1;
	double rad3 = rad1*rad2;
	double hsv = 2./3.*pi*rad3;  // half sphere volume
	double hsv_inv = 1. / hsv;  // reciprocal half sphere volume
	std::vector<uint32_t> delList;  // list of particles to delete

	// the author intends to copy this iterator
	const auto begin = particleContainer->regionIterator(regionLowCorner, regionHighCorner,
														 ParticleIterator::ALL_CELLS);  // over all cell types
	for (auto it1 = begin; it1.isValid(); ++it1) {
		uint32_t pid1 = it1->getID();
		uint32_t numParticles_hsv = 0;
		double r1[3];
		r1[0] = it1->r(0);
		r1[1] = it1->r(1);
		r1[2] = it1->r(2);

		// only consider particles within interface range
		if(INTT_VAPOR_LIQUID == _interface.type && r1[1] > _interface.range.right)
			continue;
		if(INTT_LIQUID_VAPOR == _interface.type && r1[1] < _interface.range.left)
			continue;

		for(auto it2 = begin; it2.isValid(); ++it2) {
			uint32_t pid2 = it2->getID();
			if(pid2 == pid1)
				continue;
			double r2[3];
			r2[0] = it2->r(0);
			r2[1] = it2->r(1);
			r2[2] = it2->r(2);
			if(r2[1] > r1[1])
				continue;
			double diff[3];
			diff[0] = r2[0] - r1[0];
			diff[1] = r2[1] - r1[1];
			diff[2] = r2[2] - r1[2];
			double dist = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
			if(dist <= rad2) {
				numParticles_hsv++;
			}
		}
		double rho = numParticles_hsv * hsv_inv;
		if(rho < (_densityTarget.value * _densityTarget.percent) ) {
			if(_compChange.enabled) {
				Component* comp_new = global_simulation->getEnsemble()->getComponent(_compChange.cid_ub-1);
				it1->setComponent(comp_new);
			}
			else
				delList.push_back(pid1);
		}
	}

	for(auto it = begin; it.isValid(); ++it) {
		uint32_t pid1 = it->getID();
		bool bFound = false;
		for(auto pid2:delList) {
			if(pid2 == pid1)
				bFound = true;
		}
		if(bFound) {
			particleContainer->deleteMolecule(it, false);
		}
	}

	// Perform action only once
	_bDone.afterForces = true;
}
