/*
 * MotionLimits.cpp
 *
 *  Created on: 09.03.2018
 *      Author: mheinen
 */

#include "MotionLimits.h"
#include "Domain.h"

#include <cstdint>
#include <cmath>
#include <string>

RegionML::RegionML(MotionLimits* parent, double dLowerCorner[3], double dUpperCorner[3] )
	:
	CuboidRegion(parent, dLowerCorner, dUpperCorner)
{
	uint16_t numComponents = global_simulation->getEnsemble()->getComponents()->size();
	_maxVals.reserve(numComponents);
}
RegionML::~RegionML()
{
}

void RegionML::readXML(XMLfileUnits& xmlconfig)
{
	// add regions
	uint32_t numComponents = 0;
	XMLfile::Query query = xmlconfig.query("maxvals");
	numComponents = query.card();
	global_log->info() << "RegionML: Number of components: " << numComponents << endl;
	if(numComponents < 1) {
		global_log->error() << "RegionML: No component parameters specified. Program exit ..." << endl;
		Simulation::exit(-1);
	}
	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator nodeIter;
	for(nodeIter = query.begin(); nodeIter; nodeIter++)
	{
		xmlconfig.changecurrentnode(nodeIter);

		uint32_t cid_ub;  // unity based index
		xmlconfig.getNodeValue("@cid", cid_ub);
		MaxVals mv;
		xmlconfig.getNodeValue("force", mv.F);
		xmlconfig.getNodeValue("torque", mv.M);
		xmlconfig.getNodeValue("velo", mv.v);
		xmlconfig.getNodeValue("angmom", mv.L);
		mv.F2 = mv.F* mv.F;
		mv.M2 = mv.M* mv.M;
		mv.v2 = mv.v* mv.v;
		mv.L2 = mv.L* mv.L;
		_maxVals[cid_ub] = mv;
	}  // loop over maxvals nodes
}

void RegionML::preForces(ParticleContainer* particleCont)
{
	uint32_t cid_ub;  // unity based index
	double v2, L2;
	double dPos[3];

//	RegionParticleIterator begin = particleCont->iterateRegionBegin(_dLowerCorner, _dUpperCorner, ParticleIterator::ALL_CELLS);  // ONLY_INNER_AND_BOUNDARY
//	RegionParticleIterator end = particleCont->iterateRegionEnd();
//	for(auto&& pit=begin; pit != end; ++pit)

	ParticleIterator pit;
	for( pit  = particleCont->iteratorBegin();
		 pit != particleCont->iteratorEnd();
		 ++pit )
	{
		// check if particle is inside region
		for(uint16_t d=0; d<3; ++d)
			dPos[d] = pit->r(d);
		if(PositionIsInside(dPos) == false)
			continue;

		cid_ub = pit->componentid()+1;
		if(_maxVals.find(cid_ub) == _maxVals.end() )
			continue;
		MaxVals& mv = _maxVals[cid_ub];
		// VELOCITY
		v2 = pit->v2();
		if(v2 > mv.v2)
		{
			std::cout << "WARNING: max. VELOCITY exceeded! cid=" << cid_ub << ", v=" << sqrt(v2) << ", y=" << pit->r(1) << std::endl;
			double fac = sqrt(mv.v2/v2);
			double v[3];
			for (uint16_t d=0; d<3; ++d)
			{
				v[d] = pit->v(d) * fac;
				pit->setv(d, v[d]);
			}
		}
		// ANGULAR MOMENTUM
		L2 = pit->L2();
		if(L2 > mv.L2)
		{
			std::cout << "WARNING: max. ANGULAR MOMENTUM exceeded! cid=" << cid_ub << ", L=" << sqrt(L2) << ", y=" << pit->r(1) << std::endl;
			double fac = sqrt(mv.L2/L2);
			double L[3];
			for (uint16_t d=0; d<3; ++d)
			{
				L[d] = pit->D(d) * fac;
				pit->setD(d, L[d]);
			}
		}
	}
}

void RegionML::postForces(ParticleContainer* particleCont)
{
	uint32_t cid_ub;  // unity based index
	double F2, M2, v2, L2;
	double dPos[3];

//	RegionParticleIterator begin = particleCont->iterateRegionBegin(_dLowerCorner, _dUpperCorner, ParticleIterator::ALL_CELLS);  // ONLY_INNER_AND_BOUNDARY
//	RegionParticleIterator end = particleCont->iterateRegionEnd();
//	for(auto&& pit=begin; pit != end; ++pit)

	ParticleIterator pit;
	for( pit  = particleCont->iteratorBegin();
		 pit != particleCont->iteratorEnd();
		 ++pit )
	{
		// check if particle is inside region
		for(uint16_t d=0; d<3; ++d)
			dPos[d] = pit->r(d);
		if(PositionIsInside(dPos) == false)
			continue;

		cid_ub = pit->componentid()+1;
		if(_maxVals.find(cid_ub) == _maxVals.end() )
			continue;
		MaxVals& mv = _maxVals[cid_ub];
		// FORCE
		F2 = pit->F2();
		if(F2 > mv.F2)
		{
			std::cout << "WARNING: max. FORCE exceeded! cid=" << cid_ub << ", F=" << sqrt(F2) << ", y=" << pit->r(1) << std::endl;
			double fac = sqrt(mv.F2/F2);
			double F[3];
			for (uint16_t d=0; d<3; ++d)
				F[d] = pit->F(d) * fac;
			pit->setF(F);
		}
		// TORQUE
		M2 = pit->M2();
		if(M2 > mv.M2)
		{
			std::cout << "WARNING: max. TORQUE exceeded! cid=" << cid_ub << ", M=" << sqrt(M2) << ", y=" << pit->r(1) << std::endl;
			double fac = sqrt(mv.M2/M2);
			double M[3];
			for (uint16_t d=0; d<3; ++d)
				M[d] = pit->M(d) * fac;
			pit->setM(M);
		}
		// VELOCITY
		v2 = pit->v2();
		if(v2 > mv.v2)
		{
			std::cout << "WARNING: max. VELOCITY exceeded! cid=" << cid_ub << ", v=" << sqrt(v2) << ", y=" << pit->r(1) << std::endl;
			double fac = sqrt(mv.v2/v2);
			double v[3];
			for (uint16_t d=0; d<3; ++d)
			{
				v[d] = pit->v(d) * fac;
				pit->setv(d, v[d]);
			}
		}
		// ANGULAR MOMENTUM
		L2 = pit->L2();
		if(L2 > mv.L2)
		{
			std::cout << "WARNING: max. ANGULAR MOMENTUM exceeded! cid=" << cid_ub << ", L=" << sqrt(L2) << ", y=" << pit->r(1) << std::endl;
			double fac = sqrt(mv.L2/L2);
			double L[3];
			for (uint16_t d=0; d<3; ++d)
			{
				L[d] = pit->D(d) * fac;
				pit->setD(d, L[d]);
			}
		}
	}
}

// class MotionLimits : public ControlInstance
MotionLimits::MotionLimits(Domain* domain, DomainDecompBase* domainDecomp)
	:
	ControlInstance(domain, domainDecomp)
{
}
MotionLimits::~MotionLimits()
{
}

void MotionLimits::readXML(XMLfileUnits& xmlconfig)
{
	// control
	_control.start = 0; _control.stop = 1000000;
	xmlconfig.getNodeValue("control/start", _control.start);
	xmlconfig.getNodeValue("control/stop", _control.stop);

	// add regions
	uint32_t numRegions = 0;
	uint32_t nRegID = 0;
	XMLfile::Query query = xmlconfig.query("region");
	numRegions = query.card();
	global_log->info() << "MotionLimits: Number of control regions: " << numRegions << endl;
	if(numRegions < 1) {
		global_log->error() << "MotionLimits: No region parameters specified. Program exit ..." << endl;
		Simulation::exit(-1);
	}
	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator nodeIter;
	for(nodeIter = query.begin(); nodeIter; nodeIter++)
	{
		xmlconfig.changecurrentnode(nodeIter);
		double lc[3];
		double uc[3];
		std::string strVal[3];
		std::string strControlType;

		// coordinates
		xmlconfig.getNodeValue("coords/lcx", lc[0]);
		xmlconfig.getNodeValue("coords/lcy", lc[1]);
		xmlconfig.getNodeValue("coords/lcz", lc[2]);
		xmlconfig.getNodeValue("coords/ucx", strVal[0]);
		xmlconfig.getNodeValue("coords/ucy", strVal[1]);
		xmlconfig.getNodeValue("coords/ucz", strVal[2]);
		// read upper corner
		for(uint8_t d=0; d<3; ++d)
			uc[d] = (strVal[d] == "box") ? GetDomain()->getGlobalLength(d) : atof(strVal[d].c_str() );

		global_log->info() << "MotionLimits->region["<<nRegID<<"]: lower corner: " << lc[0] << ", " << lc[1] << ", " << lc[2] << endl;
		global_log->info() << "MotionLimits->region["<<nRegID<<"]: upper corner: " << uc[0] << ", " << uc[1] << ", " << uc[2] << endl;

		// add regions
		RegionML* region = new RegionML(this, lc, uc);
		this->AddRegion(region);

		// read region parameters from XML file
		region->readXML(xmlconfig);
		nRegID++;

	}  // loop over region nodes
}

void MotionLimits::AddRegion(RegionML* region)
{
	_regions.push_back(region);
}

void MotionLimits::preForces(ParticleContainer* particleCont)
{
	for(auto&& reg : _regions)
		reg->preForces(particleCont);
}
void MotionLimits::postForces(ParticleContainer* particleCont)
{
	for(auto&& reg : _regions)
		reg->postForces(particleCont);
}

