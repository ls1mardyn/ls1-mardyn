/*
 * MotionLimits.h
 *
 *  Created on: 09.03.2018
 *      Author: mheinen
 */

#ifndef MOTIONLIMITS_H_
#define MOTIONLIMITS_H_

#include "utils/Region.h"
#include <cstdint>
#include <unordered_map>
#include <vector>

struct ControlVars
{
	uint64_t start, stop;
};

struct MaxVals
{
	double F, F2, M, M2, v, v2, L, L2;
};

class XMLfileUnits;
class ParticleContainer;
class MotionLimits;
class RegionML : public CuboidRegion
{
public:
	RegionML(MotionLimits* parent, double dLowerCorner[3], double dUpperCorner[3] );
	virtual ~RegionML();

	void readXML(XMLfileUnits& xmlconfig);
	void preForces(ParticleContainer* particleCont);
	void postForces(ParticleContainer* particleCont);

private:
	std::unordered_map<uint32_t, MaxVals> _maxVals;
};

class MotionLimits : public ControlInstance
{
public:
	MotionLimits(Domain* domain, DomainDecompBase* domainDecomp);
	~MotionLimits();

	virtual std::string GetShortName() {return "MoLim";}
	void readXML(XMLfileUnits& xmlconfig);
	void AddRegion(RegionML* region);
	void preForces(ParticleContainer* particleCont);
	void postForces(ParticleContainer* particleCont);

private:
	ControlVars _control;
	std::vector<RegionML*> _regions;
};



#endif /* MOTIONLIMITS_H_ */
