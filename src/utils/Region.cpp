/*
 * Region.cpp
 *
 *  Created on: 19.08.2016
 *      Author: mheinen
 */

#include "Region.h"
#include "utils/ObserverBase.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"

// class Region

Region::Region(ControlInstance* parent) : _parent(parent)
{
}

Region::~Region()
{
}


// class CuboidRegion

CuboidRegion::CuboidRegion(ControlInstance* parent) : Region(parent)
{
	for(unsigned short d=0; d<3; d++)
	{
		_dLowerCorner[d] = 0.;
		_dUpperCorner[d] = 1.;
	}
}

CuboidRegion::CuboidRegion(ControlInstance* parent, double dLC[3], double dUC[3] ) : Region(parent)
{
	for(unsigned short d=0; d<3; d++)
	{
		_dLowerCorner[d] = dLC[d];
		_dUpperCorner[d] = dUC[d];
	}
}

CuboidRegion::~CuboidRegion()
{
}

CuboidRegionObs::CuboidRegionObs(ControlInstance* parent) : CuboidRegion(parent)
{
}

CuboidRegionObs::CuboidRegionObs(ControlInstance* parent, double dLC[3], double dUC[3]) : CuboidRegion(parent, dLC, dUC)
{
}

CuboidRegionObs::~CuboidRegionObs()
{
}

// ObserverBase methods
void CuboidRegionObs::set(double dMidpointLeft, double dMidpointRight)
{
	// update lower corner
	_dLowerCorner[0] = _dDistToRefCoords[0] + _bMaskMidpointLeft[0] * dMidpointLeft + _bMaskMidpointRight[0] * dMidpointRight;
	_dLowerCorner[1] = _dDistToRefCoords[1] + _bMaskMidpointLeft[1] * dMidpointLeft + _bMaskMidpointRight[1] * dMidpointRight;
	_dLowerCorner[2] = _dDistToRefCoords[2] + _bMaskMidpointLeft[2] * dMidpointLeft + _bMaskMidpointRight[2] * dMidpointRight;

	// update upper corner
	_dUpperCorner[0] = _dDistToRefCoords[3] + _bMaskMidpointLeft[3] * dMidpointLeft + _bMaskMidpointRight[3] * dMidpointRight;
	_dUpperCorner[1] = _dDistToRefCoords[4] + _bMaskMidpointLeft[4] * dMidpointLeft + _bMaskMidpointRight[4] * dMidpointRight;
	_dUpperCorner[2] = _dDistToRefCoords[5] + _bMaskMidpointLeft[5] * dMidpointLeft + _bMaskMidpointRight[5] * dMidpointRight;
}

void CuboidRegionObs::PrepareAsObserver(const uint32_t refCoords[6])
{
	_dDistToRefCoords[0] = _dLowerCorner[0];
	_dDistToRefCoords[1] = _dLowerCorner[1];
	_dDistToRefCoords[2] = _dLowerCorner[2];
	_dDistToRefCoords[3] = _dUpperCorner[0];
	_dDistToRefCoords[4] = _dUpperCorner[1];
	_dDistToRefCoords[5] = _dUpperCorner[2];

	_bMaskMidpointLeft[0] = refCoords[0] == 1;
	_bMaskMidpointLeft[1] = refCoords[1] == 1;
	_bMaskMidpointLeft[2] = refCoords[2] == 1;
	_bMaskMidpointLeft[3] = refCoords[3] == 1;
	_bMaskMidpointLeft[4] = refCoords[4] == 1;
	_bMaskMidpointLeft[5] = refCoords[5] == 1;

	_bMaskMidpointRight[0] = refCoords[0] == 2;
	_bMaskMidpointRight[1] = refCoords[1] == 2;
	_bMaskMidpointRight[2] = refCoords[2] == 2;
	_bMaskMidpointRight[3] = refCoords[3] == 2;
	_bMaskMidpointRight[4] = refCoords[4] == 2;
	_bMaskMidpointRight[5] = refCoords[5] == 2;
}

std::ostream& operator<<( std::ostream& os, Region& region)
{
	region.Print(os);
	return os;
}
