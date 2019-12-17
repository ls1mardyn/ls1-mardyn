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
#include "plugins/NEMD/DistControl.h"

// class Region

Region::Region(ControlInstance* parent)
	: _parent(parent),
	_nID(0),
	_nType(0)
{
}

Region::~Region()
{
}


// class CuboidRegion

CuboidRegion::CuboidRegion(ControlInstance* parent)
	: Region(parent),
	_nSubdivisionOpt(SDOPT_UNKNOWN)
{
	for(auto d=0; d<3; d++) {
		_dLowerCorner.at(d) = 0.;
		_dUpperCorner.at(d) = 1.;
	}
}

CuboidRegion::CuboidRegion(ControlInstance* parent, double dLC[3], double dUC[3] )
	: Region(parent),
	_nSubdivisionOpt(SDOPT_UNKNOWN)
{
	for(auto d=0; d<3; d++) {
		_dLowerCorner.at(d) = dLC[d];
		_dUpperCorner.at(d) = dUC[d];
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
void CuboidRegionObs::update(SubjectBase* subject)
{
	DistControl* distControl = dynamic_cast<DistControl*>(subject);
	double dMidpointLeft = distControl->GetInterfaceMidLeft();
	double dMidpointRight = distControl->GetInterfaceMidRight();

	// update lower corner
	_dLowerCorner.at(0) = _dDistToRefCoords.at(0) + _bMaskMidpointLeft.at(0) * dMidpointLeft + _bMaskMidpointRight.at(0) * dMidpointRight;
	_dLowerCorner.at(1) = _dDistToRefCoords.at(1) + _bMaskMidpointLeft.at(1) * dMidpointLeft + _bMaskMidpointRight.at(1) * dMidpointRight;
	_dLowerCorner.at(2) = _dDistToRefCoords.at(2) + _bMaskMidpointLeft.at(2) * dMidpointLeft + _bMaskMidpointRight.at(2) * dMidpointRight;

	// update upper corner
	_dUpperCorner.at(0) = _dDistToRefCoords.at(3) + _bMaskMidpointLeft.at(3) * dMidpointLeft + _bMaskMidpointRight.at(3) * dMidpointRight;
	_dUpperCorner.at(1) = _dDistToRefCoords.at(4) + _bMaskMidpointLeft.at(4) * dMidpointLeft + _bMaskMidpointRight.at(4) * dMidpointRight;
	_dUpperCorner.at(2) = _dDistToRefCoords.at(5) + _bMaskMidpointLeft.at(5) * dMidpointLeft + _bMaskMidpointRight.at(5) * dMidpointRight;
}

void CuboidRegionObs::PrepareAsObserver(const std::vector<uint32_t>& refCoords)
{
	_dDistToRefCoords.at(0) = _dLowerCorner.at(0);
	_dDistToRefCoords.at(1) = _dLowerCorner.at(1);
	_dDistToRefCoords.at(2) = _dLowerCorner.at(2);
	_dDistToRefCoords.at(3) = _dUpperCorner.at(0);
	_dDistToRefCoords.at(4) = _dUpperCorner.at(1);
	_dDistToRefCoords.at(5) = _dUpperCorner.at(2);

	_bMaskMidpointLeft.at(0) = refCoords.at(0) == 1;
	_bMaskMidpointLeft.at(1) = refCoords.at(1) == 1;
	_bMaskMidpointLeft.at(2) = refCoords.at(2) == 1;
	_bMaskMidpointLeft.at(3) = refCoords.at(3) == 1;
	_bMaskMidpointLeft.at(4) = refCoords.at(4) == 1;
	_bMaskMidpointLeft.at(5) = refCoords.at(5) == 1;

	_bMaskMidpointRight.at(0) = refCoords.at(0) == 2;
	_bMaskMidpointRight.at(1) = refCoords.at(1) == 2;
	_bMaskMidpointRight.at(2) = refCoords.at(2) == 2;
	_bMaskMidpointRight.at(3) = refCoords.at(3) == 2;
	_bMaskMidpointRight.at(4) = refCoords.at(4) == 2;
	_bMaskMidpointRight.at(5) = refCoords.at(5) == 2;
}

std::ostream& operator<<( std::ostream& os, Region& region)
{
	region.Print(os);
	return os;
}
