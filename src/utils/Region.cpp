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












// class SphericalRegion
SphericalRegion::SphericalRegion(ControlInstance* parent)
	: Region(parent)
{
	for(auto d=0; d<3; d++) {
		_dCenter.at(d) = 20;
	}
	_dRadius = 8;
	_dRadiusSquared = 64;
}

SphericalRegion::SphericalRegion(ControlInstance* parent, double dCtr[3], double dRadius)
	: Region(parent)
	{
	for(auto d=0; d<3; d++) {
		_dCenter.at(d) = dCtr[d];
	}
	_dRadius = dRadius;
	_dRadiusSquared = pow(dRadius,2);
}

SphericalRegion::~SphericalRegion()
{
}

SphericalRegionObs::SphericalRegionObs(ControlInstance* parent) : SphericalRegion(parent)
{
}

SphericalRegionObs::SphericalRegionObs(ControlInstance* parent, double dCtr[3], double dRadius ) : SphericalRegion(parent, dCtr, dRadius)
{
}

SphericalRegionObs::~SphericalRegionObs()
{
}

// ObserverBase methods
void SphericalRegionObs::update(SubjectBase* subject)
{
	// #TODO: something may have to happen here?!

	// DistControl* distControl = dynamic_cast<DistControl*>(subject);
	// double dMidpointLeft = distControl->GetInterfaceMidLeft();
	// double dMidpointRight = distControl->GetInterfaceMidRight();

	// // update lower corner
	// _dLowerCorner.at(0) = _dDistToRefCoords.at(0) + _bMaskMidpointLeft.at(0) * dMidpointLeft + _bMaskMidpointRight.at(0) * dMidpointRight;
	// _dLowerCorner.at(1) = _dDistToRefCoords.at(1) + _bMaskMidpointLeft.at(1) * dMidpointLeft + _bMaskMidpointRight.at(1) * dMidpointRight;
	// _dLowerCorner.at(2) = _dDistToRefCoords.at(2) + _bMaskMidpointLeft.at(2) * dMidpointLeft + _bMaskMidpointRight.at(2) * dMidpointRight;

	// // update upper corner
	// _dUpperCorner.at(0) = _dDistToRefCoords.at(3) + _bMaskMidpointLeft.at(3) * dMidpointLeft + _bMaskMidpointRight.at(3) * dMidpointRight;
	// _dUpperCorner.at(1) = _dDistToRefCoords.at(4) + _bMaskMidpointLeft.at(4) * dMidpointLeft + _bMaskMidpointRight.at(4) * dMidpointRight;
	// _dUpperCorner.at(2) = _dDistToRefCoords.at(5) + _bMaskMidpointLeft.at(5) * dMidpointLeft + _bMaskMidpointRight.at(5) * dMidpointRight;
}

void SphericalRegionObs::PrepareAsObserver(const std::vector<uint32_t>& refCoords)
{
	// #TODO: something may have to happen here?!
	
	// _dDistToRefCoords.at(0) = _dLowerCorner.at(0);
	// _dDistToRefCoords.at(1) = _dLowerCorner.at(1);
	// _dDistToRefCoords.at(2) = _dLowerCorner.at(2);
	// _dDistToRefCoords.at(3) = _dUpperCorner.at(0);
	// _dDistToRefCoords.at(4) = _dUpperCorner.at(1);
	// _dDistToRefCoords.at(5) = _dUpperCorner.at(2);

	// _bMaskMidpointLeft.at(0) = refCoords.at(0) == 1;
	// _bMaskMidpointLeft.at(1) = refCoords.at(1) == 1;
	// _bMaskMidpointLeft.at(2) = refCoords.at(2) == 1;
	// _bMaskMidpointLeft.at(3) = refCoords.at(3) == 1;
	// _bMaskMidpointLeft.at(4) = refCoords.at(4) == 1;
	// _bMaskMidpointLeft.at(5) = refCoords.at(5) == 1;

	// _bMaskMidpointRight.at(0) = refCoords.at(0) == 2;
	// _bMaskMidpointRight.at(1) = refCoords.at(1) == 2;
	// _bMaskMidpointRight.at(2) = refCoords.at(2) == 2;
	// _bMaskMidpointRight.at(3) = refCoords.at(3) == 2;
	// _bMaskMidpointRight.at(4) = refCoords.at(4) == 2;
	// _bMaskMidpointRight.at(5) = refCoords.at(5) == 2;
}







// class SphereComplementRegion


SphereComplementRegion::SphereComplementRegion(ControlInstance* parent)
	: CuboidRegion(parent), SphericalRegion(parent)
{
	for(auto d=0; d<3; d++) {
		_dLowerCorner.at(d) = 0.;
		_dUpperCorner.at(d) = 40.;
		_dCenter.at(d) = 20;
	}
	_dRadius = 8;
	_dRadiusSquared = 64;
}

SphereComplementRegion::SphereComplementRegion(ControlInstance* parent, double dLC[3], double dUC[3], double dCtr[3], double dRadius)
	: CuboidRegion(parent), SphericalRegion(parent)
	{
	for(auto d=0; d<3; d++) {
		_dLowerCorner.at(d) = dLC[d];
		_dUpperCorner.at(d) = dUC[d];
		_dCenter.at(d) = dCtr[d];
	}
	_dRadius = dRadius;
	_dRadiusSquared = std::pow(dRadius,2);
}

SphereComplementRegion::~SphereComplementRegion()
{
}

SphereComplementRegionObs::SphereComplementRegionObs(ControlInstance* parent) : SphereComplementRegion(parent)
{
}

SphereComplementRegionObs::SphereComplementRegionObs(ControlInstance* parent, double dLC[3], double dUC[3], double dCtr[3], double dRadius ) : SphereComplementRegion(parent, dLC, dUC, dCtr, dRadius)
{
}

SphereComplementRegionObs::~SphereComplementRegionObs()
{
}

// ObserverBase methods
void SphereComplementRegionObs::update(SubjectBase* subject)
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

void SphereComplementRegionObs::PrepareAsObserver(const std::vector<uint32_t>& refCoords)
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
// class SphereComplementRegion



std::ostream& operator<<( std::ostream& os, Region& region)
{
	region.Print(os);
	return os;
}
