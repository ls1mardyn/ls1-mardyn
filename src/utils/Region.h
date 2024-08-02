/*
 * Region.h
 *
 *  Created on: 19.08.2016
 *      Author: mheinen
 */

#pragma once

#include "utils/ObserverBase.h"
#include <string>
#include <ostream>
#include <vector>
#include <array>
#include <cstdint>
#include <cmath>

enum SubdivisionOption
{
	SDOPT_UNKNOWN = 0,
	SDOPT_BY_NUM_SLABS = 1,
	SDOPT_BY_SLAB_WIDTH = 2,
};

class Domain;
class DomainDecompBase;

class ControlInstance
{
public:
	ControlInstance() {}
	virtual ~ControlInstance() {}
	virtual std::string getShortName() = 0;
};

// class Region
class Region
{
protected:
	Region(ControlInstance* parent);
	virtual ~Region();

public:
	unsigned short GetID() {return _nID;}
	int GetType() {return _nType;}
	ControlInstance* GetParent() {return _parent;}
	virtual void Print(std::ostream& os) = 0;

protected:
	ControlInstance* _parent;
	unsigned short _nID;
private:
	int _nType;

};  // class Region


// class CuboidRegion

class CuboidRegion : public Region
{
public:
	CuboidRegion(ControlInstance* parent);
	CuboidRegion(ControlInstance* parent, double dLC[3], double dUC[3] );
	virtual ~CuboidRegion();

	std::array<double,3> GetLowerCorner() {return _dLowerCorner;}
	std::array<double,3> GetUpperCorner() {return _dUpperCorner;}
	void SetLowerCorner(std::array<double,3> dLC) { _dLowerCorner = dLC; }
	void SetUpperCorner(std::array<double,3> dUC) { _dUpperCorner = dUC; }
	double GetWidth(const uint16_t& nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
	void GetRange(const uint16_t& nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner.at(nDim); dRangeEnd = _dUpperCorner.at(nDim);}
	bool PositionIsInside(const uint16_t& nDim, const double& dPos) {return (dPos > _dLowerCorner.at(nDim) ) && (dPos < _dUpperCorner.at(nDim) );}
	bool PositionIsInside(double* dPos) {
		if		( !(dPos[0] > _dLowerCorner.at(0) && dPos[0] < _dUpperCorner.at(0) ) ) return false;
		else if ( !(dPos[1] > _dLowerCorner.at(1) && dPos[1] < _dUpperCorner.at(1) ) ) return false;
		else if	( !(dPos[2] > _dLowerCorner.at(2) && dPos[2] < _dUpperCorner.at(2) ) ) return false;
		else	return true;
	}
	virtual void Print(std::ostream& os)
	{
		os << "----------------------------------------------------------------" << std::endl;
		os << "ID: " << _nID << std::endl;
		os << "width: " << this->GetWidth(0) << " " << this->GetWidth(1) << " " << this->GetWidth(2) << std::endl;
		std::array<double,3> lc = this->GetLowerCorner();
		std::array<double,3> uc = this->GetUpperCorner();
		os << "lowerCorner: " << lc[0] << " " << lc[1] << " " << lc[2] << std::endl;
		os << "upperCorner: " << uc[0] << " " << uc[1] << " " << uc[2] << std::endl;
		os << "----------------------------------------------------------------" << std::endl;
	}
	double GetVolume()
	{
		double dVolume = 1.;
		for(uint8_t dim=0; dim<3; ++dim)
			dVolume *= this->GetWidth(dim);
		return dVolume;
	}

protected:
	std::array<double,3> _dLowerCorner;
	std::array<double,3> _dUpperCorner;

	int _nSubdivisionOpt;

};  // class CuboidRegion

// class CuboidRegionObs

class CuboidRegionObs : public CuboidRegion, public ObserverBase
{
public:
	CuboidRegionObs(ControlInstance* parent);
	CuboidRegionObs(ControlInstance* parent, double dLC[3], double dUC[3] );
	virtual ~CuboidRegionObs();

	// ObserverBase methods
	void update(SubjectBase* subject) override;
	void PrepareAsObserver(const std::vector<uint32_t>& refCoords);

private:
	// observer
	std::array<double,6> _dDistToRefCoords;
	std::array<bool,6> _bMaskMidpointLeft;
	std::array<bool,6> _bMaskMidpointRight;

};  // class CuboidRegionObs


// class SphericalRegion
class SphericalRegion : public Region
{
public:
	SphericalRegion(ControlInstance* parent);
	SphericalRegion(ControlInstance* parent, double dCtr[3], double radius);
	virtual ~SphericalRegion();

	std::array<double,3> GetCenter() {return _dCenter;}
	double GetRadius() {return _dRadius;}
	double GetRadiusSquared() {return _dRadiusSquared;}
	void SetCenter(std::array<double,3> dCtr) { _dCenter = dCtr; }
	void SetRadius(double radius) { _dRadius = radius; _dRadiusSquared = std::pow(radius, 2);}
	// double GetWidth(const uint16_t& nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
	// void GetRange(const uint16_t& nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner.at(nDim); dRangeEnd = _dUpperCorner.at(nDim);}
	// bool PositionIsInside(const uint16_t& nDim, const double& dPos) {return (dPos > _dLowerCorner.at(nDim) ) && (dPos < _dUpperCorner.at(nDim) );}
	bool PositionIsInside(double* dPos) {
		double distanceFromCenterSquared = 0;
		for(int i = 0; i<3; i++){
			distanceFromCenterSquared += dPos[i]-_dCenter[i];
		}
		if(distanceFromCenterSquared > _dRadiusSquared) return false;
		else return true;
	}
	virtual void Print(std::ostream& os)
	{
		os << "----------------------------------------------------------------" << std::endl;
		os << "ID: " << _nID << std::endl;
		os << "regionType: SphericalRegion";
		os << "radius: " << this->GetRadius() << std::endl;
		std::array<double,3> center = this->GetCenter();
		os << "center: " << center[0] << " " << center[1] << " " << center[2] << std::endl;
		os << "----------------------------------------------------------------" << std::endl;
	}
	double GetVolume()
	{
		return 4.18879020479 * std::pow(this->GetRadius(),3);  // pi*4/3 = 4.18879020479
	}

protected:
	std::array<double,3> _dCenter;
	double _dRadius;
	double _dRadiusSquared;
};  // class SphericalRegion








// class SphericalRegionObs

class SphericalRegionObs : public SphericalRegion, public ObserverBase
{
public:
	SphericalRegionObs(ControlInstance* parent);
	SphericalRegionObs(ControlInstance* parent, double dCtr[3], double radius);
	virtual ~SphericalRegionObs();

	// ObserverBase methods
	void update(SubjectBase* subject) override;
	void PrepareAsObserver(const std::vector<uint32_t>& refCoords);

private:
	//#TODO: do I need the following:?
	// // observer
	// std::array<double,6> _dDistToRefCoords;
	// std::array<bool,6> _bMaskMidpointLeft;
	// std::array<bool,6> _bMaskMidpointRight;
};  // class SphericalRegionObs


// class SphereComplementRegion
class SphereComplementRegion : public CuboidRegion, public SphericalRegion
{
public:
	SphereComplementRegion(ControlInstance* parent);
	SphereComplementRegion(ControlInstance* parent, double dLC[3], double dUC[3], double dCtr[3], double radius);
	virtual ~SphereComplementRegion();

	// std::array<double,3> GetLowerCorner() {return _dLowerCorner;} //should be inherited!?
	// std::array<double,3> GetUpperCorner() {return _dUpperCorner;}
	// std::array<double,3> GetCenter() {return _dCenter;}
	// double GetRadius() {return _dRadius;}
	// double GetRadiusSquared() { return _dRadiusSquared; }
	// void SetLowerCorner(std::array<double,3> dLC) { _dLowerCorner = dLC; }
	// void SetUpperCorner(std::array<double,3> dUC) { _dUpperCorner = dUC; }
	// void SetCenter(std::array<double,3> dCtr) { _dCenter = dCtr; }
	// void SetRadius(double radius) { _dRadius = radius; _dRadiusSquared = std::pow(radius, 2);}

	// double GetWidth(const uint16_t& nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
	// void GetRange(const uint16_t& nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner.at(nDim); dRangeEnd = _dUpperCorner.at(nDim);}
	// bool PositionIsInside(const uint16_t& nDim, const double& dPos) {return (dPos > _dLowerCorner.at(nDim) ) && (dPos < _dUpperCorner.at(nDim) );}

	bool PositionIsInside(double* dPos) {
		return (CuboidRegion::PositionIsInside(dPos) && !SphericalRegion::PositionIsInside(dPos));
	}
	virtual void Print(std::ostream& os)
	{
		os << "----------------------------------------------------------------" << std::endl;
		os << "ID: " << SphericalRegion::_nID << std::endl;
		os << "regionType: SphereComplementRegion (Box - Sphere)";
		std::array<double,3> lc = this->GetLowerCorner();
		std::array<double,3> uc = this->GetUpperCorner();
		os << "lowerCorner: " << lc[0] << " " << lc[1] << " " << lc[2] << std::endl;
		os << "upperCorner: " << uc[0] << " " << uc[1] << " " << uc[2] << std::endl;
		os << "radiusRemoved: " << this->GetRadius() << " " << std::endl;
		os << "----------------------------------------------------------------" << std::endl;
	}
	double GetVolume()
	{
		return (CuboidRegion::GetVolume() - SphericalRegion::GetVolume());
	}

protected:
	std::array<double,3> _dLowerCorner;
	std::array<double,3> _dUpperCorner;
	std::array<double,3> _dCenter;
	double _dRadius;
	double _dRadiusSquared;
};  // class SphereComplementRegion


// class SphereComplementRegionObs

class SphereComplementRegionObs : public SphereComplementRegion, public ObserverBase
{
public:
	SphereComplementRegionObs(ControlInstance* parent);
	SphereComplementRegionObs(ControlInstance* parent, double dLC[3], double dUC[3], double dCtr[3], double dRadius);
	virtual ~SphereComplementRegionObs();

	// ObserverBase methods
	void update(SubjectBase* subject) override;
	void PrepareAsObserver(const std::vector<uint32_t>& refCoords);

private:
	// observer
	std::array<double,6> _dDistToRefCoords;
	std::array<bool,6> _bMaskMidpointLeft;
	std::array<bool,6> _bMaskMidpointRight;

};  // class SphereComplementRegionObs


std::ostream& operator<<( std::ostream& os, Region& region);

