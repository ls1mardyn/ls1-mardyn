/*
 * Region.h
 *
 *  Created on: 19.08.2016
 *      Author: mheinen
 */

#ifndef REGION_H_
#define REGION_H_

#include "utils/ObserverBase.h"
#include <string>
#include <ostream>
#include <vector>
#include <array>
#include <cstdint>

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

	double GetLowerCorner(const uint16_t& nDim) {return _dLowerCorner.at(nDim);}
	double GetUpperCorner(const uint16_t& nDim) {return _dUpperCorner.at(nDim);}
	double* GetLowerCorner() {return _dLowerCorner.data();}
	double* GetUpperCorner() {return _dUpperCorner.data();}
	void SetLowerCorner(const uint16_t& nDim, const double& dVal) {_dLowerCorner.at(nDim) = dVal;}
	void SetUpperCorner(const uint16_t& nDim, const double& dVal) {_dUpperCorner.at(nDim) = dVal;}
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
		double* lc = this->GetLowerCorner();
		double* uc = this->GetUpperCorner();
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

std::ostream& operator<<( std::ostream& os, Region& region);

#endif  // REGION_H_
