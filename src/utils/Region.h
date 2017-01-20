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
	ControlInstance(Domain* domain, DomainDecompBase* domainDecomp) : _domain(domain), _domainDecomp(domainDecomp) {}
	~ControlInstance() {}
    virtual std::string GetShortName() = 0;

    Domain* GetDomain() {return _domain;}
    DomainDecompBase* GetDomainDecomposition() {return _domainDecomp;}

private:
    Domain* _domain;
    DomainDecompBase* _domainDecomp;
};

class Region
{
protected:
    Region(ControlInstance* parent);
    ~Region();

public:
    unsigned short GetID() {return _nID;}
    int GetType() {return _nType;}
    ControlInstance* GetParent() {return _parent;}

protected:
    unsigned short _nID;
    ControlInstance* _parent;
private:
    int _nType;

};  // class Region


// class CuboidRegion

class CuboidRegion : public Region
{
public:
	CuboidRegion(ControlInstance* parent);
	CuboidRegion(ControlInstance* parent, double dLC[3], double dUC[3] );
    ~CuboidRegion();

    double GetLowerCorner(unsigned short nDim) {return _dLowerCorner[nDim];}
    double GetUpperCorner(unsigned short nDim) {return _dUpperCorner[nDim];}
    void GetLowerCorner(double* dLC) {dLC[0]=_dLowerCorner[0]; dLC[1]=_dLowerCorner[1]; dLC[2]=_dLowerCorner[2];}
    void GetUpperCorner(double* dLC) {dLC[0]=_dUpperCorner[0]; dLC[1]=_dUpperCorner[1]; dLC[2]=_dUpperCorner[2];}
    double* GetLowerCorner() {return _dLowerCorner;}
    double* GetUpperCorner() {return _dUpperCorner;}
    void SetLowerCorner(unsigned short nDim, double dVal) {_dLowerCorner[nDim] = dVal;}
    void SetUpperCorner(unsigned short nDim, double dVal) {_dUpperCorner[nDim] = dVal;}
    double GetWidth(unsigned short nDim) {return _dUpperCorner[nDim] - _dLowerCorner[nDim];}
    void GetRange(unsigned short nDim, double& dRangeBegin, double& dRangeEnd) {dRangeBegin = _dLowerCorner[nDim]; dRangeEnd = _dUpperCorner[nDim];}
    bool PositionIsInside(unsigned short nDim, double dPos) {return dPos > _dLowerCorner[nDim] && dPos < _dUpperCorner[nDim];}
    bool PositionIsInside(double* dPos)
    {
    	if		( !(dPos[0] > _dLowerCorner[0] && dPos[0] < _dUpperCorner[0]) ) return false;
    	else if ( !(dPos[1] > _dLowerCorner[1] && dPos[1] < _dUpperCorner[1]) ) return false;
    	else if	( !(dPos[2] > _dLowerCorner[2] && dPos[2] < _dUpperCorner[2]) ) return false;
    	else	return true;
    }

protected:
    double _dLowerCorner[3];
    double _dUpperCorner[3];

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
	virtual void set(double dMidpointLeft, double dMidpointRight);
	void PrepareAsObserver(const unsigned short refCoords[6]);

private:
    // observer
    double _dDistToRefCoords[6];
    bool _bMaskMidpointLeft[6];
    bool _bMaskMidpointRight[6];

};  // class CuboidRegionObs

#endif  // REGION_H_
