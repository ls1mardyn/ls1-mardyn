/*
 * MEX.cpp
 *
 *  Created on: 16.02.2015
 *      Author: mheinen
 */

#include "NEMD/DriftControl.h"
#include "NEMD/DistControl.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
//#include <cmath>
//#include <list>
//#include <map>
#include "Domain.h"
//#include "NEMD/ParticleInsertion.h"
//#include "utils/Random.h"
#include "utils/xmlfileUnits.h"
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
//#include <iterator>  // std::advance

using namespace std;

// init static ID --> instance counting
unsigned short drc::ControlRegion::_nStaticID = 0;

// class ControlRegion

drc::ControlRegion::ControlRegion( ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3] )
: CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
{
	// ID
	_nID = ++_nStaticID;
}

drc::ControlRegion::ControlRegion( ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double dDirection[3], const double& dTargetVal)
: CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
{
	// ID
	_nID = ++_nStaticID;

	// prepare drift vector
	this->PrepareDriftVector(dDirection, dTargetVal);

    // target component ID
    _nTargetComponentID = nTargetComponentID;
}

void drc::ControlRegion::readXML(XMLfileUnits& xmlconfig)
{
	xmlconfig.getNodeValue("target/componentID", _nTargetComponentID);
	xmlconfig.getNodeValue("target/direction/x", _dDriftDirection[0]);
	xmlconfig.getNodeValue("target/direction/y", _dDriftDirection[1]);
	xmlconfig.getNodeValue("target/direction/z", _dDriftDirection[2]);
	xmlconfig.getNodeValue("target/value", _dDriftVelocityTargetVal);

	// prepare drift vector
	this->PrepareDriftVector();
}

drc::ControlRegion::~ControlRegion()
{

}


void drc::ControlRegion::CalcGlobalValues()
{
	// domain decomposition
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();

#ifdef ENABLE_MPI

    MPI_Allreduce( _dDriftVelocityLocal, _dDriftVelocityGlobal, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _dEkinLocal, _dEkinGlobal, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

    for(unsigned short d=0; d<3; ++d)
    {
        _dDriftVelocityGlobal[d] = _dDriftVelocityLocal[d];
    }
#endif

//    cout << "_dDriftVelocityLocal[0] = " << _dDriftVelocityLocal[0] << endl;
//    cout << "_dDriftVelocityLocal[1] = " << _dDriftVelocityLocal[1] << endl;
//    cout << "_dDriftVelocityLocal[2] = " << _dDriftVelocityLocal[2] << endl;
//
//    cout << "_dDriftVelocityGlobal[0] = " << _dDriftVelocityGlobal[0] << endl;
//    cout << "_dDriftVelocityGlobal[1] = " << _dDriftVelocityGlobal[1] << endl;
//    cout << "_dDriftVelocityGlobal[2] = " << _dDriftVelocityGlobal[2] << endl;

    // collective Communication
    domainDecomp.collCommInit(1);
    domainDecomp.collCommAppendUnsLong(_nNumMoleculesLocal);
    domainDecomp.collCommAllreduceSum();
    _nNumMoleculesGlobal = domainDecomp.collCommGetUnsLong();
    domainDecomp.collCommFinalize();

//    cout << "_nNumMoleculesLocal = " << _nNumMoleculesLocal << endl;
//    cout << "_nNumMoleculesGlobal = " << _nNumMoleculesGlobal << endl;

    // velocity sum has to be divided by number of molecules taken into account
    if(_nNumMoleculesGlobal > 0)
    {
        for(unsigned short d = 0; d<3; ++d)
        {
            _dDriftVelocityGlobal[d] = _dDriftVelocityGlobal[d] / (double) (_nNumMoleculesGlobal);
        }
    }

    // calculate vector to add to each molecule to obtain target drift velocity of system inside control region
    if(0. == _dDriftVelocityTargetVal)
    {
        for(unsigned short d = 0; d<3; ++d)
        {
            _dAddVector[d] = _dDriftVelocityGlobal[d] * -1.;
        }
    }
    else
    {
        for(unsigned short d = 0; d<3; ++d)
        {
            _dAddVector[d] = _dDriftVelocityTargetVector[d] - _dDriftVelocityGlobal[d];
        }
    }

//    cout << "_dDriftVelocityGlobal = " << _dDriftVelocityGlobal[0] << ", " << _dDriftVelocityGlobal[1] << ", " << _dDriftVelocityGlobal[2] << endl;
//    cout << "_dDriftVelocityTargetVector = " << _dDriftVelocityTargetVector[0] << ", " << _dDriftVelocityTargetVector[1] << ", " << _dDriftVelocityTargetVector[2] << endl;
//    cout << "_dAddVector = " << _dAddVector[0] << ", " << _dAddVector[1] << ", " << _dAddVector[2] << endl;

}

void drc::ControlRegion::MeasureDrift(Molecule* mol)
{
    // check if molecule inside control region
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    // sum up drift and kin. energy
    for(unsigned short d = 0; d<3; ++d)
    {
        double v = mol->v(d);

        _dDriftVelocityLocal[d] += v;
        _dEkinLocal[d] += v*v;  // * mol->mass();
    }

    // count num molecules
    _nNumMoleculesLocal++;
}

void drc::ControlRegion::CalcScaleFactor()
{
    double dEkin_new = 0;
    double dEkin_old = 0;
    double dDeltaEkin = 0;

    // sum up x, y, z dimension of kin. energy values
    for(unsigned short d = 0; d<3; ++d)
    {
        dEkin_old += _dEkinGlobal[d];

        // (v^2_target - v^2_drift) * N * m yields the kinetic energy added or removed from system
        // dDeltaEkin is positive, when energy added
        // dDeltaEkin is negative, when energy removed
        dDeltaEkin += _dDriftVelocityTargetVector[d] * _dDriftVelocityTargetVector[d] - _dDriftVelocityGlobal[d] * _dDriftVelocityGlobal[d];
    }

    dDeltaEkin = dDeltaEkin * _nNumMoleculesGlobal;
    dEkin_new = dEkin_old +dDeltaEkin;

//    cout << "dDeltaEkin = " << dDeltaEkin << endl;
//    cout << "dEkin_new = " << dEkin_new << endl;
//    cout << "dEkin_old = " << dEkin_old << endl;

    double f = 1.;

    if(dEkin_new > 0.)
        f = dEkin_old / dEkin_new;

    _dScaleFactor = sqrt(f);

    // no scaling
    _dScaleFactor = 1.;
}

void drc::ControlRegion::ControlDrift(Molecule* mol)
{
    if(_nNumMoleculesGlobal < 1)
        return;

    // check componentID
        if(mol->componentid()+1 != _nTargetComponentID)  // program intern componentID starts with 0
            return;

    // check if molecule is inside
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    // change velocity vector
    for(unsigned short d = 0; d<3; ++d)
    {
        double v = mol->v(d);

//        cout << "v[" << d << "] = " << v << endl;
//
//        cout << "_dAddVector[" << d << "] = " << _dAddVector[d] << endl;

        // control drift
        mol->setv( d, (v + _dAddVector[d]) * _dScaleFactor);

//        cout << "v[" << d << "] = " << mol->v(d) << endl;
    }
}

void drc::ControlRegion::ResetLocalValues()
{
    // reset local values
    for(unsigned short d = 0; d<3; ++d)
    {
        _dDriftVelocityLocal[d] = 0.;
        _dEkinLocal[d] = 0.;
    }

    _nNumMoleculesLocal = 0;

    // scale factor
    _dScaleFactor = 1.;
}

void drc::ControlRegion::PrepareDriftVector()
{
	// normalize direction vector
	double dLength = 0.;

	for(unsigned short d=0; d<3; ++d)
		dLength += _dDriftDirection[d] * _dDriftDirection[d];
	dLength = sqrt(dLength);

	// drift velocity target direction, target vector
	for(unsigned short d=0; d<3; ++d)
	{
		_dDriftDirection[d] = _dDriftDirection[d] / dLength;
		_dDriftVelocityTargetVector[d] = _dDriftDirection[d] * _dDriftVelocityTargetVal;
	}

	//    cout << "_dDriftDirection = " << _dDriftDirection[0] << ", " << _dDriftDirection[1] << ", " << _dDriftDirection[2] << endl;
	//    cout << "_dDriftVelocityTargetVector = " << _dDriftVelocityTargetVector[0] << ", " << _dDriftVelocityTargetVector[1] << ", " << _dDriftVelocityTargetVector[2] << endl;

	// vector to add to each molecule to obtain target drift velocity of system inside control region
	for(unsigned short d=0; d<3; ++d)
		_dAddVector[d] = 0.;
}

void drc::ControlRegion::PrepareDriftVector(const double dDirection[3], const double& dTargetVal)
{
	// normalize direction vector
	double dLength = 0.;

	for(unsigned short d=0; d<3; ++d)
	{
		dLength += dDirection[d] * dDirection[d];
	}
	dLength = sqrt(dLength);

	// drift velocity target direction, target vector
	for(unsigned short d=0; d<3; ++d)
	{
		_dDriftDirection[d] = dDirection[d] / dLength;
		_dDriftVelocityTargetVector[d] = _dDriftDirection[d] * dTargetVal;
	}

	//    cout << "_dDriftDirection = " << _dDriftDirection[0] << ", " << _dDriftDirection[1] << ", " << _dDriftDirection[2] << endl;
	//    cout << "_dDriftVelocityTargetVector = " << _dDriftVelocityTargetVector[0] << ", " << _dDriftVelocityTargetVector[1] << ", " << _dDriftVelocityTargetVector[2] << endl;

	// drift velocity target value
	_dDriftVelocityTargetVal = dTargetVal;

	// vector to add to each molecule to obtain target drift velocity of system inside control region
	for(unsigned short d=0; d<3; ++d)
	{
		_dAddVector[d] = 0.;
	}
}


// class DriftControl

DriftControl::DriftControl( Domain* domain, DomainDecompBase* domainDecomp)
		: ControlInstance(domain, domainDecomp)
{
}

DriftControl::DriftControl( Domain* domain, DomainDecompBase* domainDecomp, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop)
                          : ControlInstance(domain, domainDecomp)
{
    // control frequency
    _nControlFreq = nControlFreq;

    // start/stop timestep
    _nStart = nStart;
    _nStop  = nStop;
}

DriftControl::~DriftControl()
{

}

void DriftControl::readXML(XMLfileUnits& xmlconfig)
{
	// control
	xmlconfig.getNodeValue("control/start", _nStart);
	xmlconfig.getNodeValue("control/frequency", _nControlFreq);
	xmlconfig.getNodeValue("control/stop", _nStop);

	// add regions
	uint32_t numRegions = 0;
	uint32_t nRegID = 0;
	XMLfile::Query query = xmlconfig.query("region");
	numRegions = query.card();
	global_log->info() << "DriftControl: Number of sampling regions: " << numRegions << endl;
	if(numRegions < 1) {
		global_log->warning() << "DriftControl: No region parameters specified. Program exit ..." << endl;
		Simulation::exit(-1);
	}
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator outputRegionIter;
	for( outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++ )
	{
		xmlconfig.changecurrentnode( outputRegionIter );
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

		global_log->info() << "DriftControl->region["<<nRegID<<"]: lower corner: " << lc[0] << ", " << lc[1] << ", " << lc[2] << endl;
		global_log->info() << "DriftControl->region["<<nRegID<<"]: upper corner: " << uc[0] << ", " << uc[1] << ", " << uc[2] << endl;

		// add regions
		drc::ControlRegion* region = new drc::ControlRegion(this, lc, uc);
		this->AddRegion(region);

		// observer mechanism
		uint32_t refCoordsID[6] = {0, 0, 0, 0, 0, 0};
		xmlconfig.getNodeValue("coords/lcx@refcoordsID", refCoordsID[0]);
		xmlconfig.getNodeValue("coords/lcy@refcoordsID", refCoordsID[1]);
		xmlconfig.getNodeValue("coords/lcz@refcoordsID", refCoordsID[2]);
		xmlconfig.getNodeValue("coords/ucx@refcoordsID", refCoordsID[3]);
		xmlconfig.getNodeValue("coords/ucy@refcoordsID", refCoordsID[4]);
		xmlconfig.getNodeValue("coords/ucz@refcoordsID", refCoordsID[5]);

		bool bIsObserver = (refCoordsID[0]+refCoordsID[1]+refCoordsID[2]+refCoordsID[3]+refCoordsID[4]+refCoordsID[5]) > 0;

		if(true == bIsObserver)
		{
			region->PrepareAsObserver(refCoordsID);

			if(global_simulation->GetDistControl() != NULL)
				global_simulation->GetDistControl()->registerObserver(region);
			else
			{
				global_log->error() << "DriftControl->region["<<region->GetID()<<"]: Initialization of feature DistControl is needed before! Program exit..." << endl;
				exit(-1);
			}
		}

		region->readXML(xmlconfig);
		nRegID++;

	}  // for( outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++ )
}

void DriftControl::AddRegion(drc::ControlRegion* region)
{
    _vecControlRegions.push_back(region);
}

void DriftControl::MeasureDrift(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // measure drift in each control region
    std::vector<drc::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->MeasureDrift(mol);
    }
}

void DriftControl::CalcGlobalValues(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<drc::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->CalcGlobalValues();
    }
}

void DriftControl::CalcScaleFactors(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc scale factors for control regions
    std::vector<drc::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->CalcScaleFactor();
    }
}

void DriftControl::ControlDrift(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // control drift of all regions
    std::vector<drc::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ControlDrift(mol);
    }
}

void DriftControl::Init(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // reset local values
    std::vector<drc::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ResetLocalValues();
    }
}




