/*
 * TemperatureControl.cpp
 *
 *  Created on: 27.05.2015
 *      Author: mheinen
 */

#include "thermostats/TemperatureControl.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "utils/Region.h"
#include "utils/DynAlloc.h"
#include "utils/FileUtils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <ostream>

using namespace std;

// init static ID --> instance counting
unsigned short tec::ControlRegion::_nStaticID = 0;

// class tec::ControlRegion

tec::ControlRegion::ControlRegion( ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nComp,
                                   double* dTargetTemperature, double dTemperatureExponent, std::string strTransDirections,
                                   int nTemperatureControlType, unsigned long nStartAdjust, unsigned long nStopAdjust, unsigned long nAdjustFreq)
: CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
{
	// ID
	_nID = ++_nStaticID;

    // init subdivision
    _nNumSlabs = 0;
    _dSlabWidthInit = 0.;

    _nTargetComponentID = nComp;
    _dTargetTemperature[0] = dTargetTemperature[0];
    _dTargetTemperature[1] = dTargetTemperature[1];

    _dTemperatureExponent = dTemperatureExponent;

    // type of temperature control constant/gradient/adjust
    _nTemperatureControlType = nTemperatureControlType;

    // temperature control adjust
	_nStartAdjust = nStartAdjust;
	_nStopAdjust  = nStopAdjust;
	_nAdjustFreq  = nAdjustFreq;

	// init datastructure pointers
	this->InitDataStructurePointers();

    if(TCT_TEMPERATURE_ADJUST == _nTemperatureControlType)
    {
		double dAdjustFreq = (double) (_nAdjustFreq);
		double dSteps = (_nStopAdjust - _nStartAdjust) / dAdjustFreq;
		_dDeltaTemperatureAdjust = (dTargetTemperature[1] - dTargetTemperature[0]) / dSteps;
		_dTargetTemperatureActual = _dTargetTemperature[0];
    }
    else if(TCT_TEMPERATURE_GRADIENT_LOWER == _nTemperatureControlType)
    {
		double dAdjustFreq = (double) (_nAdjustFreq);
		double dSteps = (_nStopAdjust - _nStartAdjust) / dAdjustFreq;

		double dT_lower = min(dTargetTemperature[0], dTargetTemperature[1]);
		double dT_upper = max(dTargetTemperature[0], dTargetTemperature[1]);

		_dDeltaTemperatureAdjust = (dT_lower - dT_upper) / dSteps;
		_dTargetTemperatureActual = dT_upper;
    }
    else if(TCT_TEMPERATURE_GRADIENT_RAISE == _nTemperatureControlType)
    {
		double dAdjustFreq = (double) (_nAdjustFreq);
		double dSteps = (_nStopAdjust - _nStartAdjust) / dAdjustFreq;

		double dT_lower = min(dTargetTemperature[0], dTargetTemperature[1]);
		double dT_upper = max(dTargetTemperature[0], dTargetTemperature[1]);

		_dDeltaTemperatureAdjust = (dT_upper - dT_lower) / dSteps;
		_dTargetTemperatureActual = dT_lower;
    }
    else
    {
		_dDeltaTemperatureAdjust = 0.;
		_dTargetTemperatureActual = _dTargetTemperature[0];
    }

    // create accumulator object dependent on which translatoric directions should be thermostated (xyz)
    if(strTransDirections == "x")
    {
        _accumulator = new AccumulatorXQ();
        _nNumThermostatedTransDirections = 1;
    }
    else if(strTransDirections == "y")
    {
        _accumulator = new AccumulatorYQ();
        _nNumThermostatedTransDirections = 1;
    }
    else if(strTransDirections == "z")
    {
        _accumulator = new AccumulatorZQ();
        _nNumThermostatedTransDirections = 1;
    }
    else if(strTransDirections == "xy")
    {
        _accumulator = new AccumulatorXYQ();
        _nNumThermostatedTransDirections = 2;
    }
    else if(strTransDirections == "xz")
    {
        _accumulator = new AccumulatorXZQ();
        _nNumThermostatedTransDirections = 2;
    }
    else if(strTransDirections == "yz")
    {
        _accumulator = new AccumulatorYZQ();
        _nNumThermostatedTransDirections = 2;
    }
    else if(strTransDirections == "xyz")
    {
        _accumulator = new AccumulatorXYZQ();
        _nNumThermostatedTransDirections = 3;
    }
    else
        _accumulator = NULL;

}


tec::ControlRegion::~ControlRegion()
{

}

void tec::ControlRegion::InitDataStructurePointers()
{
    _nNumMoleculesLocal  = NULL;
    _nNumMoleculesGlobal = NULL;
    _nRotDOFLocal  = NULL;
    _nRotDOFGlobal = NULL;

    _d2EkinTransLocal  = NULL;
    _d2EkinTransGlobal = NULL;
	_d2EkinRotLocal  = NULL;
    _d2EkinRotGlobal = NULL;

	_dBetaTransGlobal = NULL;
	_dBetaRotGlobal = NULL;

	_dTargetTemperatureVec = NULL;

    // heat supply
	_nNumMoleculesSumLocal  = NULL;
	_nNumMoleculesSumGlobal = NULL;

	_dDelta2EkinTransSumLocal  = NULL;
	_dDelta2EkinTransSumGlobal = NULL;
	_dDelta2EkinRotSumLocal  = NULL;
	_dDelta2EkinRotSumGlobal = NULL;
}

void tec::ControlRegion::AllocateDataStructuresT()
{
	// allocate more slabs as initially needed as a reserve for the case that the
	// temperature control region grows during simulation
    double dBoxWidthY = _parent->GetDomain()->getGlobalLength(1);
    _nNumSlabsReserve = (unsigned int) ( ceil(dBoxWidthY / this->GetWidth(1) * _nNumSlabs) );

#ifndef NDEBUG
    cout << "this->GetWidth(1) = " << this->GetWidth(1) << endl;
    cout << "dBoxWidthY = " << dBoxWidthY << endl;
    cout << "_nNumSlabs = " << _nNumSlabs << endl;
    cout << "_nNumSlabsReserve = " << _nNumSlabsReserve << endl;
#endif

	// allocate temperature data structures
	AllocateUnsLongArray(_nNumMoleculesLocal, _nNumSlabsReserve);
	AllocateUnsLongArray(_nNumMoleculesGlobal, _nNumSlabsReserve);
	AllocateUnsLongArray(_nRotDOFLocal, _nNumSlabsReserve);
	AllocateUnsLongArray(_nRotDOFGlobal, _nNumSlabsReserve);

	AllocateDoubleArray(_d2EkinTransLocal, _nNumSlabsReserve);
	AllocateDoubleArray(_d2EkinTransGlobal, _nNumSlabsReserve);
	AllocateDoubleArray(_d2EkinRotLocal, _nNumSlabsReserve);
	AllocateDoubleArray(_d2EkinRotGlobal, _nNumSlabsReserve);

	AllocateDoubleArray(_dBetaTransGlobal, _nNumSlabsReserve);
	AllocateDoubleArray(_dBetaRotGlobal, _nNumSlabsReserve);

	// target temperature vector to maintain a temperature gradient
	AllocateDoubleArray(_dTargetTemperatureVec, _nNumSlabsReserve);
}

void tec::ControlRegion::InitDataStructuresT()
{
//	cout << "tec::ControlRegion::_nNumSlabsReserve = " << _nNumSlabsReserve << endl;

	// init temperature datastructures
    for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
    {
        _nNumMoleculesLocal[s]  = 0;
        _nNumMoleculesGlobal[s] = 0;
        _nRotDOFLocal[s]        = 0;
        _nRotDOFGlobal[s]       = 0;

        _d2EkinTransLocal[s]  = 0.;
        _d2EkinTransGlobal[s] = 0.;
        _d2EkinRotLocal[s]    = 0.;
        _d2EkinRotGlobal[s]   = 0.;

        _dBetaTransGlobal[s] = 0.;
        _dBetaRotGlobal[s]   = 0.;
    }

	// target temperature vector to maintain a temperature gradient
    if( TCT_TEMPERATURE_GRADIENT       == _nTemperatureControlType ||
    	TCT_TEMPERATURE_GRADIENT_LOWER == _nTemperatureControlType ||
    	TCT_TEMPERATURE_GRADIENT_RAISE == _nTemperatureControlType )
    	this->AdjustTemperatureGradient();
}

void tec::ControlRegion::AllocateDataStructuresDEkin()
{
	// allocate more slabs as initially needed as a reserve for the case that the
	// temperature control region grows during simulation
    double dBoxWidthY = _parent->GetDomain()->getGlobalLength(1);
    _nNumSlabsDEkinReserve = (unsigned int) ( ceil(dBoxWidthY / this->GetWidth(1) * _nNumSlabs) );

#ifndef NDEBUG
	cout << "_nNumSlabsDEkinReserve = " << _nNumSlabsDEkinReserve << endl;
#endif

	// allocate Delta Ekin data structures
	AllocateUnsLongArray(_nNumMoleculesSumLocal, _nNumSlabsDEkinReserve);
	AllocateUnsLongArray(_nNumMoleculesSumGlobal, _nNumSlabsDEkinReserve);

	AllocateDoubleArray(_dDelta2EkinTransSumLocal, _nNumSlabsDEkinReserve);
	AllocateDoubleArray(_dDelta2EkinTransSumGlobal, _nNumSlabsDEkinReserve);
	AllocateDoubleArray(_dDelta2EkinRotSumLocal, _nNumSlabsDEkinReserve);
	AllocateDoubleArray(_dDelta2EkinRotSumGlobal, _nNumSlabsDEkinReserve);
}

void tec::ControlRegion::InitDataStructuresDEkin()
{
	// init Delta Ekin datastructures
    for(unsigned int s = 0; s<_nNumSlabsDEkinReserve; ++s)
    {
        _nNumMoleculesSumLocal[s]  = 0;
        _nNumMoleculesSumGlobal[s]  = 0;

        _dDelta2EkinTransSumLocal[s]  = 0.;
        _dDelta2EkinTransSumGlobal[s]  = 0.;
		_dDelta2EkinRotSumLocal[s]  = 0.;
		_dDelta2EkinRotSumGlobal[s]  = 0.;
    }

	// write out files (headers only)
    this->WriteHeaderDeltaEkin();
}

void tec::ControlRegion::PrepareDataStructures()
{
	this->AllocateDataStructuresT();
	this->InitDataStructuresT();
	this->AllocateDataStructuresDEkin();
	this->InitDataStructuresDEkin();
}

void tec::ControlRegion::CalcGlobalValues(unsigned long simstep)
{
#ifdef ENABLE_MPI

    // ToDo: communicate _nNumSlabs is enough?? or have to communicate _nNumSlabsReserve (whole data structure)
    MPI_Allreduce( _nNumMoleculesLocal, _nNumMoleculesGlobal, _nNumSlabsReserve, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _nRotDOFLocal, _nRotDOFGlobal, _nNumSlabsReserve, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce( _d2EkinTransLocal, _d2EkinTransGlobal, _nNumSlabsReserve, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _d2EkinRotLocal, _d2EkinRotGlobal, _nNumSlabsReserve, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
    {
        _nNumMoleculesGlobal[s] = _nNumMoleculesLocal[s];
        _nRotDOFGlobal[s] = _nRotDOFLocal[s];

        _d2EkinTransGlobal[s] = _d2EkinTransLocal[s];
        _d2EkinRotGlobal[s] = _d2EkinRotLocal[s];
    }
#endif

    // calc betaTrans, betaRot
    double dTargetTemperature = _dTargetTemperature[0];

    switch(_nTemperatureControlType)
    {

    case TCT_CONSTANT_TEMPERATURE:

		for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
		{
			if( _nNumMoleculesGlobal[s] < 1 )
				_dBetaTransGlobal[s] = 1.;
			else
				_dBetaTransGlobal[s] = pow(_nNumThermostatedTransDirections * _nNumMoleculesGlobal[s] * dTargetTemperature / _d2EkinTransGlobal[s], _dTemperatureExponent);

			if( _nRotDOFGlobal[s] < 1 )
				_dBetaRotGlobal[s] = 1.;
			else
				_dBetaRotGlobal[s] = pow( _nRotDOFGlobal[s] * dTargetTemperature / _d2EkinRotGlobal[s], _dTemperatureExponent);
		}
		break;


    case TCT_TEMPERATURE_GRADIENT:
    case TCT_TEMPERATURE_GRADIENT_LOWER:
    case TCT_TEMPERATURE_GRADIENT_RAISE:

    	if( _nStartAdjust < simstep && _nStopAdjust >= simstep &&
    	    0 == simstep % _nAdjustFreq && TCT_TEMPERATURE_GRADIENT != _nTemperatureControlType)
    	{
    		_dTargetTemperatureActual += _dDeltaTemperatureAdjust;
    		this->AdjustTemperatureGradient();
    	}

		for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
		{
			if( _nNumMoleculesGlobal[s] < 1 )
				_dBetaTransGlobal[s] = 1.;
			else
				_dBetaTransGlobal[s] = pow(_nNumThermostatedTransDirections * _nNumMoleculesGlobal[s] * _dTargetTemperatureVec[s] / _d2EkinTransGlobal[s], _dTemperatureExponent);

			if( _nRotDOFGlobal[s] < 1 )
				_dBetaRotGlobal[s] = 1.;
			else
				_dBetaRotGlobal[s] = pow( _nRotDOFGlobal[s] * _dTargetTemperatureVec[s] / _d2EkinRotGlobal[s], _dTemperatureExponent);
		}
		break;

    case TCT_TEMPERATURE_ADJUST:

    	if( _nStartAdjust < simstep && _nStopAdjust >= simstep &&
    	    0 == simstep % _nAdjustFreq)
    	{
    		_dTargetTemperatureActual += _dDeltaTemperatureAdjust;
//    		cout << "_dTargetTemperatureActual = " << _dTargetTemperatureActual << endl;
    	}

		for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
		{
			if( _nNumMoleculesGlobal[s] < 1 )
				_dBetaTransGlobal[s] = 1.;
			else
				_dBetaTransGlobal[s] = pow(_nNumThermostatedTransDirections * _nNumMoleculesGlobal[s] * _dTargetTemperatureActual / _d2EkinTransGlobal[s], _dTemperatureExponent);

			if( _nRotDOFGlobal[s] < 1 )
				_dBetaRotGlobal[s] = 1.;
			else
				_dBetaRotGlobal[s] = pow( _nRotDOFGlobal[s] * _dTargetTemperatureActual / _d2EkinRotGlobal[s], _dTemperatureExponent);
		}
		break;

    default:

    	cout << "No valid TemperatureControl type! Programm exit...";
    	exit(-1);
    }


//    cout << "_nNumMoleculesGlobal[0] = " << _nNumMoleculesGlobal[0] << endl;
//    cout << "_dBetaTransGlobal[0] = " << _dBetaTransGlobal[0] << endl;
//    cout << "_dTargetTemperature = " << _dTargetTemperature << endl;
//    cout << "_d2EkinTransGlobal[0] = " << _d2EkinTransGlobal[0] << endl;
//
//    cout << "_nRotDOFGlobal[0] = " << _nRotDOFGlobal[0] << endl;
//    cout << "_dBetaRotGlobal[0] = " << _dBetaRotGlobal[0] << endl;
//    cout << "_d2EkinRotGlobal[0] = " << _d2EkinRotGlobal[0] << endl;

}

void tec::ControlRegion::PrepareSubdivision()
{
	// catch errors
	if( (0 == _nNumSlabs && 0. == _dSlabWidthInit) || (0 != _nNumSlabs && 0. != _dSlabWidthInit) )
	{
		global_log->error() << "ERROR in tec::ControlRegion::PrepareSubdivision(): Neither _dSlabWidthInit nor _nNumSlabs was set correctly! Programm exit..." << endl;
		exit(-1);
	}

    // calc slab width
	if(0. == _dSlabWidthInit)
    	_dSlabWidthInit = this->GetWidth(1) / ( (double)(_nNumSlabs) );

	this->UpdateSlabParameters();

    // heat supply data structure TODO: heat supply data structure --> independent number of slabs
//    _dSlabWidthDeltaEkin = _parent->GetDomain()->getGlobalLength(1) / ( (double)(_nNumSlabsDeltaEkin) );

}

void tec::ControlRegion::MeasureKineticEnergy(Molecule* mol)
{
    // check componentID
    if(mol->componentid()+1 != _nTargetComponentID && 0 != _nTargetComponentID)  // program intern componentID starts with 0
        return;

    // check if molecule inside control region
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    unsigned int nPosIndex;
    unsigned int nIndexMax = _nNumSlabs - 1;

    // calc position index
    double* dLowerCorner = this->GetLowerCorner();
    double dPosRelative = mol->r(1) - dLowerCorner[1];

    nPosIndex = (unsigned int) floor(dPosRelative / _dSlabWidth);

    // ignore outer (halo) molecules
    if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
        return;

    // sum up transl. kinetic energy (2x)
/*
    double vx = mol->v(0);
//    double vy = mol->v(1);
    double vz = mol->v(2);
    double m  = mol->mass();

//    _d2EkinTransLocal += m*(vx*vx + vy*vy);
    _d2EkinTransLocal[nPosIndex] += m*(vx*vx + vz*vz);
*/

    _d2EkinTransLocal[nPosIndex] += _accumulator->CalcKineticEnergyContribution(mol);

    // sum up rot. kinetic energy (2x)
    double dDummy = 0.;

    mol->calculate_mv2_Iw2(dDummy, _d2EkinRotLocal[nPosIndex] );

    // count num molecules
    _nNumMoleculesLocal[nPosIndex]++;

    // count rotational DOF
    _nRotDOFLocal[nPosIndex] += mol->component()->getRotationalDegreesOfFreedom();
}


void tec::ControlRegion::ControlTemperature(Molecule* mol)
{
	if(this->GetID() == 2 && mol->r(1) > 120.)
		int cid = mol->componentid();

    // check componentID
    if(mol->componentid()+1 != _nTargetComponentID && 0 != _nTargetComponentID)  // program intern componentID starts with 0
        return;

//    if(mol->componentid()+1 == 2)
//    	cout << "region = " << this->GetID() << ": _nTargetComponentID = " << _nTargetComponentID << endl;

	// check if molecule is inside
	for(uint8_t d=0; d<3; ++d)
		if( !(PositionIsInside(d, mol->r(d) ) ) ) return;

    unsigned int nPosIndex;
    unsigned int nIndexMax = _nNumSlabs - 1;

    // calc position index
    double* dLowerCorner = this->GetLowerCorner();
    double dPosRelative = mol->r(1) - dLowerCorner[1];

    nPosIndex = (unsigned int) floor(dPosRelative / _dSlabWidth);

    // ignore outer (halo) molecules
    if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
        return;

	// calc scale factors
    double vcorr = 2. - 1. / _dBetaTransGlobal[nPosIndex];
    double Dcorr = 2. - 1. / _dBetaRotGlobal[nPosIndex];

    // add kinetic energy
	double dAdded2EkinTrans = 0.;
	double dAdded2EkinRot = 0.;
	_accumulator->ScaleVelocityComponents(mol, vcorr, Dcorr, dAdded2EkinTrans, dAdded2EkinRot);
	_dDelta2EkinTransSumLocal[nPosIndex] += dAdded2EkinTrans;
	_dDelta2EkinRotSumLocal[nPosIndex] += dAdded2EkinRot;

    // count number of molecules, where energy was added
    _nNumMoleculesSumLocal[nPosIndex]++;
//    cout << "nPosIndex = " << nPosIndex << endl;
//    cout << "_nNumMoleculesSumLocal[nPosIndex] = " << _nNumMoleculesSumLocal[nPosIndex] << endl;

	// sample kinetic energy, added to whole system
	dynamic_cast<TemperatureControl*>(_parent)->SampleDeltaEkin(mol, dAdded2EkinTrans, dAdded2EkinRot);
}

void tec::ControlRegion::ResetLocalValues()
{
    // reset local values
    for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
    {
        _nNumMoleculesLocal[s] = 0;
        _nRotDOFLocal[s] = 0;

        _d2EkinTransLocal[s] = 0.;
        _d2EkinRotLocal[s] = 0.;

        _dBetaTransGlobal[s] = 1.;
        _dBetaRotGlobal[s] = 1.;
    }
}

void tec::ControlRegion::UpdateSlabParameters()
{
    double dWidth = this->GetWidth(1);

    _nNumSlabs = round(dWidth / _dSlabWidthInit);
    _dSlabWidth =  dWidth / ( (double)(_nNumSlabs) );

	// heat supply data structure TODO: heat supply data structure --> independent number of slabs
	_dSlabWidthDeltaEkin = _dSlabWidth;
	_nNumSlabsDeltaEkin = _nNumSlabs;

    // target temperature vector to maintain a temperature gradient
    if( TCT_TEMPERATURE_GRADIENT       == _nTemperatureControlType ||
        TCT_TEMPERATURE_GRADIENT_LOWER == _nTemperatureControlType ||
        TCT_TEMPERATURE_GRADIENT_RAISE == _nTemperatureControlType )
        	this->AdjustTemperatureGradient();
}


void tec::ControlRegion::CalcGlobalValuesDeltaEkin()
{

#ifdef ENABLE_MPI

    MPI_Reduce( _nNumMoleculesSumLocal,    _nNumMoleculesSumGlobal,    _nNumSlabsDEkinReserve, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( _dDelta2EkinTransSumLocal, _dDelta2EkinTransSumGlobal, _nNumSlabsDEkinReserve, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( _dDelta2EkinRotSumLocal,   _dDelta2EkinRotSumGlobal,   _nNumSlabsDEkinReserve, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);

#else
    for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
    {
        _dDelta2EkinTransSumGlobal[s] = _dDelta2EkinTransSumLocal[s];
		_dDelta2EkinRotSumGlobal[s]   = _dDelta2EkinRotSumLocal[s];
    }
#endif

#ifndef NDEBUG
	cout << "Rank: " << _parent->GetDomainDecomposition()->getRank() << ", region: " << this->GetID() << ":_nNumMoleculesSumLocal[0] = " << _nNumMoleculesSumLocal[0] << endl;
	cout << "Rank: " << _parent->GetDomainDecomposition()->getRank() << ", region: " << this->GetID() << ":_nNumMoleculesSumGlobal[0] = " << _nNumMoleculesSumGlobal[0] << endl;
#endif
}


void tec::ControlRegion::ResetValuesDeltaEkin()
{
    // reset local values
    for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
    {
        _nNumMoleculesSumLocal[s] = 0;
        _dDelta2EkinTransSumLocal[s] = 0.;
		_dDelta2EkinRotSumLocal[s] = 0.;
    }

}

void tec::ControlRegion::WriteHeaderDeltaEkin()
{
	// domain decomposition
    DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();

    // write header
    stringstream outputstreamDeltaEkin;
    stringstream outputstreamN;

    std::stringstream filenamestreamDeltaEkin;
    std::stringstream filenamestreamN;

    filenamestreamDeltaEkin << "TemperatureControl_dEkin_region" << this->GetID() << ".dat";
    filenamestreamN         << "TemperatureControl_N_region" << this->GetID() << ".dat";

    string strFilenameDeltaEkin = filenamestreamDeltaEkin.str();
    string strFilenameN = filenamestreamN.str();

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif

    outputstreamDeltaEkin << "         simstep";
    outputstreamN         << "         simstep";

    for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
    {
        outputstreamDeltaEkin << "        slab[" << std::setfill('0') << std::setw(4) << s << "]";
        outputstreamN         << "        slab[" << std::setfill('0') << std::setw(4) << s << "]";  // TODO: add width for slab number
    }
    outputstreamDeltaEkin << endl;
    outputstreamN         << endl;

    ofstream fileoutDeltaEkin(strFilenameDeltaEkin.c_str(), ios::out);
    fileoutDeltaEkin << outputstreamDeltaEkin.str();
    fileoutDeltaEkin.close();

    ofstream fileoutN(strFilenameN.c_str(), ios::out);
    fileoutN << outputstreamN.str();
    fileoutN.close();

#ifdef ENABLE_MPI
    }
#endif
}

void tec::ControlRegion::WriteDataDeltaEkin(unsigned long simstep)
{
	// domain decomposition
    DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();

    // calc global values
    this->CalcGlobalValuesDeltaEkin();

    // reset local values
    this->ResetValuesDeltaEkin();

    // writing .dat-files
    std::stringstream outputstreamDeltaEkin;
    std::stringstream outputstreamN;

    std::stringstream filenamestreamDeltaEkin;
    std::stringstream filenamestreamN;

    filenamestreamDeltaEkin << "TemperatureControl_dEkin_region" << this->GetID() << ".dat";
    filenamestreamN         << "TemperatureControl_N_region"     << this->GetID() << ".dat";

    string strFilenameDeltaEkin = filenamestreamDeltaEkin.str();
    string strFilenameN         = filenamestreamN.str();


    #ifdef ENABLE_MPI
        int rank = domainDecomp->getRank();
        // int numprocs = domainDecomp->getNumProcs();
        if (rank== 0)
        {
    #endif

            // simstep
            outputstreamDeltaEkin << std::setw(16) << simstep;
            outputstreamN         << std::setw(16) << simstep;

            // DeltaEkin data
            for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
            {
                outputstreamDeltaEkin << std::setw(16) << fixed << std::setprecision(3) << _dDelta2EkinTransSumGlobal[s];
                outputstreamDeltaEkin << std::setw(16) << fixed << std::setprecision(3) << _dDelta2EkinRotSumGlobal[s];
                outputstreamN         << std::setw(16) << fixed << std::setprecision(3) << _nNumMoleculesSumGlobal[s];
            }
            outputstreamDeltaEkin << endl;
            outputstreamN         << endl;

            ofstream fileoutDeltaEkin(strFilenameDeltaEkin.c_str(), ios::app);
            fileoutDeltaEkin << outputstreamDeltaEkin.str();
            fileoutDeltaEkin.close();

            ofstream fileoutN(strFilenameN.c_str(), ios::app);
            fileoutN << outputstreamN.str();
            fileoutN.close();

    #ifdef ENABLE_MPI
        }
    #endif
}

void tec::ControlRegion::AdjustTemperatureGradient()
{
	double T[2];
	T[0] = _dTargetTemperature[0];
	T[1] = _dTargetTemperature[1];

    if(TCT_TEMPERATURE_GRADIENT_RAISE == _nTemperatureControlType)
    {
    	if(T[1] > T[0])
    		T[1] = _dTargetTemperatureActual;
    	else
    		T[0] = _dTargetTemperatureActual;
    }
    else if(TCT_TEMPERATURE_GRADIENT_LOWER == _nTemperatureControlType)
    {
    	if(T[1] > T[0])
    		T[0] = _dTargetTemperatureActual;
    	else
    		T[1] = _dTargetTemperatureActual;
    }

	double dT = T[1] - T[0];
	double dSlabsMinusOne = (double)(_nNumSlabs-1);

//	cout << "--------------------------" << endl;
//	cout << "region: " << this->GetID() << endl;
	for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
	{
		_dTargetTemperatureVec[s] = T[0] + dT/dSlabsMinusOne*s;
//		cout << "[" << s << "]: " << dT << " | " << dT/dSlabsMinusOne*s << " | " << T[0] << " | " << _dTargetTemperatureVec[s] << "" << endl;

	}
//	cout << "--------------------------" << endl;
}


// class TemperatureControl

TemperatureControl::TemperatureControl(Domain* domain, DomainDecompBase* domainDecomp, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop)
: ControlInstance(domain, domainDecomp)
{
    // control frequency
    _nControlFreq = nControlFreq;

    // start/stop timestep
    _nStart = nStart;
    _nStop  = nStop;

    // init heat supply variables
    _bWriteDataDeltaEkin = false;
    _nNumSlabsDeltaEkin = 1;
    _nWriteFreqDeltaEkin = 1000;
}

TemperatureControl::~TemperatureControl()
{

}

void TemperatureControl::AddRegion(tec::ControlRegion* region)
{
    _vecControlRegions.push_back(region);
}

void TemperatureControl::MeasureKineticEnergy(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // measure drift in each control region
    std::vector<tec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->MeasureKineticEnergy(mol);
    }
}

void TemperatureControl::CalcGlobalValues(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<tec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->CalcGlobalValues(simstep);
    }
}


void TemperatureControl::ControlTemperature(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // control drift of all regions
    std::vector<tec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ControlTemperature(mol);

//        cout << *(*it);
    }
}

void TemperatureControl::InitControl(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // reset local values
    std::vector<tec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ResetLocalValues();
    }
}

void TemperatureControl::DoLoopsOverMolecules(ParticleContainer* particleContainer, unsigned long simstep)
{
    // respect start/stop
    if(this->GetStart() <= simstep && this->GetStop() > simstep)
    {
//        cout << "Thermostat ON!" << endl;

        Molecule* tM;

        // init temperature control
        this->InitControl(simstep);

        for( tM  = particleContainer->begin();
             tM != particleContainer->end();
             tM  = particleContainer->next() )
        {
            // measure kinetic energy
            this->MeasureKineticEnergy(tM, simstep);

//          cout << "id = " << tM->id() << ", (vx,vy,vz) = " << tM->v(0) << ", " << tM->v(1) << ", " << tM->v(2) << endl;
        }


        // calc global values
        this->CalcGlobalValues(simstep);


        for( tM  = particleContainer->begin();
             tM != particleContainer->end();
             tM  = particleContainer->next() )
        {
            // control temperature
            this->ControlTemperature(tM, simstep);

//          cout << "id = " << tM->id() << ", (vx,vy,vz) = " << tM->v(0) << ", " << tM->v(1) << ", " << tM->v(2) << endl;
        }


        // write out heat supply
        if(true == _bWriteDataDeltaEkin)
        	this->WriteDataDeltaEkin(simstep);
    }
}

void TemperatureControl::SetDeltaEkinParameters(const paramLineHeat &paramLine)
{
	double dBoxWidth = this->GetDomain()->getGlobalLength(1);

	_nWriteFreqDeltaEkin = paramLine.nWriteFreq;
	if(paramLine.slabsKey == "num")
	{
		_nNumSlabsDeltaEkin = atoi(paramLine.slabsVal.c_str() );
		double dNum = (double)(_nNumSlabsDeltaEkin);
		_dSlabWidthDeltaEkin = dBoxWidth / dNum;
	}
	else if(paramLine.slabsKey == "width")
	{
		double dSlabWidth = atof(paramLine.slabsVal.c_str() );
		_nNumSlabsDeltaEkin = round(dBoxWidth / dSlabWidth);
		double dNum = (double)(_nNumSlabsDeltaEkin);
		_dSlabWidthDeltaEkin = dBoxWidth / dNum;
	}
	_nWriteFreqDeltaEkinRegions = paramLine.nWriteFreqRegions;

	_bWriteDataDeltaEkin = true;

	// prepare data structures
	uint8_t nNumComponents = this->GetDomain()->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components
	std::vector<unsigned long> vn(_nNumSlabsDeltaEkin, 0);
	std::vector<double> vd(_nNumSlabsDeltaEkin, 0);
	for(uint8_t cid=0; cid<nNumComponents; cid++)
	{
		_nNumMoleculesSumLocal.push_back(vn);
		_nNumMoleculesSumGlobal.push_back(vn);
		_dDelta2EkinTransSumLocal.push_back(vd);
		_dDelta2EkinTransSumGlobal.push_back(vd);
		_dDelta2EkinRotSumLocal.push_back(vd);
		_dDelta2EkinRotSumGlobal.push_back(vd);
	}

	double dPos = 0.5*_dSlabWidthDeltaEkin;
	_dSlabCoords.resize(_nNumSlabsDeltaEkin);

	for(uint32_t s=0; s<_nNumSlabsDeltaEkin; s++)
	{
		_dSlabCoords.at(s) = dPos;
		dPos += _dSlabWidthDeltaEkin;
	}
}

void TemperatureControl::SampleDeltaEkin(Molecule* mol, const double &dAdded2EkinTrans, const double &dAdded2EkinRot)
{
	// calc position index for DeltaEkin profile
	double dPosY = mol->r(1);
	unsigned int nPosIndex = (unsigned int) floor(dPosY / _dSlabWidthDeltaEkin);
	unsigned int nIndexMax = _nNumSlabsDeltaEkin - 1;
//	cout << "nPosIndex = " << nPosIndex << endl;
//	cout << "nIndexMax = " << nIndexMax << endl;

	// ignore outer (halo) molecules
	if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
		return;

	_dDelta2EkinTransSumLocal.at(0).at(nPosIndex) += dAdded2EkinTrans;
	_dDelta2EkinRotSumLocal  .at(0).at(nPosIndex) += dAdded2EkinRot;
	_nNumMoleculesSumLocal   .at(0).at(nPosIndex)++;

	unsigned int cid = mol->componentid()+1;

	_dDelta2EkinTransSumLocal.at(cid).at(nPosIndex) += dAdded2EkinTrans;
	_dDelta2EkinRotSumLocal  .at(cid).at(nPosIndex) += dAdded2EkinRot;
	_nNumMoleculesSumLocal   .at(cid).at(nPosIndex)++;
}

void TemperatureControl::CalcGlobalValuesDeltaEkin()
{
	DomainDecompBase* domainDecomp = this->GetDomainDecomposition();
	uint8_t nNumComponents = this->GetDomain()->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

	uint8_t nNumDataStructures = 3;
	unsigned int nNumSlabs = _nNumMoleculesSumLocal.at(0).size();
	unsigned long nNumVals = nNumComponents * nNumDataStructures * nNumSlabs;

	// init coll. comm.
	domainDecomp->collCommInit(nNumVals);
	// append local values
	for(uint8_t cid = 0; cid < nNumComponents; ++cid)
	{
		for(auto& pos : _nNumMoleculesSumLocal.at(cid) )
		{
			domainDecomp->collCommAppendUnsLong(pos);
		}
		for(auto& pos : _dDelta2EkinTransSumLocal.at(cid) )
		{
			domainDecomp->collCommAppendDouble(pos);
		}
		for(auto& pos : _dDelta2EkinRotSumLocal.at(cid) )
		{
			domainDecomp->collCommAppendDouble(pos);
		}
	}
	// reduce
	domainDecomp->collCommAllreduceSum();
	// collect global values
	for(uint8_t cid = 0; cid < nNumComponents; ++cid)
	{
		for(auto& pos : _nNumMoleculesSumGlobal.at(cid) )
		{
			pos = domainDecomp->collCommGetUnsLong();
		}
		for(auto& pos : _dDelta2EkinTransSumGlobal.at(cid) )
		{
			pos = domainDecomp->collCommGetDouble();
		}
		for(auto& pos : _dDelta2EkinRotSumGlobal.at(cid) )
		{
			pos = domainDecomp->collCommGetDouble();
		}
	}
	// finalize
	domainDecomp->collCommFinalize();
}

void TemperatureControl::ResetLocalValuesDeltaEkin()
{
	uint8_t nNumComponents = this->GetDomain()->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

	// append local values
	for(uint8_t cid = 0; cid < nNumComponents; ++cid)
	{
		std::fill(_nNumMoleculesSumLocal.at(cid).begin(), _nNumMoleculesSumLocal.at(cid).end(), 0);
		std::fill(_dDelta2EkinTransSumLocal.at(cid).begin(), _dDelta2EkinTransSumLocal.at(cid).end(), 0);
		std::fill(_dDelta2EkinRotSumLocal.at(cid).begin(), _dDelta2EkinRotSumLocal.at(cid).end(), 0);
	}
}

void TemperatureControl::WriteDataDeltaEkin(unsigned long simstep)
{
	if(simstep % _nWriteFreqDeltaEkinRegions == 0)
	{
		// reset local values
		std::vector<tec::ControlRegion*>::iterator it;

		for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
		{
			(*it)->WriteDataDeltaEkin(simstep);
		}
	}

	if(simstep % _nWriteFreqDeltaEkin == 0)
	{
		// calc global values
		this->CalcGlobalValuesDeltaEkin();

		// reset local values
//		this->ResetLocalValuesDeltaEkin();

		// writing .dat-files
		std::stringstream outputstream;
		std::stringstream filenamestream;
		filenamestream << "TemperatureControl_dEkin" << "_TS" << fill_width('0', 9) << simstep << ".dat";  // << "_r" << this->GetDomainDecomposition()->getRank();

#ifdef ENABLE_MPI
		int rank = this->GetDomainDecomposition()->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank == 0)
		{
#endif

		uint8_t nNumComponents = this->GetDomain()->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

		// header
		outputstream << "                     pos";

		// temperature/density
		for(uint8_t c = 0; c < nNumComponents; ++c)
		{
			outputstream << "              U_trans[" << c << "]";
			outputstream << "                U_rot[" << c << "]";
		}
		outputstream << endl;

		// data
		for(unsigned int s = 0; s < _nNumSlabsDeltaEkin; ++s)
		{
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dSlabCoords.at(s);

			for(uint8_t c = 0; c < nNumComponents; ++c)
			{
				// kin. energy componentwise
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDelta2EkinTransSumGlobal.at(c).at(s);
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDelta2EkinRotSumGlobal  .at(c).at(s);
			}
			outputstream << endl;
		}

		// Datei zum schreiben oeffnen, daten schreiben
		ofstream fileout(filenamestream.str().c_str(), ios::out);
		fileout << outputstream.str();
		fileout.close();

#ifdef ENABLE_MPI
		}
#endif
	}
}


void TemperatureControl::PrepareRegionDataStructures()
{
    std::vector<tec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->PrepareDataStructures();
    }
}

void TemperatureControl::PrepareRegionSubdivisions()
{
    std::vector<tec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->PrepareSubdivision();
    }
}

