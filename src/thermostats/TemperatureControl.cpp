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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;


// class ControlRegionT

ControlRegionT::ControlRegionT(TemperatureControl* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nNumSlabs, unsigned int nComp, double dTargetTemperature, double dTemperatureExponent, std::string strTransDirections)
{
    // store parent pointer
    _parent = parent;

    // region span
    for(unsigned short d=0; d<3; ++d)
    {
        _dLowerCorner[d] = dLowerCorner[d];
        _dUpperCorner[d] = dUpperCorner[d];
    }

    _nTargetComponentID = nComp;
    _dTargetTemperature = dTargetTemperature;

    _dTemperatureExponent = dTemperatureExponent;

    // number of slabs
    _nNumSlabs = nNumSlabs;

    // calc slab width
    _dSlabWidthInit = this->GetWidth(1) / ( (double)(_nNumSlabs) );
    _dSlabWidth = _dSlabWidthInit;

    // init data structures
    this->Init();

    // create accumulator object dependent on which translatoric directions should be thermostated (xyz)
    if(strTransDirections == "x")
    {
        _accumulator = new AccumulatorX();
        _nNumThermostatedTransDirections = 1;
    }
    else if(strTransDirections == "y")
    {
        _accumulator = new AccumulatorY();
        _nNumThermostatedTransDirections = 1;
    }
    else if(strTransDirections == "z")
    {
        _accumulator = new AccumulatorZ();
        _nNumThermostatedTransDirections = 1;
    }
    else if(strTransDirections == "xy")
    {
        _accumulator = new AccumulatorXY();
        _nNumThermostatedTransDirections = 2;
    }
    else if(strTransDirections == "xz")
    {
        _accumulator = new AccumulatorXZ();
        _nNumThermostatedTransDirections = 2;
    }
    else if(strTransDirections == "yz")
    {
        _accumulator = new AccumulatorYZ();
        _nNumThermostatedTransDirections = 2;
    }
    else if(strTransDirections == "xyz")
    {
        _accumulator = new AccumulatorXYZ();
        _nNumThermostatedTransDirections = 3;
    }
    else
        _accumulator = NULL;

}


ControlRegionT::~ControlRegionT()
{

}

void ControlRegionT::Init()
{
    // allocate more slabs as initially needed as a reserve for the case that the
    // temperature control region grows during simulation
    unsigned int nNumSlabsReserve;
    double dBoxWidthY = _parent->GetDomain()->getGlobalLength(1);
    nNumSlabsReserve = (unsigned int) ( ceil(dBoxWidthY / this->GetWidth(1) ) );

    cout << "nNumSlabsReserve = " << nNumSlabsReserve << endl;

    _nNumMoleculesLocal  = new unsigned long[nNumSlabsReserve];
    _nNumMoleculesGlobal = new unsigned long[nNumSlabsReserve];
    _nRotDOFLocal  = new unsigned long[nNumSlabsReserve];
    _nRotDOFGlobal = new unsigned long[nNumSlabsReserve];

    _d2EkinTransLocal  = new double[nNumSlabsReserve];
    _d2EkinTransGlobal = new double[nNumSlabsReserve];
    _d2EkinRotLocal  = new double[nNumSlabsReserve];
    _d2EkinRotGlobal = new double[nNumSlabsReserve];

    _dBetaTransGlobal = new double[nNumSlabsReserve];
    _dBetaRotGlobal   = new double[nNumSlabsReserve];

    for(unsigned int s = 0; s<nNumSlabsReserve; ++s)
    {
        _nNumMoleculesLocal[s]  = 0;
        _nNumMoleculesGlobal[s] = 0;
        _nRotDOFLocal[s]  = 0;
        _nRotDOFGlobal[s] = 0;

        _d2EkinTransLocal[s]  = 0.;
        _d2EkinTransGlobal[s] = 0.;
        _d2EkinRotLocal[s]  = 0.;
        _d2EkinRotGlobal[s] = 0.;

        _dBetaTransGlobal[s] = 0.;
        _dBetaRotGlobal[s]   = 0.;
    }
}

void ControlRegionT::CalcGlobalValues(DomainDecompBase* domainDecomp)
{
#ifdef ENABLE_MPI

    MPI_Allreduce( _nNumMoleculesLocal, _nNumMoleculesGlobal, _nNumSlabs, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _nRotDOFLocal, _nRotDOFGlobal, _nNumSlabs, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce( _d2EkinTransLocal, _d2EkinTransGlobal, _nNumSlabs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _d2EkinRotLocal, _d2EkinRotGlobal, _nNumSlabs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    for(unsigned int s = 0; s<_nNumSlabs; ++s)
    {
        _nNumMoleculesGlobal[s] = _nNumMoleculesLocal[s];
        _nRotDOFGlobal[s] = _nRotDOFLocal[s];

        _d2EkinTransGlobal[s] = _d2EkinTransLocal[s];
        _d2EkinRotGlobal[s] = _d2EkinRotLocal[s];
    }
#endif

    // calc betaTrans, betaRot
    for(unsigned int s = 0; s<_nNumSlabs; ++s)
    {
        if( _nNumMoleculesGlobal[s] < 1 )
            _dBetaTransGlobal[s] = 1.;
        else
            _dBetaTransGlobal[s] = pow(_nNumThermostatedTransDirections * _nNumMoleculesGlobal[s] * _dTargetTemperature / _d2EkinTransGlobal[s], _dTemperatureExponent);

        if( _nRotDOFGlobal[s] < 1 )
            _dBetaRotGlobal[s] = 1.;
        else
            _dBetaRotGlobal[s] = pow( _nRotDOFGlobal[s] * _dTargetTemperature / _d2EkinRotGlobal[s], _dTemperatureExponent);
    }

//    cout << "_nNumMoleculesGlobal = " << _nNumMoleculesGlobal << endl;
//    cout << "_dBetaTransGlobal = " << _dBetaTransGlobal << endl;
//    cout << "_dTargetTemperature = " << _dTargetTemperature << endl;
//    cout << "_d2EkinRotGlobal = " << _d2EkinRotGlobal << endl;
//
//    cout << "_nRotDOFGlobal = " << _nRotDOFGlobal << endl;
//    cout << "_dBetaRotGlobal = " << _dBetaRotGlobal << endl;
//    cout << "_d2EkinRotGlobal = " << _d2EkinRotGlobal << endl;

}

void ControlRegionT::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp)
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


void ControlRegionT::ControlTemperature(Molecule* mol)
{
    // check componentID
    if(mol->componentid()+1 != _nTargetComponentID && 0 != _nTargetComponentID)  // program intern componentID starts with 0
        return;

    // check if molecule is inside
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

    if(_nNumMoleculesGlobal[nPosIndex] < 1)
        return;


    // scale velocity
    double vcorr = 2. - 1. / _dBetaTransGlobal[nPosIndex];
    double Dcorr = 2. - 1. / _dBetaRotGlobal[nPosIndex];

/*
    mol->setv(0, mol->v(0) * vcorr);
//    mol->setv(1, mol->v(1) * vcorr);
    mol->setv(2, mol->v(2) * vcorr);
*/

    _accumulator->ScaleVelocityComponents(mol, vcorr);

    mol->scale_D(Dcorr);
}

void ControlRegionT::ResetLocalValues()
{
    // reset local values
    for(unsigned int s = 0; s<_nNumSlabs; ++s)
    {
        _nNumMoleculesLocal[s] = 0;
        _nRotDOFLocal[s] = 0;

        _d2EkinTransLocal[s] = 0.;
        _d2EkinRotLocal[s] = 0.;

        _dBetaTransGlobal[s] = 1.;
        _dBetaRotGlobal[s] = 1.;
    }
}

void ControlRegionT::UpdateSlabParameters()
{
    double dWidth = this->GetWidth(1);
//    unsigned int nNumSlabsOld = _nNumSlabs;

    _nNumSlabs = round(dWidth / _dSlabWidthInit);
    _dSlabWidth =  dWidth / ( (double)(_nNumSlabs) );

/*
    // number of slabs cannot increase, otherwise data structures have to be reallocated
    if(_nNumSlabs > nNumSlabsOld)
    {
        _nNumSlabs = nNumSlabsOld;
        _dSlabWidth =  dWidth / ( (double)(_nNumSlabs) );

        cout << "WARNING! Temperature reason increased!" << endl;
    }
*/
}



// class TemperatureControl

TemperatureControl::TemperatureControl(Domain* domain, DomainDecompBase* domainDecomp, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop)
{
    // store pointer to Domain and DomainDecomposition
    _domain = domain;
    _domainDecomp = domainDecomp;

    // control frequency
    _nControlFreq = nControlFreq;

    // start/stop timestep
    _nStart = nStart;
    _nStop  = nStop;
}

TemperatureControl::~TemperatureControl()
{

}

void TemperatureControl::AddRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nNumSlabs, unsigned int nComp, double dTargetTemperature, double dTemperatureExponent, std::string strTransDirections)
{
    _vecControlRegions.push_back(ControlRegionT(this, dLowerCorner, dUpperCorner, nNumSlabs, nComp, dTargetTemperature, dTemperatureExponent, strTransDirections) );
}

void TemperatureControl::MeasureKineticEnergy(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // measure drift in each control region
    std::vector<ControlRegionT>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).MeasureKineticEnergy(mol, domainDecomp);
    }
}

void TemperatureControl::CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<ControlRegionT>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).CalcGlobalValues(domainDecomp);
    }
}


void TemperatureControl::ControlTemperature(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // control drift of all regions
    std::vector<ControlRegionT>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).ControlTemperature(mol);
    }
}

void TemperatureControl::Init(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // reset local values
    std::vector<ControlRegionT>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).ResetLocalValues();
    }
}

void TemperatureControl::DoLoopsOverMolecules(DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, unsigned long simstep)
{
    // respect start/stop
    if(this->GetStart() <= simstep && this->GetStop() > simstep)
    {
//        cout << "Thermostat ON!" << endl;

        Molecule* tM;

        // init temperature control
        this->Init(simstep);

        for( tM  = particleContainer->begin();
             tM != particleContainer->end();
             tM  = particleContainer->next() )
        {
            // measure kinetic energy
            this->MeasureKineticEnergy(tM, domainDecomposition, simstep);

//          cout << "id = " << tM->id() << ", (vx,vy,vz) = " << tM->v(0) << ", " << tM->v(1) << ", " << tM->v(2) << endl;
        }

        // calc global values
        this->CalcGlobalValues(domainDecomposition, simstep);


        for( tM  = particleContainer->begin();
             tM != particleContainer->end();
             tM  = particleContainer->next() )
        {
            // control temperature
            this->ControlTemperature(tM, simstep);

//          cout << "id = " << tM->id() << ", (vx,vy,vz) = " << tM->v(0) << ", " << tM->v(1) << ", " << tM->v(2) << endl;
        }
    }
}



