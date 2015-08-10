/*
 * DensityControl.cpp
 *
 *  Created on: 29.05.2015
 *      Author: mheinen
 */

#include "NEMD/DensityControl.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
//#include <cmath>
//#include <list>
//#include <map>
#include "Domain.h"
//#include "NEMD/ParticleInsertion.h"
//#include "utils/Random.h"
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
//#include <iterator>  // std::advance

#include <cstdlib>

using namespace std;


// class ControlRegionD

ControlRegionD::ControlRegionD(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double& dTargetDensity)
{
    // region span
    for(unsigned short d=0; d<3; ++d)
    {
        _dLowerCorner[d] = dLowerCorner[d];
        _dUpperCorner[d] = dUpperCorner[d];
    }

    // calc volume
    _dInvertVolume = 1. / (this->GetWidth(0) * this->GetWidth(1) * this->GetWidth(2) );

    // target density
    _dTargetDensity = dTargetDensity;

    // target component ID (inert gas)
    _nTargetComponentID = nTargetComponentID;

    // store parent pointer
    _parent = parent;

    // init process relevance
    _bProcessIsRelevant = true;

    // init MPI
    this->InitMPI();
}


ControlRegionD::~ControlRegionD()
{

}

void ControlRegionD::InitMPI()
{
    DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();
//    int numprocs = domainDecomp->getNumProcs();
    int nRank = domainDecomp->getRank();

    // get number of relevant processes
    int nRelevant;
    int nNumRelevantGlobal;

    double bbMin[3];
    double bbMax[3];

    domainDecomp->getBoundingBoxMinMax(_parent->GetDomain(), bbMin, bbMax);

    if( (bbMin[1] < _dLowerCorner[1] && bbMax[1] < _dLowerCorner[1]) ||  (bbMin[1] > _dUpperCorner[1] && bbMax[1] > _dUpperCorner[1]) )
        nRelevant = 0;
    else
        nRelevant = 1;

    // collective Communication
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendInt(nRelevant);
    domainDecomp->collCommAllreduceSum();
    nNumRelevantGlobal = domainDecomp->collCommGetInt();
    domainDecomp->collCommFinalize();

    // cout << "nNumRelevantGlobal = " << nNumRelevantGlobal << endl;

    // allocate rank array
    _ranks = new int[nNumRelevantGlobal];

    int nUnregistered = 1;

    for(int r=nNumRelevantGlobal-1; r>=0; --r)
    {
        int nRankLocal = (nRank+1) * nRelevant * nUnregistered;
        int nRankGlobal;

       // cout << "nRankLocal = " << nRankLocal << endl;

        MPI_Allreduce( &nRankLocal, &nRankGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        _ranks[r] = nRankGlobal-1;

        if(nRankGlobal-1 == nRank)
            nUnregistered = 0;
    }

    // check / set process relevance
    bool bProcessIsRelevant = false;

    for(int r=0; r<nNumRelevantGlobal; ++r)
    {
      //  cout << "_ranks[" << r << "] = " << _ranks[r] << endl;

        if(nRank == _ranks[r])
            bProcessIsRelevant = true;
    }

    if(true == bProcessIsRelevant)
        _bProcessIsRelevant = true;
    else
        _bProcessIsRelevant = false;

    MPI_Group          group, newgroup;

    MPI_Comm_group(MPI_COMM_WORLD, &group);

    MPI_Group_incl(group, nNumRelevantGlobal, _ranks, &newgroup);
    MPI_Comm_create(MPI_COMM_WORLD, newgroup, &_newcomm);

    int groupRank;

    MPI_Group_rank(group, &groupRank);

  //  cout << "[" << nRank << "]: " << "groupRank = " << groupRank << endl;

    MPI_Group_rank(newgroup, &groupRank);

  //  cout << "[" << nRank << "]: " << "groupRank = " << groupRank << endl;

    int nRankLocal = nRank;
    int nRankGlobal;

    if (_newcomm != MPI_COMM_NULL)
    {
        MPI_Allreduce( &nRankLocal, &nRankGlobal, 1, MPI_INT, MPI_MAX, _newcomm);
  //      cout << "[" << nRank << "]: " << "NEW COMM! nRankGlobal = " << nRankGlobal << endl;
    }
    else
    {
 //       cout << "new comm is NULL!" << endl;
    }
}

void ControlRegionD::CalcGlobalValues(DomainDecompBase* domainDecomp)
{
#ifdef ENABLE_MPI
    if (_newcomm == MPI_COMM_NULL)
        return;

    MPI_Allreduce( &_nNumMoleculesLocal, &_nNumMoleculesGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, _newcomm);

    _dDensityGlobal = _nNumMoleculesGlobal * _dInvertVolume;
#else
    _dDensityGlobal = _nNumMoleculesLocal * _dInvertVolume;
#endif

}

void ControlRegionD::UpdateGlobalDensity(DomainDecompBase* domainDecomp, bool bDeleteMolecule)
{
    double dDensityLocal = _dDensityGlobal;
    int nDeletionsLocal = 0 + bDeleteMolecule;
    int nDeletionsGlobal = 0;

    // collective Communication
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendInt(nDeletionsLocal);
    domainDecomp->collCommAllreduceSum();
    nDeletionsGlobal = domainDecomp->collCommGetInt();
    domainDecomp->collCommFinalize();

    // update density
    dDensityLocal -= nDeletionsGlobal * _dInvertVolume;

    // collective Communication
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendDouble(dDensityLocal);
    domainDecomp->collCommAllreduceSum();
    _dDensityGlobal = domainDecomp->collCommGetDouble();
    domainDecomp->collCommFinalize();
}

void ControlRegionD::MeasureDensity(Molecule* mol, DomainDecompBase* domainDecomp)
{
    // check if molecule inside control region
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    // count num molecules
    _nNumMoleculesLocal++;
}


void ControlRegionD::ControlDensity(DomainDecompBase* domainDecomp, Molecule* mol, Simulation* simulation, bool& bDeleteMolecule)
{
    /*
    // DEBUG
    cout << "_dLowerCorner = " << _dLowerCorner[0] << ", " << _dLowerCorner[1] << ", " << _dLowerCorner[2] << endl;
    cout << "_dUpperCorner = " << _dUpperCorner[0] << ", " << _dUpperCorner[1] << ", " << _dUpperCorner[2] << endl;
    // DEBUG
*/

    /*
    // check component ID: if inert gas --> do nothing (return)
    if(mol->componentid()+1 == 3)  // program intern component ID starts with 0
        return;
*/

//    int nRank = domainDecomp->getRank();
//    double dPosY = mol->r(1);

    // check if molecule is inside
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

//    cout << "[" << nRank << "]: " << "id = " << mol->id() << ", y = " << dPosY;
//    cout << ", lcY = " << this->GetLowerCorner()[1] << ", ucY = " << this->GetUpperCorner()[1] << endl;

    // if(_dDensityGlobal <= _dTargetDensity)  // exchange identity, dapt velocity

//    cout << "[" << nRank << "]: " << "id = " << mol->id() << ", cid = " << mol->componentid()+1 << endl;

    if(mol->componentid()+1 != 3)  // program intern component ID starts with 0
    {
        // store old mass / velocity
        double m1 = mol->mass();

        // inertgas component
        vector<Component>& vComponents = *(simulation->getEnsemble()->components());
        vector<Component>::iterator cit;
        Component* comp3;  // inertgas

        cit = vComponents.begin();
        cit += 2;
        comp3 = &(*(cit));

        mol->setComponent(comp3);

        // calc new velocity
        double m2 = mol->mass();
        double frac_m1_m2 = m1 / m2;

        mol->scale_v(frac_m1_m2);
    }

    /*
    else  // remove molecule
    {
        bDeleteMolecule = true;
    }
    */

    /* initialize random seed: */
    srand (time(NULL) + mol->id() );

    /* generate secret number between 0 and 99999: */
    int nRand = rand() % 100000;

    double dRand = (double) (nRand / 100000.);

//    cout << "dRand = " << dRand << endl;
//
//    cout << "_dDensityGlobal = " << _dDensityGlobal << endl;
//    cout << "_dTargetDensity = " << _dTargetDensity << endl;

    double dPercentToTakeOut = (_dDensityGlobal - _dTargetDensity) / _dDensityGlobal;

//    cout << "dPercentToTakeOut = " << dPercentToTakeOut << endl;

    if(dPercentToTakeOut > 0. && dRand < abs(dPercentToTakeOut) )
        bDeleteMolecule = true;
    else
        bDeleteMolecule = false;
}

void ControlRegionD::ResetLocalValues()
{
    _nNumMoleculesLocal = 0;
}



// class DensityControl

DensityControl::DensityControl(DomainDecompBase* domainDecomp, Domain* domain, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop)
{
    // control frequency
    _nControlFreq = nControlFreq;

    // start/stop timestep
    _nStart = nStart;
    _nStop  = nStop;

    // store domain / domainDecomp pointer
    _domain = domain;
    _domainDecomp = domainDecomp;
}

DensityControl::~DensityControl()
{
}

void DensityControl::AddRegion(double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double& dTargetDensity)
{
    ControlRegionD region(this, dLowerCorner, dUpperCorner, nTargetComponentID, dTargetDensity);

    _vecControlRegions.push_back(region);

    // check / set process relevance
    std::vector<ControlRegionD>::iterator it;
    bool bProcessIsRelevant = false;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        if ( true == (*it).ProcessIsRelevant() );
            bProcessIsRelevant = true;
    }

    if(true == bProcessIsRelevant)
        _bProcessIsRelevant = true;
    else
        _bProcessIsRelevant = false;
}

void DensityControl::MeasureDensity(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // measure drift in each control region
    std::vector<ControlRegionD>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).MeasureDensity(mol, domainDecomp);
    }
}

void DensityControl::CalcGlobalValues(DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<ControlRegionD>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).CalcGlobalValues(domainDecomp);
    }
}

void DensityControl::UpdateGlobalDensities(DomainDecompBase* domainDecomp, unsigned long simstep, bool bDeleteMolecule)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<ControlRegionD>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).UpdateGlobalDensity(domainDecomp, bDeleteMolecule);
    }
}

void DensityControl::ControlDensity(DomainDecompBase* domainDecomp, Molecule* mol, Simulation* simulation, unsigned long simstep, bool& bDeleteMolecule)
{
    if(simstep % _nControlFreq != 0)
        return;

    // control drift of all regions
    std::vector<ControlRegionD>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).ControlDensity(domainDecomp, mol, simulation, bDeleteMolecule);
    }
}

void DensityControl::Init(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // reset local values
    std::vector<ControlRegionD>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).ResetLocalValues();

        // update volume
        (*it).UpdateVolume();
    }
}




