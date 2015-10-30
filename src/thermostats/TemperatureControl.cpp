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
#include <iomanip>

using namespace std;


// class ControlRegionT

ControlRegionT::ControlRegionT( TemperatureControl* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nNumSlabs, unsigned int nComp,
                                double dTargetTemperature, double dTemperatureExponent, std::string strTransDirections, unsigned short nRegionID, unsigned int nNumSlabsDeltaEkin )
{
    // store parent pointer
    _parent = parent;

    // region ID
    _nRegionID = nRegionID;

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

    // heat supply data structure
    _nNumSlabsDeltaEkin = nNumSlabsDeltaEkin;
    _dSlabWidthDeltaEkin = _parent->GetDomain()->getGlobalLength(1) / ( (double)(_nNumSlabsDeltaEkin) );

    // init data structures
    this->Init();

    // create accumulator object dependent on which translatoric directions should be thermostated (xyz)
/*
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
*/

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


ControlRegionT::~ControlRegionT()
{

}

void ControlRegionT::Init()
{
    // allocate more slabs as initially needed as a reserve for the case that the
    // temperature control region grows during simulation
    double dBoxWidthY = _parent->GetDomain()->getGlobalLength(1);
    _nNumSlabsReserve = (unsigned int) ( ceil(dBoxWidthY / this->GetWidth(1) * _nNumSlabs) );

#ifdef DEBUG
    cout << "_nNumSlabsReserve = " << _nNumSlabsReserve << endl;
#endif

    _nNumMoleculesLocal  = new unsigned long[_nNumSlabsReserve];
    _nNumMoleculesGlobal = new unsigned long[_nNumSlabsReserve];
    _nRotDOFLocal  = new unsigned long[_nNumSlabsReserve];
    _nRotDOFGlobal = new unsigned long[_nNumSlabsReserve];

    _d2EkinTransLocal  = new double[_nNumSlabsReserve];
    _d2EkinTransGlobal = new double[_nNumSlabsReserve];
    _d2EkinRotLocal  = new double[_nNumSlabsReserve];
    _d2EkinRotGlobal = new double[_nNumSlabsReserve];

    _dBetaTransGlobal = new double[_nNumSlabsReserve];
    _dBetaRotGlobal   = new double[_nNumSlabsReserve];

    for(unsigned int s = 0; s<_nNumSlabsReserve; ++s)
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


    // mheinen 2015-10-27 --> CALC_HEAT_SUPPLY

    _nNumMoleculesSumLocal = new unsigned long[_nNumSlabsDeltaEkin];
    _nNumMoleculesSumGlobal = new unsigned long[_nNumSlabsDeltaEkin];

    _dDelta2EkinTransSumLocal = new double[_nNumSlabsDeltaEkin];
    _dDelta2EkinTransSumGlobal = new double[_nNumSlabsDeltaEkin];

    for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
    {
        _nNumMoleculesSumLocal[s]  = 0;
        _nNumMoleculesSumGlobal[s]  = 0;

        _dDelta2EkinTransSumLocal[s]  = 0.;
        _dDelta2EkinTransSumGlobal[s]  = 0.;
    }

    // write out files (headers only)
    this->WriteHeaderDeltaEkin(_parent->GetDomainDecomposition(), _parent->GetDomain() );

    // <-- CALC_HEAT_SUPPLY
}

void ControlRegionT::CalcGlobalValues(DomainDecompBase* domainDecomp)
{
#ifdef ENABLE_MPI

    // ToDo: communicate _nNumSlabs is enough?? or have to communicate _nNumSlabsReserve (whole data structure)
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

//    cout << "_nNumMoleculesGlobal[0] = " << _nNumMoleculesGlobal[0] << endl;
//    cout << "_dBetaTransGlobal[0] = " << _dBetaTransGlobal[0] << endl;
//    cout << "_dTargetTemperature = " << _dTargetTemperature << endl;
//    cout << "_d2EkinTransGlobal[0] = " << _d2EkinTransGlobal[0] << endl;
//
//    cout << "_nRotDOFGlobal[0] = " << _nRotDOFGlobal[0] << endl;
//    cout << "_dBetaRotGlobal[0] = " << _dBetaRotGlobal[0] << endl;
//    cout << "_d2EkinRotGlobal[0] = " << _d2EkinRotGlobal[0] << endl;

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

    // calc position index for DeltaEkin profile
    double dPosY = mol->r(1);
    nPosIndex = (unsigned int) floor(dPosY / _dSlabWidthDeltaEkin);

    // add kinetic energy
    _accumulator->ScaleVelocityComponents(mol, vcorr, _dDelta2EkinTransSumLocal[nPosIndex] );
    mol->scale_D(Dcorr);

    // ToDo: add rotational contribution of added kinetic energy

    // count number of molecules, where energy was added
    _nNumMoleculesSumLocal[nPosIndex]++;
}

void ControlRegionT::ResetLocalValues()
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


void ControlRegionT::CalcGlobalValuesDeltaEkin()
{

#ifdef ENABLE_MPI

    MPI_Allreduce( _nNumMoleculesSumLocal, _nNumMoleculesSumGlobal, _nNumSlabsDeltaEkin, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _dDelta2EkinTransSumLocal, _dDelta2EkinTransSumGlobal, _nNumSlabsDeltaEkin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
    {
        _dDelta2EkinTransSumGlobal[s] = _dDelta2EkinTransSumLocal[s];
    }
#endif

}


void ControlRegionT::ResetValuesDeltaEkin()
{
    // reset local values
    for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
    {
        _nNumMoleculesSumLocal[s] = 0;
        _dDelta2EkinTransSumLocal[s] = 0.;
    }

}

void ControlRegionT::WriteHeaderDeltaEkin(DomainDecompBase* domainDecomp, Domain* domain)
{
    // write header
    stringstream outputstream;

    std::stringstream filenamestreamDeltaEkin;
    filenamestreamDeltaEkin << "TemperatureControl_dEkin_region" << this->GetID() << ".dat";

    string strFilename = filenamestreamDeltaEkin.str();

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif

    outputstream << "         simstep";

    for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
    {
        outputstream << "        slab[" << std::setfill('0') << std::setw(4) << s << "]";  // TODO: add width for slab number
    }
    outputstream << endl;

    ofstream fileout(strFilename.c_str(), ios::out);
    fileout << outputstream.str();
    fileout.close();

#ifdef ENABLE_MPI
    }
#endif
}

void ControlRegionT::WriteDataDeltaEkin(DomainDecompBase* domainDecomp, unsigned long simstep)
{
    // calc global values
    this->CalcGlobalValuesDeltaEkin();

    // reset local values
    this->ResetValuesDeltaEkin();

    // writing .rpf-files
    std::stringstream outputstreamDeltaEkin;
    std::stringstream filenamestreamDeltaEkin;

    filenamestreamDeltaEkin << "TemperatureControl_dEkin_region" << this->GetID() << ".dat";
    string strFilename = filenamestreamDeltaEkin.str();


    #ifdef ENABLE_MPI
        int rank = domainDecomp->getRank();
        // int numprocs = domainDecomp->getNumProcs();
        if (rank== 0)
        {
    #endif

            // simstep
            outputstreamDeltaEkin << std::setw(16) << simstep;

            // DeltaEkin data
            for(unsigned int s = 0; s<_nNumSlabsDeltaEkin; ++s)
            {
                outputstreamDeltaEkin << std::setw(16) << fixed << std::setprecision(3) << _dDelta2EkinTransSumGlobal[s];
            }
            outputstreamDeltaEkin << endl;

            ofstream fileout(strFilename.c_str(), ios::app);
            fileout << outputstreamDeltaEkin.str();
            fileout.close();

    #ifdef ENABLE_MPI
        }
    #endif
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
    _vecControlRegions.push_back(ControlRegionT(this, dLowerCorner, dUpperCorner, nNumSlabs, nComp, dTargetTemperature, dTemperatureExponent, strTransDirections, GetNumRegions()+1, _nNumSlabsDeltaEkin ) );
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


        // write out heat supply
        this->WriteDataDeltaEkin(domainDecomposition, simstep);
    }
}

void TemperatureControl::WriteDataDeltaEkin(DomainDecompBase* domainDecomp, unsigned long simstep)
{
    if(simstep % _nWriteFreqDeltaEkin != 0)
        return;

    // reset local values
    std::vector<ControlRegionT>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it).WriteDataDeltaEkin(domainDecomp, simstep);
    }
}



