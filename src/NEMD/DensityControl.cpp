/*
 * DensityControl.cpp
 *
 *  Created on: 29.05.2015
 *      Author: mheinen
 */

#include "NEMD/DensityControl.h"
#include "NEMD/DistControl.h"
#include "NEMD/MettDeamon.h"
#include "NEMD/NEMD.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "molecules/Quaternion.h"
#include "Simulation.h"
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

#include <cstdlib>
#include <cstdint>
#include <limits>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

using namespace std;

// init static ID --> instance counting
unsigned short dec::ControlRegion::_nStaticID = 0;

// class dec::ControlRegion

dec::ControlRegion::ControlRegion(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3] )
		: CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
{
	// ID
	_nID = ++_nStaticID;

	// init process relevance
	_bProcessIsRelevant = true;

	// init rank array
	_ranks = NULL;

	// reset local values
	_nDeletedNumMoleculesLocal = 0;

	for(unsigned int d=0; d<3; ++d)
	{
		_dDeletedVelocityLocal[d] = 0.;
		_dDeletedVelocitySquaredLocal[d] = 0.;
		_dDeletedEkinLocal[d] = 0.;
	}

	this->WriteHeaderDeletedMolecules();

	// Connection to MettDeamon
	_mettDeamon = NULL;
	_bMettDeamonConnected = false;
}

dec::ControlRegion::ControlRegion(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3],
		unsigned int nTargetComponentID, const double& dTargetDensity)
		: CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
{
	// ID
	_nID = ++_nStaticID;

    // target density
    _dTargetDensity = dTargetDensity;

    // target component ID (inert gas)
    _nTargetComponentID = nTargetComponentID;

    // init process relevance
    _bProcessIsRelevant = true;

    // init rank array
    _ranks = NULL;

    // reset local values
    _nDeletedNumMoleculesLocal = 0;

    for(unsigned int d=0; d<3; ++d)
    {
        _dDeletedVelocityLocal[d] = 0.;
        _dDeletedVelocitySquaredLocal[d] = 0.;
        _dDeletedEkinLocal[d] = 0.;
    }

    this->WriteHeaderDeletedMolecules();

	// Connection to MettDeamon
	_mettDeamon = NULL;
	_bMettDeamonConnected = false;
}


dec::ControlRegion::~ControlRegion()
{
}

void dec::ControlRegion::readXML(XMLfileUnits& xmlconfig)
{
	bool bRet = true;
	bool bVal = false;
	std::string strName;
	std::string oldpath = xmlconfig.getcurrentnodepath();
	bRet = bRet && xmlconfig.changecurrentnode("option/[@name='connect-mettdeamon']");
	if(true == bRet)
		bRet = bRet && xmlconfig.getNodeValue(".", bVal);
	if(true == bRet && true == bVal)
		_bMettDeamonConnected = true;

	xmlconfig.getNodeValue("target/componentID", _nTargetComponentID);
	xmlconfig.getNodeValue("target/density", _dTargetDensity);
}

void dec::ControlRegion::CheckBounds()
{
	// check if control region is outside of simulation volume
	double dBoxLengthY = _parent->GetDomain()->getGlobalLength(1);

	if ( (_dLowerCorner[1] > dBoxLengthY && _dUpperCorner[1] > dBoxLengthY) || (_dLowerCorner[1] < 0. && _dUpperCorner[1] < 0.) )
	{
		cout << "dec::ControlRegion::ControlRegion: Control region dimensions are outside of simulation volume! Programm exit..." << endl;
		exit(-1);
	}
}

void dec::ControlRegion::Init()
{
	// update region volume
	this->UpdateVolume();
}

void dec::ControlRegion::InitMPI()
{
	// domain decomposition
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
    if(_ranks != NULL)
    	delete _ranks;

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

    MPI_Group group, newgroup;

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

void dec::ControlRegion::CalcGlobalValues()
{
	// In vacuum evaporation simulation global density of control region has not to be calculated
	if( 0. == _dTargetDensity)
		return;

#ifdef ENABLE_MPI
//    if (_newcomm == MPI_COMM_NULL)
//        return;

    MPI_Allreduce( &_nNumMoleculesLocal, &_nNumMoleculesGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); // _newcomm);

    _dDensityGlobal = _nNumMoleculesGlobal * _dInvertVolume;
#else
    _dDensityGlobal = _nNumMoleculesLocal * _dInvertVolume;
#endif

}

void dec::ControlRegion::UpdateGlobalDensity(bool bDeleteMolecule)
{
	// domain decomposition
	DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();

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

void dec::ControlRegion::MeasureDensity(Molecule* mol)
{
	// In vacuum evaporation simulation global density of control region has not to be calculated
	if( 0. == _dTargetDensity)
		return;

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


void dec::ControlRegion::ControlDensity(Molecule* mol, Simulation* simulation, bool& bDeleteMolecule)
{
	// check componentID
	if(mol->componentid()+1 != _nTargetComponentID && 0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;

//    int nRank = domainDecomp->getRank();
//    double dPosY = mol->r(1);

	// check if molecule is inside
	for(uint8_t d=0; d<3; ++d)
		if( !(PositionIsInside(d, mol->r(d) ) ) ) return;

	uint32_t flagsNEMD = dynamic_cast<DensityControl*>(this->GetParent() )->GetFlagsNEMD();
	// change identity (component) of molecules Ac --> N2 (aceton to nitrogen) --> NEMD_CHANGE_COMPONENT_AC_TO_N2
	if(flagsNEMD & NEMD_CHANGE_COMPONENT_AC_TO_N2)
	{
		std::array<uint8_t, 3> arrChangComps;
		arrChangComps = {2, 1, 2};

		uint8_t cid = (uint8_t)mol->componentid();

		std::vector<Component>* ptrComps = simulation->getEnsemble()->getComponents();
		if(arrChangComps.at(cid) != cid)
		{
			Component* compOld = mol->component();
			Component* compNew = &(ptrComps->at(arrChangComps.at(cid) ) );

			// rotation
			double U_rot = mol->U_rot();
	#ifndef NDEBUG
			cout << "U_rot = " << U_rot << endl;
	#endif
			double L[3];
			double Ipa[3];
			double U_rot_FG[3];

			Ipa[0] = compNew->I11();
			Ipa[1] = compNew->I22();
			Ipa[2] = compNew->I33();

			U_rot_FG[0] = U_rot * 1./3.;  // 0.5 * 2./3.
			U_rot_FG[1] = U_rot * 1./3.;  // 0.5 * 2./3.
			U_rot_FG[2] = U_rot * 0.0;

			L[0] = sqrt(U_rot_FG[0] * 2. * Ipa[0] );
			L[1] = sqrt(U_rot_FG[1] * 2. * Ipa[1] );
			L[2] = sqrt(U_rot_FG[2] * 2. * Ipa[2] );
	#ifndef NDEBUG
			cout << "L[0] = " << L[0] << endl;
			cout << "L[1] = " << L[1] << endl;
			cout << "L[2] = " << L[2] << endl;
	#endif
			mol->setD(0, L[0] );
			mol->setD(1, L[1] );
			mol->setD(2, L[2] );

			Quaternion q(1., 0., 0., 0.);
			mol->setq(q);

			mol->setComponent(compNew);
	//		mol->clearFM();  // <-- necessary?
	#ifndef NDEBUG
			cout << "Changed cid of molecule " << mol->id() << " from: " << (int32_t)cid << " to: " << mol->componentid() << endl;
	#endif
			double mr = compOld->m()/compNew->m();
			double mrf = sqrt(mr);
			mol->scale_v(mrf);

			U_rot = mol->U_rot();
	#ifndef NDEBUG
			cout << "U_rot = " << U_rot << endl;
	#endif
		}
	}
	// <-- NEMD_CHANGE_COMPONENT_AC_TO_N2

    if( 0. == _dTargetDensity)
    {
        bDeleteMolecule = true;
    }
    else
    {
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

	// sample deleted molecules data
	if(true == bDeleteMolecule)
	{
		_nDeletedNumMoleculesLocal++;
		_dDeletedEkinLocal[0] += mol->U_kin();
		_dDeletedEkinLocal[1] += mol->U_trans();
		_dDeletedEkinLocal[2] += mol->U_rot();

		for(unsigned short d = 0; d<3; ++d)
		{
			double v = mol->v(d);
			_dDeletedVelocityLocal[d] += v;
			_dDeletedVelocitySquaredLocal[d] += v*v;
		}

		//connection to MettDeamon
		if(true == _bMettDeamonConnected)
			_mettDeamon->IncrementDeletedMoleculesLocal();
	}
}

void dec::ControlRegion::ResetLocalValues()
{
    _nNumMoleculesLocal = 0;
}

void dec::ControlRegion::WriteHeaderDeletedMolecules()
{
	// domain decomposition
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();

#ifdef ENABLE_MPI
int rank = domainDecomp.getRank();
// int numprocs = domainDecomp->getNumProcs();
if (rank != 0)
    return;
#endif

    // write header
    stringstream outputstream;
    stringstream sstrFilename;    
    sstrFilename << "DensityControl_del-mol-data_region" << this->GetID() << ".dat";

    outputstream << "         simstep";
    outputstream << "         numMols";
    outputstream << "                   U_kin";
    outputstream << "                 U_trans";
    outputstream << "                   U_rot";
    outputstream << "                      vx";
    outputstream << "                      vy";
    outputstream << "                      vz";
    outputstream << "                     vx2";
    outputstream << "                     vy2";
    outputstream << "                     vz2";
    outputstream << endl;

    ofstream fileout(sstrFilename.str().c_str(), ios::out);
    fileout << outputstream.str();
    fileout.close();
}

void dec::ControlRegion::WriteDataDeletedMolecules(unsigned long simstep)
{
	// domain decomposition
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();

    // calc global values 
#ifdef ENABLE_MPI

    MPI_Reduce( &_nDeletedNumMoleculesLocal,    &_nDeletedNumMoleculesGlobal,    1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  _dDeletedEkinLocal,             _dDeletedEkinGlobal,            3, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  _dDeletedVelocityLocal,         _dDeletedVelocityGlobal,        3, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  _dDeletedVelocitySquaredLocal,  _dDeletedVelocitySquaredGlobal, 3, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);

#else
    _nDeletedNumMoleculesGlobal = _nDeletedNumMoleculesLocal;

    for(unsigned int d=0; d<3; ++d)
    {
    	_dDeletedEkinGlobal[d] = _dDeletedEkinLocal[d];
        _dDeletedVelocityGlobal[d] = _dDeletedVelocityLocal[d];
        _dDeletedVelocitySquaredGlobal[d] = _dDeletedVelocitySquaredLocal[d];
    }
#endif

    // reset local values
    _nDeletedNumMoleculesLocal = 0;

    for(unsigned int d=0; d<3; ++d)
    {
        _dDeletedEkinLocal[d] = 0.;
        _dDeletedVelocityLocal[d] = 0.;
        _dDeletedVelocitySquaredLocal[d] = 0.;
    }

    // write out data
    #ifdef ENABLE_MPI
    int rank = domainDecomp.getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank != 0)
        return;
    #endif

    stringstream outputstream;
    stringstream sstrFilename;
    sstrFilename << "DensityControl_del-mol-data_region" << this->GetID() << ".dat";

    outputstream << std::setw(16) << simstep;
    outputstream << std::setw(16) << _nDeletedNumMoleculesGlobal;

    for(unsigned int d=0; d<3; ++d)
        outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDeletedEkinGlobal[d];

    for(unsigned int d=0; d<3; ++d)
        outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDeletedVelocityGlobal[d];

    for(unsigned int d=0; d<3; ++d)
        outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDeletedVelocitySquaredGlobal[d];

    outputstream << endl;

    ofstream fileout(sstrFilename.str().c_str(), ios::app);
    fileout << outputstream.str();
    fileout.close();
}
        
// class DensityControl

DensityControl::DensityControl(DomainDecompBase* domainDecomp, Domain* domain)
	: ControlInstance(domain, domainDecomp),
		_nStart(0),
		_nControlFreq(1),
		_nStop(1000000000),
		_nWriteFreqDeleted(1000),
		_bProcessIsRelevant(true),
		_flagsNEMD(0)
{
}

DensityControl::DensityControl(DomainDecompBase* domainDecomp, Domain* domain,
		unsigned long nControlFreq, unsigned long nStart, unsigned long nStop)
: ControlInstance(domain, domainDecomp), _flagsNEMD(0)
{
    // control frequency
    _nControlFreq = nControlFreq;

    // start/stop timestep
    _nStart = nStart;
    _nStop  = nStop;

    // deleted molecules data
    _nWriteFreqDeleted = 1000;
}

DensityControl::~DensityControl()
{
}

void DensityControl::readXML(XMLfileUnits& xmlconfig)
{
	// flags, TODO: implement mechanism to add multiple flags by logic | operation
	std::string strFlags = "none";
	xmlconfig.getNodeValue("@flags", strFlags);
	if ("NEMD_CHANGE_COMPONENT_AC_TO_N2" == strFlags)
	{
		_flagsNEMD = (_flagsNEMD | NEMD_CHANGE_COMPONENT_AC_TO_N2);
		std::string str = "false";
		if(true == (_flagsNEMD & NEMD_CHANGE_COMPONENT_AC_TO_N2) )
			str = "true";
		global_log->info() << "DensityControl: Change component AC --> N2:"
			" " << str << endl;
	}

	// control
	xmlconfig.getNodeValue("control/start", _nStart);
	xmlconfig.getNodeValue("control/controlfreq", _nControlFreq);
	xmlconfig.getNodeValue("control/writefreq", _nWriteFreqDeleted);
	xmlconfig.getNodeValue("control/stop", _nStop);

	// add regions
	uint32_t numRegions = 0;
	uint32_t nRegID = 0;
	XMLfile::Query query = xmlconfig.query("region");
	numRegions = query.card();
	global_log->info() << "DensityControl: Number of sampling regions: " << numRegions << endl;
	if(numRegions < 1) {
		global_log->warning() << "DensityControl: No region parameters specified. Program exit ..." << endl;
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

		global_log->info() << "DensityControl->region["<<nRegID<<"]: lower corner: " << lc[0] << ", " << lc[1] << ", " << lc[2] << endl;
		global_log->info() << "DensityControl->region["<<nRegID<<"]: upper corner: " << uc[0] << ", " << uc[1] << ", " << uc[2] << endl;

		// add regions
		dec::ControlRegion* region = new dec::ControlRegion(this, lc, uc);
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
				global_log->error() << "DensityControl->region["<<region->GetID()<<"]: Initialization of feature DistControl is needed before! Program exit..." << endl;
				exit(-1);
			}
		}

		// read region parameters from XML file
		region->readXML(xmlconfig);
		nRegID++;

	}  // for( outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++ )
}

void DensityControl::AddRegion(dec::ControlRegion* region)
{
    _vecControlRegions.push_back(region);

    // check / set process relevance
    std::vector<dec::ControlRegion*>::iterator it;
    bool bProcessIsRelevant = false;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        if ( true == (*it)->ProcessIsRelevant() )
            bProcessIsRelevant = true;
    }

    if(true == bProcessIsRelevant)
        _bProcessIsRelevant = true;
    else
        _bProcessIsRelevant = false;
}

void DensityControl::MeasureDensity(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // measure drift in each control region
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->MeasureDensity(mol);
    }
}

void DensityControl::CalcGlobalValues(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->CalcGlobalValues();
    }
}

void DensityControl::UpdateGlobalDensities(unsigned long simstep, bool bDeleteMolecule)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->UpdateGlobalDensity(bDeleteMolecule);
    }
}

void DensityControl::ControlDensity(Molecule* mol, Simulation* simulation, unsigned long simstep, bool& bDeleteMolecule)
{
    if(simstep % _nControlFreq != 0)
        return;

    // control drift of all regions
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ControlDensity(mol, simulation, bDeleteMolecule);
    }
}

void DensityControl::CheckRegionBounds()
{
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->CheckBounds();
    }
}

void DensityControl::Init(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // reset local values
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ResetLocalValues();

        // Init
        (*it)->Init();
    }
}

void DensityControl::InitMPI()
{
    // reset local values
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->InitMPI();
    }
}

void DensityControl::WriteDataDeletedMolecules(unsigned long simstep)
{
    if(simstep % _nWriteFreqDeleted != 0)
        return;

    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->WriteDataDeletedMolecules(simstep);
    }
}

// Connection to MettDeamon
void DensityControl::ConnectMettDeamon(MettDeamon* mettDeamon)
{
	for(auto&& ri : _vecControlRegions)
	{
		if(true == ri->MettDeamonConnected() )
			ri->ConnectMettDeamon(mettDeamon);
	}
}
