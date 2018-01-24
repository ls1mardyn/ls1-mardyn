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
#include "NEMD/ParticleInsertion.h"
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
#include <numeric>
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
	:
	CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
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

	// init identity change vector
	uint8_t nNumComponents = global_simulation->getEnsemble()->getComponents()->size()+1;
	_vecChangeCompIDs.resize(nNumComponents);
	std::iota (std::begin(_vecChangeCompIDs), std::end(_vecChangeCompIDs), 0);

	// target vector
	_compVars.resize(nNumComponents);
	for(uint8_t cid=0; cid<nNumComponents; ++cid)
	{
		_compVars.at(cid).compID = cid;
		_compVars.at(cid).proxyID = cid;
		_compVars.at(cid).numMolecules.actual.local = 0;
		_compVars.at(cid).numMolecules.actual.global = 0;
		_compVars.at(cid).numMolecules.target.local = 0;
		_compVars.at(cid).numMolecules.target.global = 0;
		_compVars.at(cid).insertion = nullptr;
	}
}

dec::ControlRegion::ControlRegion(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3],
		unsigned int nTargetComponentID, const double& dTargetDensity)
		: CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
{
	// ID
	_nID = ++_nStaticID;

    // target density
//    _dTargetDensity = dTargetDensity;

    // target component ID (inert gas)
//    _nTargetComponentID = nTargetComponentID;

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
	int nVal = 0;
	std::string strName;
	std::string oldpath = xmlconfig.getcurrentnodepath();
	bRet = bRet && xmlconfig.changecurrentnode("connect/[@name='MettDeamon']");
	if(true == bRet)
		bRet = bRet && xmlconfig.getNodeValue(".", nVal);
	if(true == bRet && nVal > 0)
	{
		_bMettDeamonConnected = true;
		_nMettDeamonInstanceIndex = (uint8_t)(nVal-1);
	}
	xmlconfig.changecurrentnode(oldpath);

	// set target variables for all components
	if(xmlconfig.changecurrentnode("targets")) {
		uint8_t numNodes = 0;
		XMLfile::Query query = xmlconfig.query("target");
		numNodes = query.card();
		global_log->info() << "Number of targets: " << (uint32_t)numNodes << endl;
		if(numNodes < 1) {
			global_log->error() << "No targets defined in XML-config file. Program exit ..." << endl;
			Simulation::exit(-1);
		}
		string oldpath = xmlconfig.getcurrentnodepath();
		XMLfile::Query::const_iterator nodeIter;
		for( nodeIter = query.begin(); nodeIter; nodeIter++ ) {
			xmlconfig.changecurrentnode(nodeIter);
			double density = -1.0;
			int32_t proxyID = -1;
			int32_t cid = -1;
			xmlconfig.getNodeValue("@cid", cid);
			if(-1 == cid) {
				global_log->error() << "Missing attribute 'cid' in target. Program exit ..." << endl;
				Simulation::exit(-1);
			}
			else
				_compVars.at(cid).compID = cid;
			xmlconfig.getNodeValue("density", density);
			_compVars.at(cid).density.target.local = _compVars.at(cid).density.target.global = density;
			xmlconfig.getNodeValue("proxyID", proxyID);
			if(-1 == proxyID)
				_compVars.at(cid).proxyID = proxyID;

			// particle insertion
			std::string strType = "unknown";
			xmlconfig.getNodeValue("insertion@type", strType);
			cout << "Insertiontype: " << strType << endl;
			if("BubbleMethod" == strType)
			{
				_compVars.at(cid).insertion = new BubbleMethod(this);
				string oldpath = xmlconfig.getcurrentnodepath();
				xmlconfig.changecurrentnode("insertion");
				_compVars.at(cid).insertion->readXML(xmlconfig);
				xmlconfig.changecurrentnode(oldpath);
			}
		}
		xmlconfig.changecurrentnode(oldpath);
		xmlconfig.changecurrentnode("..");
	}

	cout << "current path: " << xmlconfig.getcurrentnodepath() << endl;

	// change identity of molecules by component ID
	if(xmlconfig.changecurrentnode("changes")) {
		uint8_t numChanges = 0;
		XMLfile::Query query = xmlconfig.query("change");
		numChanges = query.card();
		global_log->info() << "Number of components to change: " << (uint32_t)numChanges << endl;
		if(numChanges < 1) {
			global_log->error() << "No component change defined in XML-config file. Program exit ..." << endl;
			Simulation::exit(-1);
		}
		string oldpath = xmlconfig.getcurrentnodepath();
		XMLfile::Query::const_iterator changeIter;
		for( changeIter = query.begin(); changeIter; changeIter++ ) {
			xmlconfig.changecurrentnode(changeIter);
			uint32_t nFrom, nTo;
			nFrom = nTo = 1;
			xmlconfig.getNodeValue("from", nFrom);
			xmlconfig.getNodeValue("to", nTo);
			_vecChangeCompIDs.at(nFrom-1) = nTo-1;
		}
		xmlconfig.changecurrentnode(oldpath);
		xmlconfig.changecurrentnode("..");
	}
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

void dec::ControlRegion::Init(DomainDecompBase* domainDecomp)
{
	// update region volume
	this->UpdateVolume(domainDecomp);

	// DEBUG
//	cout << "Rank[" << domainDecomp->getRank() << "]: _dVolume.local  = " << _dVolume.local << endl;
//	cout << "Rank[" << domainDecomp->getRank() << "]: _dVolume.global = " << _dVolume.global << endl;
	// DEBUG

	for(auto&& comp : _compVars)
	{
		comp.numMolecules.target.local  = comp.density.target.local  * _dVolume.local;
		comp.numMolecules.target.global = comp.density.target.global * _dVolume.global;
	}
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

void dec::ControlRegion::MeasureDensity(Molecule* mol)
{
	// In vacuum evaporation simulation global density of control region has not to be calculated
	if( 0. == _compVars.at(0).density.target.global)
		return;

    // check if molecule inside control region
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

	// count number of molecules
	uint32_t cid = mol->componentid()+1;
	uint64_t mid = mol->id();
	_compVars.at(0).numMolecules.actual.local++;
	_compVars.at(cid).numMolecules.actual.local++;

	// store particle IDs
	_compVars.at(cid).particleIDs.push_back(mid);
}

void dec::ControlRegion::CalcGlobalValues()
{
	// In vacuum evaporation simulation global density of control region has not to be calculated
	if( 0. == _compVars.at(0).density.target.global)
		return;

#ifdef ENABLE_MPI
//    if (_newcomm == MPI_COMM_NULL)
//        return;

	for(auto&& comp : _compVars)
	{
		MPI_Allreduce( &comp.numMolecules.actual.local, &comp.numMolecules.actual.global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); // _newcomm);
		comp.density.actual.global = comp.numMolecules.actual.global * _dInvertVolume.global;
		comp.density.actual.local  = comp.numMolecules.actual.local  * _dInvertVolume.local;
	}
#else
	for(auto&& comp : _compVars)
		comp.density.actual.local = comp.density.actual.global = comp.numMolecules.actual.local * _dInvertVolume.local;
#endif

	// calc spread between target and actual values
	for(auto&& comp : _compVars)
	{
		comp.numMolecules.spread.local = comp.numMolecules.actual.local - comp.numMolecules.target.local;
		comp.numMolecules.spread.global = comp.numMolecules.actual.global - comp.numMolecules.target.global;
		comp.density.spread.local = comp.density.actual.local - comp.density.target.local;
		comp.density.spread.global = comp.density.actual.global - comp.density.target.global;
	}

	// change or add molecules??
	this->CheckState();

	if(true == _bMettDeamonConnected)
	{
//		_mettDeamon->StoreDensity(_dDensityGlobal);
//		_mettDeamon->StoreValuesCV(_dTargetDensity, this->GetVolume() );  // TODO: move this, so its called only once
	}

	// particle insertion
	if(_nState == CRS_INSERT_MOLECULES)
	{
		uint32_t cid = this->NextInsertionID();
		ParticleInsertion* insertion = _compVars.at(cid).insertion;
		if(nullptr != insertion)
		{
			if(this->InsertionAllIdle() == true)
				insertion->setState(BMS_SELECT_RANDOM_MOLECULE);
			insertion->preLoopAction();
		}
	}
}

void dec::ControlRegion::CreateDeletionLists()
{
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	uint8_t numComps = _compVars.size();
	commVar<uint64_t> numDel[numComps];
	if(_compVars.at(0).numMolecules.spread.global > 0)
		numDel[0].global = _compVars.at(0).numMolecules.spread.global;
	else
	{
		for(auto&& comp:_compVars)
			comp.deletionList.clear();
		return;
	}
	commVar<int64_t> positiveSpread[numComps];
	commVar<int64_t> positiveSpreadSumOverComp;
	positiveSpreadSumOverComp.local  = 0;
	positiveSpreadSumOverComp.global = 0;

	for(uint8_t cid=1; cid<numComps; ++cid)
	{
		commVar<int64_t> spread;
		// sum over components
		spread.global = _compVars.at(cid).numMolecules.spread.global;
		if(spread.global > 0)
			positiveSpreadSumOverComp.global += spread.global;
		// sum over processes
		spread.local = _compVars.at(cid).numMolecules.spread.local;
		if(spread.local > 0)
			positiveSpread[cid].local = spread.local;
		else
			positiveSpread[cid].local = 0;

	#ifdef ENABLE_MPI
		MPI_Allreduce( &positiveSpread[cid].local, &positiveSpread[cid].global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	#else
		positiveSpread[cid].global = positiveSpread[cid].local;
	#endif
	}

	commVar<double> dInvPositiveSpread[numComps];
	commVar<double> dInvPositiveSpreadSumOverComp;
	dInvPositiveSpreadSumOverComp.global = 1./( (double) (positiveSpreadSumOverComp.global) );

	for(uint8_t cid=1; cid<numComps; ++cid)
	{
		commVar<int64_t> spread;
		// global
		spread.global = _compVars.at(cid).numMolecules.spread.global;
		if(spread.global > 0)
			numDel[cid].global = round(spread.global * dInvPositiveSpreadSumOverComp.global * numDel[0].global);
		else
			numDel[cid].global = 0;
		// local
		dInvPositiveSpread[cid].global = 1./( (double) (positiveSpread[cid].global) );
		spread.local = _compVars.at(cid).numMolecules.spread.local;
		cout << "Rank[" << domainDecomp.getRank() << "]: spread.local["<<(uint32_t)(cid)<<"] = " << spread.local << endl;
		if(spread.local > 0)
		{
			cout << "Rank[" << domainDecomp.getRank() << "]: numDel["<<(uint32_t)(cid)<<"].local = " << numDel[cid].local << endl;
			numDel[cid].local = round(spread.local * dInvPositiveSpread[cid].global * numDel[cid].global);
			this->select_rnd_elements(_compVars.at(cid).particleIDs, _compVars.at(cid).deletionList, numDel[cid].local);
		}
		else
			numDel[cid].local = 0;

		//DEBUG
		cout << "Rank[" << domainDecomp.getRank() << "]_compVars.at("<<(uint32_t)(cid)<<").deletionList:";
		for(auto mid:_compVars.at(cid).deletionList)
		{
			cout << " " << mid;
		}
		cout << endl;
		//DEBUG
	}
}

void dec::ControlRegion::ControlDensity(Molecule* mol, Simulation* simulation)
{
	ParticleContainer* particleCont = simulation->getMolecules();
//	// particle insertion
//	{
//		if(nullptr != _insertion)
//			_insertion->insideLoopAction(mol);
//	}

	// check componentID
	uint32_t cid = mol->componentid();
	/*
	if(cid+1 != _nTargetComponentID && 0 != _nTargetComponentID)  // program intern componentID starts with 0
		return;
	*/

//    int nRank = domainDecomp->getRank();
//    double dPosY = mol->r(1);

	// check if molecule is inside
	for(uint8_t d=0; d<3; ++d)
		if( !(PositionIsInside(d, mol->r(d) ) ) ) return;

	// --> CHANGE_IDENTITY
	if(cid != _vecChangeCompIDs.at(cid))
	{
		std::vector<Component>* ptrComps = simulation->getEnsemble()->getComponents();
		Component* compOld = mol->component();
		Component* compNew = &(ptrComps->at(_vecChangeCompIDs.at(cid) ) );

		uint8_t numRotDOF_old = compOld->getRotationalDegreesOfFreedom();
		uint8_t numRotDOF_new = compNew->getRotationalDegreesOfFreedom();
		double dUkinOld = mol->U_kin();
		double dUkinPerDOF = dUkinOld / (3 + numRotDOF_old);

		// rotation
		double U_rot = mol->U_rot();
		global_log->info() << "U_rot_old = " << U_rot << endl;

		double L[3];
		double Ipa[3];

		Ipa[0] = compNew->I11();
		Ipa[1] = compNew->I22();
		Ipa[2] = compNew->I33();

		for(uint8_t dim=0; dim<3; ++dim)
		{
			L[dim] = sqrt(dUkinPerDOF * 2. * Ipa[dim] );
			mol->setD(dim, L[dim] );
		}

#ifndef NDEBUG
		cout << "L[0] = " << L[0] << endl;
		cout << "L[1] = " << L[1] << endl;
		cout << "L[2] = " << L[2] << endl;
#endif
		Quaternion q(1., 0., 0., 0.);
		mol->setq(q);

		mol->setComponent(compNew);
//		mol->clearFM();  // <-- necessary?
#ifndef NDEBUG
		cout << "Changed cid of molecule " << mol->id() << " from: " << (int32_t)cid << " to: " << mol->componentid() << endl;
#endif

		// update transl. kin. energy
		double dScaleFactorTrans = sqrt(6*dUkinPerDOF/compNew->m()/mol->v2() );
		mol->scale_v(dScaleFactorTrans);

		U_rot = mol->U_rot();
		global_log->info() << "U_rot_new = " << U_rot << endl;

		//connection to MettDeamon
		if(NULL != _mettDeamon)
			_mettDeamon->IncrementChangedMoleculesLocal();
	}
	// <-- CHANGE_IDENTITY

	bool bDeleteMolecule;
    if( 0.0000001 > _compVars.at(0).density.target.global)
    {
        bDeleteMolecule = true;
    }
	else if(_compVars.at(0).density.spread.global > 0)
	{
		uint64_t mid = mol->id();
		bDeleteMolecule = false;
		for(auto comp:_compVars)
		{
			for(auto did:comp.deletionList)
			{
				if(did == mid)
				{
					bDeleteMolecule = true;
				}
			}
			if(true == bDeleteMolecule)
				break;
		}
	}

	// sample deleted molecules data
	if(true == bDeleteMolecule)
	{
		particleCont->deleteMolecule(*mol, false);
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
		if(NULL != _mettDeamon)
			_mettDeamon->IncrementDeletedMoleculesLocal();
	}
}

void dec::ControlRegion::postLoopAction()
{
	for(auto comp:_compVars)
	{
		if(comp.insertion == nullptr)
			return;
		comp.insertion->postLoopAction();
	}
}

void dec::ControlRegion::postEventNewTimestepAction()
{
	for(auto comp:_compVars)
	{
		if(comp.insertion == nullptr)
			return;
		comp.insertion->postEventNewTimestepAction();
	}
}

void dec::ControlRegion::postUpdateForcesAction()
{
	for(auto comp:_compVars)
	{
		if(comp.insertion == nullptr)
			return;
		comp.insertion->postUpdateForcesAction();
	}
}

void dec::ControlRegion::ResetLocalValues()
{
	for(auto&& comp : _compVars)
	{
		comp.numMolecules.actual.local = 0;
		comp.particleIDs.clear();
	}
}

void dec::ControlRegion::UpdateVolume(DomainDecompBase* domainDecomp)
{
	// global
	_dVolume.global = this->GetWidth(0) * this->GetWidth(1) * this->GetWidth(2);
	_dInvertVolume.global = 1./_dVolume.global;

	// local
	Domain* domain = global_simulation->getDomain();
	double bbMin[3], bbMax[3];
	double lc[3], uc[3];
	domainDecomp->getBoundingBoxMinMax(domain, bbMin, bbMax);
	_dVolume.local = 1.0;
	for(uint8_t d=0; d<3; ++d)
	{
//		cout << "bb["<<(uint32_t)d<<"]: " << bbMin[d] << " - " << bbMax[d] << endl;
		lc[d] = ( (bbMin[d] > this->GetLowerCorner(d) ) ? bbMin[d] : this->GetLowerCorner(d) );
		uc[d] = ( (bbMax[d] < this->GetUpperCorner(d) ) ? bbMax[d] : this->GetUpperCorner(d) );
		_dVolume.local *= (uc[d] - lc[d]);
	}
	_dInvertVolume.local = 1./_dVolume.local;
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

bool dec::ControlRegion::InsertionAllIdle()
{
	bool bRet = true;
	for(auto comp : _compVars)
		bRet = bRet && comp.insertion->getState() == BMS_IDLE;
	return bRet;
}

template<typename T1, typename T2>
void dec::ControlRegion::select_rnd_elements(std::list<T1>& mylist, std::vector<T1>& myvec, T2 numSelect)
{
	myvec.clear();
	if(numSelect < 1)
		return;
	uint64_t numElements = mylist.size();
	uint64_t numElementsSub = numElements / numSelect;
	T2 numResidual = numElements % numSelect;

	std::vector<std::vector<T1> > mat;
	mat.resize(numSelect);
	for(T2 ci=0; ci<numSelect; ++ci)
	{
		if(ci<numResidual)
			mat.at(ci).resize(numElementsSub+1);
		else
			mat.at(ci).resize(numElementsSub);

		for(auto&& eli : mat.at(ci))
		{
			eli = mylist.front();
			mylist.pop_front();
		}

		std::srand(std::time(nullptr));
		int rnd = rand()%numElementsSub;
		myvec.push_back(mat.at(ci).at(rnd) );
	}
}


// class DensityControl

DensityControl::DensityControl(DomainDecompBase* domainDecomp, Domain* domain)
	: ControlInstance(domain, domainDecomp),
		_nStart(0),
		_nControlFreq(1),
		_nStop(1000000000),
		_nWriteFreqDeleted(1000),
		_bProcessIsRelevant(true),
		_flagsNEMD(0),
		_preForceAction(nullptr),
		_postForceAction(nullptr)
{
	_preForceAction = new PreForceAction(this);
	_postForceAction = new PostForceAction(this);
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
    bool bProcessIsRelevant = false;

	for(auto&& reg : _vecControlRegions)
	{
		if ( true == reg->ProcessIsRelevant() )
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

	for(auto&& reg : _vecControlRegions)
		reg->MeasureDensity(mol);
}

void DensityControl::CalcGlobalValues(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

	for(auto&& reg : _vecControlRegions)
		reg->CalcGlobalValues();
}

void DensityControl::CreateDeletionLists()
{
	for(auto&& reg : _vecControlRegions)
		reg->CreateDeletionLists();
}

void DensityControl::ControlDensity(Molecule* mol, Simulation* simulation)
{
    if(simulation->getSimulationStep() % _nControlFreq != 0)
        return;

	for(auto&& reg : _vecControlRegions)
		reg->ControlDensity(mol, simulation);
}

void DensityControl::postLoopAction()
{
	for(auto&& reg : _vecControlRegions)
		reg->postLoopAction();
}

void DensityControl::postEventNewTimestepAction()
{
	for(auto&& reg : _vecControlRegions)
		reg->postEventNewTimestepAction();
}

void DensityControl::postUpdateForcesAction()
{
	for(auto&& reg : _vecControlRegions)
		reg->postUpdateForcesAction();
}

void DensityControl::CheckRegionBounds()
{
	for(auto&& reg : _vecControlRegions)
		reg->CheckBounds();
}

void DensityControl::Init(DomainDecompBase* domainDecomp)
{
	for(auto&& reg : _vecControlRegions)
		reg->Init(domainDecomp);
}

void DensityControl::ResetLocalValues(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

	for(auto&& reg : _vecControlRegions)
		reg->ResetLocalValues();
}

void DensityControl::InitMPI()
{
	for(auto&& reg : _vecControlRegions)
		reg->InitMPI();
}

void DensityControl::WriteDataDeletedMolecules(unsigned long simstep)
{
    if(simstep % _nWriteFreqDeleted != 0)
        return;

	for(auto&& reg : _vecControlRegions)
		reg->WriteDataDeletedMolecules(simstep);
}

void DensityControl::preForce_action(Simulation* simulation)
{
	_preForceAction->performAction(simulation);
}
void DensityControl::postForce_action(Simulation* simulation)
{
	_postForceAction->performAction(simulation);
}

// Connection to MettDeamon
void DensityControl::ConnectMettDeamon(const std::vector<MettDeamon*>& mettDeamon)
{
	for(auto&& reg : _vecControlRegions)
		reg->ConnectMettDeamon(mettDeamon);
}

// class MainLoopAction

void MainLoopAction::performAction(Simulation* simulation)
{
	ParticleContainer* particleCont = simulation->getMolecules();
	uint64_t simstep = simulation->getSimulationStep();
	this->preFirstLoop(simstep);

	ParticleIterator pit;
	for( pit  = particleCont->iteratorBegin();
			pit != particleCont->iteratorEnd();
		 ++pit )
	{
		this->insideFirstLoop( &(*pit), simstep);
	}

	this->postFirstPreSecondLoop(simstep);

	for( pit  = particleCont->iteratorBegin();
			pit != particleCont->iteratorEnd();
		 ++pit )
	{
		this->insideSecondLoop( &(*pit), simulation);
	}

	this->postSecondLoop(simstep);
}

// class PreForceAction : public MainLoopAction
void PreForceAction::preFirstLoop(unsigned long simstep)
{
	_parent->ResetLocalValues(simstep);
	_parent->postEventNewTimestepAction();
}
void PreForceAction::insideFirstLoop(Molecule* mol, unsigned long simstep)
{
	_parent->MeasureDensity(mol, simstep);
}
void PreForceAction::postFirstPreSecondLoop(unsigned long simstep)
{
	_parent->CalcGlobalValues(simstep);
	_parent->CreateDeletionLists();
}
void PreForceAction::insideSecondLoop(Molecule* mol, Simulation* simulation)
{
	ParticleContainer* particleCont = simulation->getMolecules();
	_parent->ControlDensity(mol, simulation);
}
void PreForceAction::postSecondLoop(unsigned long simstep)
{
	_parent->WriteDataDeletedMolecules(simstep);
}

// class PostForceAction : public MainLoopAction
void PostForceAction::preFirstLoop(unsigned long simstep)
{
}
void PostForceAction::insideFirstLoop(Molecule* mol, unsigned long simstep)
{
}
void PostForceAction::postFirstPreSecondLoop(unsigned long simstep)
{
}
void PostForceAction::insideSecondLoop(Molecule* mol, Simulation* simulation)
{
}
void PostForceAction::postSecondLoop(unsigned long simstep)
{
}
