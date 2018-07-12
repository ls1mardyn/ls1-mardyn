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
#include <iomanip>
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

std::ostream& operator<<( std::ostream& os, const double vec3d[3] ) {
	os << "[" << vec3d[0] << " " << vec3d[1] << " " << vec3d[2] << "]";
	return os;
}

// init static ID --> instance counting
unsigned short dec::ControlRegion::_nStaticID = 0;

// class dec::ControlRegion

dec::ControlRegion::ControlRegion(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3] )
	:
	CuboidRegionObs(parent, dLowerCorner, dUpperCorner),
	_bVacuum(false),
	_director(nullptr)
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
	}

	// particle manipulation director
	_director = new ParticleManipDirector(this);
}

dec::ControlRegion::~ControlRegion()
{
}

void dec::ControlRegion::readXML(XMLfileUnits& xmlconfig)
{
	// MettDeamon connection
	{
		bool bRet = true;
		int nVal = 0;
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
	}

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
		else
		{
			// reserve places for relevant components in map for counting particle manipulations
			_manip_count_map.reserve(numNodes);
			_manip_count_map[0].deleted.local = _manip_count_map[0].deleted.global = 0;
			_manip_count_map[0].inserted.local = _manip_count_map[0].inserted.global = 0;
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
			else {
				_compVars.at(cid).compID = cid;
				_manip_count_map[cid].deleted.local = _manip_count_map[cid].deleted.global = 0;
				_manip_count_map[cid].inserted.local = _manip_count_map[cid].inserted.global = 0;
			}
			xmlconfig.getNodeValue("density", density);
			_compVars.at(cid).density.target.local = _compVars.at(cid).density.target.global = density;
			xmlconfig.getNodeValue("proxyID", proxyID);
			if(-1 == proxyID)
				_compVars.at(cid).proxyID = proxyID;

			/*
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
			*/
		}
		xmlconfig.changecurrentnode(oldpath);
		xmlconfig.changecurrentnode("..");

		for(auto&& ti : _manip_count_map)
		{
			ti.second.changed_to_from.reserve(numNodes);
			for(auto&& ti2 : _manip_count_map)
				ti2.second.changed_to_from[ti.first].local = ti2.second.changed_to_from[ti.first].global = 0;
		}
	}
	// Establish vacuum?
	if(_compVars.at(0).density.target.global <= 1e-15)
	{
		_bVacuum = true;
		_director->setVacuum(true);
	}

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

	// init director
	_director->readXML(xmlconfig);
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

#ifndef NDEBUG
	cout << "Rank[" << domainDecomp->getRank() << "]: _dVolume.local  = " << _dVolume.local << endl;
	cout << "Rank[" << domainDecomp->getRank() << "]: _dVolume.global = " << _dVolume.global << endl;
#endif

	for(auto&& comp : _compVars)
	{
		if(_dVolume.local < 0.)
			comp.numMolecules.target.local = 0;
		else
			comp.numMolecules.target.local  = (int64_t)(comp.density.target.local  * _dVolume.local);
//		comp.numMolecules.target.global = (int64_t)(comp.density.target.global * _dVolume.global);
#ifdef ENABLE_MPI
		MPI_Allreduce( &comp.numMolecules.target.local, &comp.numMolecules.target.global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
		comp.numMolecules.target.global = comp.numMolecules.target.local;
#endif
	}

	// write file headers
	this->WriteHeaderParticleManipCount();
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

void dec::ControlRegion::MeasureDensity(Simulation* simulation, Molecule* mol)
{
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

void dec::ControlRegion::CalcGlobalValues(Simulation* simulation)
{

#ifdef ENABLE_MPI
//    if (_newcomm == MPI_COMM_NULL)
//        return;

	for(auto&& comp : _compVars)
	{
		MPI_Allreduce( &comp.numMolecules.actual.local, &comp.numMolecules.actual.global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); // _newcomm);
	}

	for(auto&& comp : _manip_count_map)
	{
		MPI_Allreduce( &comp.second.inserted.local, &comp.second.inserted.global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); // _newcomm);
		MPI_Allreduce( &comp.second.deleted.local,  &comp.second.deleted.global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); // _newcomm);

		for(auto&& comp_from : comp.second.changed_to_from)
		{
			MPI_Allreduce( &comp_from.second.local,  &comp_from.second.global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); // _newcomm);
		}
	}

#else
	for(auto&& comp : _compVars)
		comp.density.actual.local = comp.density.actual.global = comp.numMolecules.actual.local * _dInvertVolume.local;

	for(auto&& comp : _manip_count_map)
	{
		comp.second.inserted.global = comp.second.inserted.local;
		comp.second.deleted.global = comp.second.deleted.local;

		for(auto&& comp_from : comp.second.changed_to_from)
		{
			comp_from.second.global = comp_from.second.local;
		}
	}
#endif

	// calc actual density, and spread between target and actual values
	for(auto&& comp : _compVars)
	{
		comp.density.actual.global = comp.numMolecules.actual.global * _dInvertVolume.global;
		comp.density.actual.local  = comp.numMolecules.actual.local  * _dInvertVolume.local;

		comp.numMolecules.spread.local = comp.numMolecules.actual.local - comp.numMolecules.target.local;
		comp.numMolecules.spread.global = comp.numMolecules.actual.global - comp.numMolecules.target.global;
		comp.density.spread.local = comp.density.actual.local - comp.density.target.local;
		comp.density.spread.global = comp.density.actual.global - comp.density.target.global;
	}

	// change or add molecules??
	this->CheckState();

	//connection to MettDeamon
	if(true == _bMettDeamonConnected)
	{
//		_mettDeamon->StoreDensity(_dDensityGlobal);
//		_mettDeamon->StoreValuesCV(_dTargetDensity, this->GetVolume() );  // TODO: move this, so its called only once
	}

	/* DEBUG -->
	for(auto&& comp : _compVars)
	{
		cout << "actual: " << setw(7) << comp.numMolecules.actual.global << ", ";
		cout << "target: " << setw(7) << comp.numMolecules.target.global << ", ";
		cout << "spread: " << setw(7) << comp.numMolecules.spread.global << endl;
	}
	// <-- DEBUG */

	// inform director
	_director->globalValuesCalculated(simulation);
}

void dec::ControlRegion::ManipulateParticles(Simulation* simulation, Molecule* mol, bool& bDeleteMolecule)
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

	_director->ManipulateParticles(simulation, mol, bDeleteMolecule);
	return;

//    int nRank = domainDecomp->getRank();
//    double dPosY = mol->r(1);

	// check if molecule is inside
	for(uint8_t d=0; d<3; ++d)
		if( !(PositionIsInside(d, mol->r(d) ) ) ) return;

	bDeleteMolecule = false;
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
//			for(auto did:comp.deletionList)
//			{
//				if(did == mid)
//				{
//					bDeleteMolecule = true;
//				}
//			}
//			if(true == bDeleteMolecule)
//				break;
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

void dec::ControlRegion::ManipulateParticleForces(Simulation* simulation, Molecule* mol)
{
	_director->ManipulateParticleForces(simulation, mol);
}

void dec::ControlRegion::FinalizeParticleManipulation(Simulation* simulation, MainLoopAction* action)
{
	_director->FinalizeParticleManipulation(simulation, action);
}

void dec::ControlRegion::ResetLocalValues(Simulation* simulation)
{
	for(auto&& comp : _compVars)
	{
		comp.numMolecules.actual.local = 0;
		comp.particleIDs.clear();
	}

	// inform director
	_director->localValuesReseted(simulation);
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
	for(uint16_t d=0; d<3; ++d)
	{
//		cout << "bb["<<(uint32_t)d<<"]: " << bbMin[d] << " - " << bbMax[d] << endl;
		lc[d] = ( (bbMin[d] > this->GetLowerCorner(d) ) ? bbMin[d] : this->GetLowerCorner(d) );
		uc[d] = ( (bbMax[d] < this->GetUpperCorner(d) ) ? bbMax[d] : this->GetUpperCorner(d) );
		_dVolume.local *= (uc[d] - lc[d]);
	}
	_dInvertVolume.local = 1./_dVolume.local;

	// DEBUG
//	int ownRank = domainDecomp->getRank();
//	cout << "["<<ownRank<<"]: " << "bbMin, bbMax: " << bbMin << ", " << bbMax
//			<< "; lc, uc: " << lc << ", " << uc << endl;
	// DEBUG
}

void dec::ControlRegion::WriteHeaderDeletedMolecules()
{
	// domain decomposition
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

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
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

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

void dec::ControlRegion::WriteHeaderParticleManipCount()
{
	// domain decomposition
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

#ifdef ENABLE_MPI
int rank = domainDecomp.getRank();
// int numprocs = domainDecomp->getNumProcs();
if (rank != 0)
	return;
#endif

	for(auto&& comp : _manip_count_map)
	{
		// write header
		stringstream outputstream;
		stringstream sstrFilename;
		sstrFilename << "manip-count_reg" << this->GetID() << "_cid" << comp.first << ".dat";

		outputstream << "             simstep";
		outputstream << "             deleted";
		outputstream << "            inserted";
		for(auto&& comp_from : comp.second.changed_to_from)
			outputstream << "    changed_to_from" << comp_from.first;
		outputstream << endl;

		ofstream fileout(sstrFilename.str().c_str(), ios::out);
		fileout << outputstream.str();
		fileout.close();
	}
}

void dec::ControlRegion::WriteDataParticleManipCount(unsigned long simstep)
{
	if(simstep % 100 != 0)
		return;

	// domain decomposition
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

	// write out data
	#ifdef ENABLE_MPI
	int rank = domainDecomp.getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank != 0)
		return;
	#endif

	for(auto&& comp : _manip_count_map)
	{
		stringstream outputstream;
		stringstream sstrFilename;
		sstrFilename << "manip-count_reg" << this->GetID() << "_cid" << comp.first << ".dat";

		outputstream << std::setw(20) << simstep;
		outputstream << std::setw(20) << comp.second.deleted.global;
		outputstream << std::setw(20) << comp.second.inserted.global;
		for(auto&& comp_from : comp.second.changed_to_from)
			outputstream << std::setw(20) << comp_from.second.global;

		outputstream << endl;

		ofstream fileout(sstrFilename.str().c_str(), ios::app);
		fileout << outputstream.str();
		fileout.close();
	}

	// reset local values
	// count particle manipulations
	for(auto&& comp : _manip_count_map)
	{
		comp.second.inserted.local = 0;
		comp.second.deleted.local = 0;

		for(auto&& comp_from : comp.second.changed_to_from)
		{
			comp_from.second.local = 0;
		}
	}
}

// getters
uint32_t dec::ControlRegion::getGlobalMinNumMoleculesSpreadCompID()
{
	uint32_t cidMinSpread = 0;
	int64_t spreadMin = 0;
	for(auto comp:_compVars)
	{
		if(0 == comp.compID)
			continue;
		int64_t spread = comp.numMolecules.spread.global;
		if(spread < spreadMin)
		{
			spreadMin = spread;
			cidMinSpread = comp.compID;
		}
	}
	return cidMinSpread;
}

uint32_t dec::ControlRegion::getGlobalMaxNumMoleculesSpreadCompID()
{
	uint32_t cidMaxSpread = 0;
	int64_t spreadMax = 0;
	for(auto comp:_compVars)
	{
		if(0 == comp.compID)
			continue;
		int64_t spread = comp.numMolecules.spread.global;
		if(spread > spreadMax)
		{
			spreadMax = spread;
			cidMaxSpread = comp.compID;
		}
	}
	return cidMaxSpread;
}

void dec::ControlRegion::getGlobalMinMaxNumMoleculesSpreadCompIDs(uint32_t& cidMinSpread, uint32_t& cidMaxSpread)
{
	cidMinSpread = cidMaxSpread = 0;
	int64_t spreadMin, spreadMax;
	spreadMin = spreadMax = 0;
	for(auto comp:_compVars)
	{
		if(0 == comp.compID)
			continue;
		int64_t spread = comp.numMolecules.spread.global;
		if(spread < spreadMin)
		{
			spreadMin = spread;
			cidMinSpread = comp.compID;
		}
		else if(spread > spreadMax)
		{
			spreadMax = spread;
			cidMaxSpread = comp.compID;
		}
	}
//	// DEBUG
//	if(0 == cidMinSpread || 0 == cidMaxSpread)
//		cout << "FAILED" << endl;
//	// DEBUG
}

// checks
bool dec::ControlRegion::globalCompositionBalanced()
{
	bool bBalanced = true;
	uint16_t numSpreadsNotZero = 0;
	for(auto comp:_compVars)
	{
		if(0 == comp.compID)
			continue;
		if(comp.numMolecules.spread.global != 0)
			numSpreadsNotZero++;
	}
	bBalanced = numSpreadsNotZero < 2;
	return bBalanced;
}

// inform about particle manipulations
void dec::ControlRegion::informParticleInserted(Molecule mol)
{
	uint32_t cid_zb = mol.componentid();
	cid_manip_count_map_it it;
	it = _manip_count_map.find(cid_zb+1);
	it->second.inserted.local++;
	it = _manip_count_map.find(0);
	it->second.inserted.local++;

	//connection to MettDeamon
	if(true == _bMettDeamonConnected)
	{
		if(cid_zb+1 == _mettDeamon->getFeedRateTargetComponentID() )
			_mettDeamon->IncrementInsertedMoleculesLocal();
	}
}

void dec::ControlRegion::informParticleDeleted(Molecule mol)
{
	uint32_t cid_zb = mol.componentid();
	cid_manip_count_map_it it;
	it = _manip_count_map.find(cid_zb+1);
	it->second.deleted.local++;
	it = _manip_count_map.find(0);
	it->second.deleted.local++;

	//connection to MettDeamon
	if(true == _bMettDeamonConnected)
	{
		if(cid_zb+1 == _mettDeamon->getFeedRateTargetComponentID() )
			_mettDeamon->IncrementDeletedMoleculesLocal();
	}
}

void dec::ControlRegion::informParticleChanged(Molecule from, Molecule to)
{
	uint32_t cid_from_zb = from.componentid();
	uint32_t cid_to_zb = to.componentid();
	cid_manip_count_map_it it;

	// changed to this component
	it = _manip_count_map.find(cid_to_zb+1);
	it->second.changed_to_from[0].local++;
	it->second.changed_to_from[cid_from_zb+1].local++;

	// changed to arbirrary component
	it = _manip_count_map.find(0);
	it->second.changed_to_from[0].local++;
	it->second.changed_to_from[cid_from_zb+1].local++;

	//connection to MettDeamon
	if(true == _bMettDeamonConnected)
	{
		if(cid_to_zb+1 == _mettDeamon->getFeedRateTargetComponentID() )
			_mettDeamon->IncrementChangedToMoleculesLocal();
		if(cid_from_zb+1 == _mettDeamon->getFeedRateTargetComponentID() )
			_mettDeamon->IncrementChangedFromMoleculesLocal();
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

void DensityControl::MeasureDensity(Simulation* simulation, Molecule* mol)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->MeasureDensity(simulation, mol);
}

void DensityControl::CalcGlobalValues(Simulation* simulation)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->CalcGlobalValues(simulation);
}

void DensityControl::ManipulateParticles(Simulation* simulation, Molecule* mol, bool& bDeleteMolecule)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	bDeleteMolecule = false;
	for(auto&& reg : _vecControlRegions) {
		bool bDeleteMoleculeRegion = false;
		reg->ManipulateParticles(simulation, mol, bDeleteMoleculeRegion);
		bDeleteMolecule = bDeleteMolecule || bDeleteMoleculeRegion;
	}
}

void DensityControl::ManipulateParticleForces(Simulation* simulation, Molecule* mol)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->ManipulateParticleForces(simulation, mol);
}

void DensityControl::FinalizeParticleManipulation(Simulation* simulation, MainLoopAction* action)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->FinalizeParticleManipulation(simulation, action);
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

void DensityControl::ResetLocalValues(Simulation* simulation)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->ResetLocalValues(simulation);
}

void DensityControl::InitMPI()
{
	for(auto&& reg : _vecControlRegions)
		reg->InitMPI();
}

void DensityControl::WriteDataDeletedMolecules(Simulation* simulation)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->WriteDataDeletedMolecules(simulation->getSimulationStep() );
}

void DensityControl::WriteDataParticleManipCount(Simulation* simulation)
{
	if(simulation->getSimulationStep() % _nControlFreq != 0)
		return;

	for(auto&& reg : _vecControlRegions)
		reg->WriteDataParticleManipCount(simulation->getSimulationStep() );
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
	this->preFirstLoop(simulation);

	for (ParticleIterator pit = particleCont->iterator(); pit.hasNext(); pit.next())
	{
		this->insideFirstLoop(simulation, &(*pit) );
	}

	this->postFirstPreSecondLoop(simulation);

	for (ParticleIterator pit = particleCont->iterator(); pit.hasNext(); pit.next())
	{
		bool bDeleteMolecule = false;
		this->insideSecondLoop(simulation, &(*pit), bDeleteMolecule);
		if(true == bDeleteMolecule)
		{
//			cout << "DEL molecule with mid=" << pit->id() << endl;
			pit.deleteCurrentParticle();
		}
	}

	this->postSecondLoop(simulation);
}

// class PreForceAction : public MainLoopAction
void PreForceAction::preFirstLoop(Simulation* simulation)
{
	_parent->ResetLocalValues(simulation);
}
void PreForceAction::insideFirstLoop(Simulation* simulation, Molecule* mol)
{
	_parent->MeasureDensity(simulation, mol);
}
void PreForceAction::postFirstPreSecondLoop(Simulation* simulation)
{
	_parent->CalcGlobalValues(simulation);
}
void PreForceAction::insideSecondLoop(Simulation* simulation, Molecule* mol, bool& bDeleteMolecule)
{
	_parent->ManipulateParticles(simulation, mol, bDeleteMolecule);
}
void PreForceAction::postSecondLoop(Simulation* simulation)
{
	_parent->FinalizeParticleManipulation(simulation, this);
	_parent->WriteDataDeletedMolecules(simulation);
}

// class PostForceAction : public MainLoopAction
void PostForceAction::preFirstLoop(Simulation* simulation)
{
}
void PostForceAction::insideFirstLoop(Simulation* simulation, Molecule* mol)
{
}
void PostForceAction::postFirstPreSecondLoop(Simulation* simulation)
{
}
void PostForceAction::insideSecondLoop(Simulation* simulation, Molecule* mol, bool& bDeleteMolecule)
{
	_parent->ManipulateParticleForces(simulation, mol);
}
void PostForceAction::postSecondLoop(Simulation* simulation)
{
	_parent->FinalizeParticleManipulation(simulation, this);
	_parent->WriteDataParticleManipCount(simulation);
}
