/*
 * MettDeamon.cpp
 *
 *  Created on: 03.04.2017
 *      Author: thet
 */

#include "MettDeamon.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/xmlfileUnits.h"

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

MettDeamon::MettDeamon(double cutoffRadius)
	: 	_rho_l(0.),
		_dAreaXZ(0.),
		_dInvDensityArea(0.),
		_dY(0.),
		_dYsum(0.),
		_velocityBarrier(0.),
		_dSlabWidthInit(0),
		_dSlabWidth(0),
		_cutoffRadius(cutoffRadius),
		_nUpdateFreq(0),
		_nWriteFreqRestart(0),
		_maxId(0),
		_nNumMoleculesDeletedLocal(0),
		_nNumMoleculesDeletedGlobal(0),
		_nNumMoleculesDeletedGlobalAlltime(0),
		_nNumMoleculesTooFast(0),
		_nNumMoleculesTooFastGlobal(0),
		_reservoirNumMolecules(0),
		_reservoirSlabs(0),
		_nSlabindex(0),
		_reservoirFilename("unknown"),
		_bIsRestart(false),
		_nNumValsSummation(0),
		_numDeletedMolsSum(0),
		_dDeletedMolsPerTimestep(0.),
		_dInvNumTimestepsSummation(0.),
		_bMirrorActivated(false),
		_dMirrorPosY(0.)
{
	_dAreaXZ = global_simulation->getDomain()->getGlobalLength(0) * global_simulation->getDomain()->getGlobalLength(2);

	// init restart file
	std::ofstream ofs("MettDeamonRestart.dat", std::ios::out);
	std::stringstream outputstream;
	outputstream << "     simstep" << "   slabIndex" << "                  deltaY" << std::endl;
	ofs << outputstream.str();
	ofs.close();

	// summation of deleted molecules
	_listDeletedMolecules.clear();
	_listDeletedMolecules.push_back(0);
}

MettDeamon::~MettDeamon()
{
}

void MettDeamon::readXML(XMLfileUnits& xmlconfig)
{
	// control
	xmlconfig.getNodeValue("control/updatefreq", _nUpdateFreq);
	xmlconfig.getNodeValue("control/writefreq", _nWriteFreqRestart);
	xmlconfig.getNodeValue("control/vmax", _velocityBarrier);
	xmlconfig.getNodeValue("control/numvals", _nNumValsSummation);
	_dInvNumTimestepsSummation = 1. / (double)(_nNumValsSummation*_nUpdateFreq);

	// reservoir
	xmlconfig.getNodeValue("reservoir/file", _reservoirFilename);
	xmlconfig.getNodeValue("reservoir/slabwidth", _dSlabWidthInit);

	// restart
	_bIsRestart = true;
	_bIsRestart = _bIsRestart && xmlconfig.getNodeValue("restart/slabindex", _nSlabindex);
	_bIsRestart = _bIsRestart && xmlconfig.getNodeValue("restart/deltaY", _dYsum);

	// mirror
	bool bRet = xmlconfig.getNodeValue("mirror/position", _dMirrorPosY);
	_bMirrorActivated = bRet;
}

void MettDeamon::ReadReservoir(DomainDecompBase* domainDecomp)
{
	std::ifstream ifs;
	global_log->info() << "Opening Mettdeamon Reservoirfile " << _reservoirFilename << endl;
	ifs.open( _reservoirFilename.c_str() );
	if (!ifs.is_open()) {
		global_log->error() << "Could not open Mettdeamon Reservoirfile " << _reservoirFilename << endl;
		Simulation::exit(1);
	}
	global_log->info() << "Reading Mettdeamon Reservoirfile " << _reservoirFilename << endl;

	string token;
	vector<Component>& dcomponents = *(_simulation.getEnsemble()->getComponents());
	unsigned int numcomponents = dcomponents.size();
	unsigned long maxid = 0; // stores the highest molecule ID found in the phase space file
	string ntypestring("ICRVQD");
	enum Ndatatype { ICRVQDV, ICRVQD, IRV, ICRV } ntype = ICRVQD;

	double Xlength, Ylength, Zlength, V;
	Xlength = Ylength = Zlength = V = 1.0;
	while(ifs && (token != "NumberOfMolecules") && (token != "N"))
	{
		ifs >> token;

		if(token=="Length" || token=="L")
		{
			ifs >> Xlength >> Ylength >> Zlength;
			_reservoirSlabs = Ylength/_dSlabWidthInit;  // parameter sliceWidth
			global_log->info() << "Mettdeamon: _reservoirSlabs=" << _reservoirSlabs << endl;
			_dSlabWidth = Ylength / (double)(_reservoirSlabs);
//				_reservoirLastSlice = Ylength-((_reservoirSlices-1)*_cutoffRadius);
			_reservoir.resize(_reservoirSlabs);
			V = Xlength * Ylength * Zlength;
		}
	}

	if((token != "NumberOfMolecules") && (token != "N")) {
		global_log->error() << "Expected the token 'NumberOfMolecules (N)' instead of '" << token << "'" << endl;
		Simulation::exit(1);
	}
	ifs >> _reservoirNumMolecules;
	global_log->info() << " number of Mettdeamon Reservoirmolecules: " << _reservoirNumMolecules << endl;
	_rho_l = _reservoirNumMolecules / V;
	_dInvDensityArea = 1. / (_dAreaXZ * _rho_l);
	global_log->info() << "Density of Mettdeamon Reservoir: " << _rho_l << endl;

	streampos spos = ifs.tellg();
	ifs >> token;
	if((token=="MoleculeFormat") || (token == "M"))
	{
		ifs >> ntypestring;
		ntypestring.erase( ntypestring.find_last_not_of( " \t\n") + 1 );
		ntypestring.erase( 0, ntypestring.find_first_not_of( " \t\n" ) );

		if (ntypestring == "ICRVQDV") ntype = ICRVQDV;
		else if (ntypestring == "ICRVQD") ntype = ICRVQD;
		else if (ntypestring == "ICRV") ntype = ICRV;
		else if (ntypestring == "IRV")  ntype = IRV;
		else {
			global_log->error() << "Unknown molecule format '" << ntypestring << "'" << endl;
			Simulation::exit(1);
		}
	} else {
		ifs.seekg(spos);
	}
	global_log->info() << " molecule format: " << ntypestring << endl;

	if( numcomponents < 1 ) {
		global_log->warning() << "No components defined! Setting up single one-centered LJ" << endl;
		numcomponents = 1;
		dcomponents.resize( numcomponents );
		dcomponents[0].setID(0);
		dcomponents[0].addLJcenter(0., 0., 0., 1., 1., 1., 6., false);
	}
	double x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz, Vix, Viy, Viz;
	unsigned long id;
	unsigned int componentid;

	x=y=z=vx=vy=vz=q1=q2=q3=Dx=Dy=Dz=Vix=Viy=Viz=0.;
	q0=1.;

	for( unsigned long i = 0; i < _reservoirNumMolecules; i++ )
	{
		switch ( ntype ) {
			case ICRVQDV:
				ifs >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
					>> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz >> Vix >> Viy >> Viz;
				break;
			case ICRVQD:
				ifs >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
					>> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz;
				break;
			case ICRV :
				ifs >> id >> componentid >> x >> y >> z >> vx >> vy >> vz;
				break;
			case IRV :
				ifs >> id >> x >> y >> z >> vx >> vy >> vz;
				break;
		}

		if( componentid > numcomponents ) {
			global_log->error() << "Molecule id " << id << " has wrong componentid: " << componentid << ">" << numcomponents << endl;
			Simulation::exit(1);
		}
		componentid --; // TODO: Component IDs start with 0 in the program.
		//Simon changed id to i+1
		Molecule m1 = Molecule(i+1,&dcomponents[componentid],x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz);

//				long temp= m1.r(1)/_cutoffRadius;
//				_reservoirIter = _reservoir.begin();
		uint32_t nSlabindex = floor(y / _dSlabWidth);
		m1.setr(1, y - nSlabindex*_dSlabWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)

		double bbMin[3];
		double bbMax[3];
		bool bIsInsideSubdomain = false;
		domainDecomp->getBoundingBoxMinMax(global_simulation->getDomain(), bbMin, bbMax);
		bIsInsideSubdomain = x > bbMin[0] && x < bbMax[0] && y > bbMin[1] && y < bbMax[1] && z > bbMin[2] && z < bbMax[2];

		if(true == bIsInsideSubdomain)
			_reservoir.at(nSlabindex).push_back(m1);

		componentid=m1.componentid();
		// TODO: The following should be done by the addPartice method.
		dcomponents[componentid].incNumMolecules();

		if(id > maxid) maxid = id;
		_maxId = maxid;

		// Print status message
		unsigned long iph = _reservoirNumMolecules / 100;
		if( iph != 0 && (i % iph) == 0 )
			global_log->info() << "Finished reading molecules: " << i/iph << "%\r" << flush;
	}

	global_log->info() << "Finished reading Mettdeamon Rerservoirmolecules: 100%" << endl;

	ifs.close();
	if(false == _bIsRestart)
		_nSlabindex = _reservoirSlabs-1;
}

uint64_t MettDeamon::getnNumMoleculesDeleted2( DomainDecompBase* domainDecomposition)
{
	domainDecomposition->collCommInit(1);
		domainDecomposition->collCommAppendUnsLong(_nNumMoleculesTooFast);
		domainDecomposition->collCommAllreduceSum();
		_nNumMoleculesTooFastGlobal = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
//
//		std::cout << "Particles deleted: "<< _nNumMoleculesDeletedGlobalAlltime << std::endl;
//		std::cout << "Of which were too fast: " << _nNumMoleculesTooFastGlobal << std::endl;
		return _nNumMoleculesTooFastGlobal;
}
void MettDeamon::prepare_start(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer, double cutoffRadius)
{
	this->ReadReservoir(domainDecomp);

	//ParticleContainer* _moleculeContainer;
	particleContainer->deleteOuterParticles();
	// fixed components

	Component* comp2;
	comp2 = global_simulation->getEnsemble()->getComponent(1);

	Molecule* tM;
	double dPosY;
	double dBoxY = global_simulation->getDomain()->getGlobalLength(1);
	double dLeftMirror = 2*_dSlabWidth;
	global_log->info() << "Position of MettDeamons insertion area: " << dLeftMirror << endl;

	if(true == _bIsRestart)
		return;

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd();
			++tM)
	{
		dPosY = tM->r(1);

//		if(dPosY < 2.*_cutoffRadius)
		if(dPosY < (dLeftMirror-0.5) )
		{
//			cout << "cid(old) = " << tM->componentid() << endl;
			//tM->setXYZ();
			tM->setComponent(comp2);
//			cout << "cid(new) = " << tM->componentid() << endl;
		}
		// vapor phase
		else if(dPosY < (dLeftMirror+0.5) )  // || dPosY > (dBoxY-2.* cutoffRadius) ) <-- vacuum established by feature: DensityControl
		{
			particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), false);
//			cout << "delete: dY = " << dPosY << endl;
			particleContainer->update();
			tM  = particleContainer->iteratorBegin();
		}
		else if(dPosY < (dLeftMirror+1.0) )
		{
			tM->setv(1, 0.);
		}
	}
	particleContainer->update();
}
void MettDeamon::init_positionMap(ParticleContainer* particleContainer)
{
	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM) {

		unsigned long mid = tM->id();
		unsigned int  cid = tM->componentid()+1;

		if(cid == 2)
		{
			//savevelo
			std::array<double, 6> pos;
			pos.at(0) = tM->r(0);
			pos.at(1) = tM->r(1);
			pos.at(2) = tM->r(2);
			pos.at(3) = tM->v(0);
			pos.at(4) = tM->v(1);
			pos.at(5) = tM->v(2);
//			_storePosition.insert ( std::pair<unsigned long, std::array<double, 3> >(mid, pos) );
			_storePosition[tM->id()] = pos;
		}
	}
}
void MettDeamon::preForce_action(ParticleContainer* particleContainer, double cutoffRadius)
{
	double dBoxY = global_simulation->getDomain()->getGlobalLength(1);
	double dMirrorPosLeft = 2*_dSlabWidth;
	unsigned int cid;
	std::map<unsigned long, std::array<double, 6> >::iterator it;

	Component* comp1=global_simulation->getEnsemble()->getComponent(0);
	Component* comp2=global_simulation->getEnsemble()->getComponent(1);
	Component* comp3=global_simulation->getEnsemble()->getComponent(2);

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM) {

		cid = tM->componentid()+1;  // +1: cid starts with 0
		double dY = tM->r(1);

//		if(dY > dMirrorPosLeft && dY < (dMirrorPosLeft + 2.5) )
//			tM->setComponent(comp3);
//		else if(dY > (dMirrorPosLeft + 2.5) )
//			tM->setComponent(comp1);

		if( (cid == 2) && (dY > dMirrorPosLeft) )
		{
			tM->setComponent(comp1);
			tM->setv(1, abs(tM->v(1) ) );
		}

		// reset position of molecules of component 2 (cid == 2 --> fixed molecules)
		if(cid == 2)
		{
			it = _storePosition.find(tM->id() );
			if(it != _storePosition.end() )
			{
				tM->setr(0, it->second.at(0) );
				tM->setr(1, it->second.at(1) + this->getDeltaY() );
				tM->setr(2, it->second.at(2) );
			}
			// restore velocity
			tM->setv(0,it->second.at(3) );
			tM->setv(1,it->second.at(4) );
			tM->setv(2,it->second.at(5) );
		}

		/** Vacuum will be established by DensityControl
		 *
		 *
		// delete molecules of component 1 (cid == 1), close to bounding box on the right side
		if(cid == 1)
		{
			if(dY > (dBoxY - 2*cutoffRadius) )
			{
				particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), false);
				_nNumMoleculesDeletedLocal++;
			}
		}
		*/

	}  // loop over molecules
	_dYsum += this->getDeltaY();

	if (_dYsum >= _dSlabWidth)
	{
		global_log->info() << "_dSlabWidth = " << _dSlabWidth << endl;
		global_log->info() << "_dYsum = " << _dYsum << endl;
		global_log->info() << "_nSlabindex = " << _nSlabindex << endl;
		std::vector<Molecule> currentReservoirSlab = _reservoir.at(_nSlabindex);

		for(auto mi : currentReservoirSlab)
		{
			unsigned long tempId = mi.id();
			mi.setid(_maxId + tempId);
			mi.setComponent(comp2);
			mi.setr(1, mi.r(1) + _dYsum - _dSlabWidth);
			particleContainer->addParticle(mi);
		}

		_dYsum -= _dSlabWidth;
		_nSlabindex--;
		if(_nSlabindex == -1)
		{
			_nSlabindex = _reservoirSlabs-1;
			_maxId = _maxId + _reservoirNumMolecules;
		}
	}
	particleContainer->update();
}
void MettDeamon::postForce_action(ParticleContainer* particleContainer, DomainDecompBase* domainDecomposition)
{
	unsigned long nNumMoleculesLocal = 0;
	unsigned long nNumMoleculesGlobal = 0;

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM) {

		unsigned int cid = tM->componentid()+1;
		double v2 = tM->v2();
		if(cid != 2 && v2 > _velocityBarrier*_velocityBarrier) // v2_limit
		{
			unsigned long id = tM->id();
			double dY = tM->r(1);

			cout << "cid = " << cid << endl;
			cout << "id = " << id << endl;
			cout << "dY = " << dY << endl;
			cout << "v2 = " << v2 << endl;

			particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), false);
			_nNumMoleculesDeletedLocal++;
			_nNumMoleculesTooFast++;
		}
		if(cid == 2)
		{
			tM->setv(0, 0.);
			tM->setv(1, 0.);
			tM->setv(2, 0.);
		}

		// mirror, to simulate VLE
		if(true == _bMirrorActivated)
		{
			if(tM->r(1) >= _dMirrorPosY)
				tM->setv(1, -1.*tM->v(1) );
		}
	}
	particleContainer->update();
	nNumMoleculesLocal = particleContainer->getNumberOfParticles();

	// delta y berechnen: alle x Zeitschritte
	if(global_simulation->getSimulationStep() % _nUpdateFreq == 0)
	{
		// update global number of particles / calc global number of deleted particles
		domainDecomposition->collCommInit(2);
		domainDecomposition->collCommAppendUnsLong(nNumMoleculesLocal);
		domainDecomposition->collCommAppendUnsLong(_nNumMoleculesDeletedLocal);
		domainDecomposition->collCommAllreduceSum();
		nNumMoleculesGlobal = domainDecomposition->collCommGetUnsLong();
		_nNumMoleculesDeletedGlobal = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
		_nNumMoleculesDeletedGlobalAlltime += _nNumMoleculesDeletedGlobal;
		_nNumMoleculesDeletedLocal = 0;

		// update sum and summation list
		_numDeletedMolsSum += _nNumMoleculesDeletedGlobal;
		_numDeletedMolsSum -= _listDeletedMolecules.front();
		_listDeletedMolecules.push_back(_nNumMoleculesDeletedGlobal);
		if(_listDeletedMolecules.size() > _nNumValsSummation)
			_listDeletedMolecules.pop_front();
		else
		{
			_numDeletedMolsSum = 0;
			for(auto&& vi : _listDeletedMolecules)
				_numDeletedMolsSum += vi;
		}
		_dDeletedMolsPerTimestep = _numDeletedMolsSum * _dInvNumTimestepsSummation;
		this->calcDeltaY();
		global_log->info() << "_nNumMoleculesDeletedGlobal = " << _nNumMoleculesDeletedGlobal << endl;
		global_log->info() << "_numDeletedMolsSum = " << _numDeletedMolsSum << endl;
		global_log->info() << "_dDeletedMolsPerTimestep = " << _dDeletedMolsPerTimestep << endl;
		global_log->info() << "_dY = " << _dY << endl;
	}
	else
	{
		// update global number of particles
		domainDecomposition->collCommInit(1);
		domainDecomposition->collCommAppendUnsLong(nNumMoleculesLocal);
		domainDecomposition->collCommAllreduceSum();
		nNumMoleculesGlobal = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
	}
	global_simulation->getDomain()->setglobalNumMolecules(nNumMoleculesGlobal);

	// write restart file
	this->writeRestartfile();
}

void MettDeamon::writeRestartfile()
{
	if(0 != global_simulation->getSimulationStep() % _nWriteFreqRestart)
		return;

	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();

	if(domainDecomp.getRank() != 0)
		return;

	std::ofstream ofs("MettDeamonRestart.dat", std::ios::app);
	std::stringstream outputstream;

	outputstream << setw(12) << global_simulation->getSimulationStep() << setw(12) << _nSlabindex;
	outputstream << FORMAT_SCI_MAX_DIGITS << _dYsum << std::endl;

	ofs << outputstream.str();
	ofs.close();
}


