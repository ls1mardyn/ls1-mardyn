/*
 * MettDeamon.h
 *
 *  Created on: 03.04.2017
 *      Author: thet
 */

#ifndef METTDEAMON_H_
#define METTDEAMON_H_

#include "molecules/Molecule.h"

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <cstdint>
#include <limits>

#define FORMAT_SCI_MAX_DIGITS std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)

class Domain;
class Ensemble;
class DomainDecompBase;
class ParticleContainer;
class XMLfileUnits;

class MettDeamon
{
public:
	MettDeamon(double cutoffRadius);
	~MettDeamon();

	void readXML(XMLfileUnits& xmlconfig);
	double getAreaY(){return _dAreaXZ;}
	double getDeltaY() {return _dY;}
	int getnSlabindex() {return _nSlabindex;}
	double getdYsum() {return _dYsum;}

	void setnSlabindex(int slabindex) {_nSlabindex = slabindex;}
	void setdYsum(double dY) {_dYsum = dY;}

	uint64_t getnNumMoleculesDeleted( DomainDecompBase* domainDecomposition){return _nNumMoleculesDeletedGlobalAlltime;}
	uint64_t getnNumMoleculesDeleted2( DomainDecompBase* domainDecomposition);

	void prepare_start(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer, double cutoffRadius);
	void init_positionMap(ParticleContainer* particleContainer);
	void preForce_action(ParticleContainer* particleContainer, double cutoffRadius);
	void postForce_action(ParticleContainer* particleContainer, DomainDecompBase* domainDecomposition);

	// connection to DensityControl
	void IncrementDeletedMoleculesLocal() {_nNumMoleculesDeletedLocal++;}

private:
	void ReadReservoir(DomainDecompBase* domainDecomp);
	void writeRestartfile();
	void calcDeltaY() { _dY = _dDeletedMolsPerTimestep * _dInvDensityArea; }

private:
	double _rho_l;
	double _dAreaXZ;
	double _dInvDensityArea;
	double _dY;
	double _dYsum;
	double _velocityBarrier;
	double _dSlabWidthInit;
	double _dSlabWidth;
	double _cutoffRadius;
	uint64_t _nUpdateFreq;
	uint64_t _nWriteFreqRestart;
	uint64_t _maxId;
	uint64_t _nNumMoleculesDeletedLocal;
	uint64_t _nNumMoleculesDeletedGlobal;
	uint64_t _nNumMoleculesDeletedGlobalAlltime;
	uint64_t _nNumMoleculesTooFast;
	uint64_t _nNumMoleculesTooFastGlobal;
	uint64_t _reservoirNumMolecules;
	uint64_t _reservoirSlabs;
	int32_t _nSlabindex;
	std::string _reservoirFilename;
	std::map<uint64_t, std::array<double, 6> > _storePosition;  //Map for frozen particle position storage <"id, position">
	std::vector< std::vector<Molecule> >_reservoir;
	bool _bIsRestart;  // simulation is a restart?
	std::list<uint64_t> _listDeletedMolecules;
	uint32_t _nNumValsSummation;
	uint64_t _numDeletedMolsSum;
	double _dDeletedMolsPerTimestep;
	double _dInvNumTimestepsSummation;
	bool _bMirrorActivated;
	double _dMirrorPosY;
};

#endif /* METTDEAMON_H_ */
