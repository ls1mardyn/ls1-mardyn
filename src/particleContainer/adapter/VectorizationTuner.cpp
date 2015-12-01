/*
 * VectorizationTuner.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: tchipevn
 */

#include "VectorizationTuner.h"

#include <vector>

#include "ensemble/EnsembleBase.h"
#include "molecules/Component.h"
#include "Simulation.h"
#include "particleContainer/ParticleCell.h"
#include "CellDataSoA.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "utils/Timer.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
//#include <string.h>

//preprocessing macros
//#define MASKING

VectorizationTuner::VectorizationTuner(std::string outputPrefix, unsigned int minMoleculeCnt, unsigned int maxMoleculeCnt,
		MoleculeCntIncreaseTypeEnum moleculeCntIncreaseType, double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor):
	_outputPrefix(outputPrefix), _minMoleculeCnt(minMoleculeCnt), _maxMoleculeCnt(maxMoleculeCnt), _moleculeCntIncreaseType(moleculeCntIncreaseType),
	_cellProcessor(cellProcessor), _cutoffRadius(cutoffRadius), _LJCutoffRadius(LJCutoffRadius), _flopCounterBigRc(NULL), _flopCounterNormalRc(NULL){

}

VectorizationTuner::VectorizationTuner(double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor):
	_outputPrefix("Mardyn"), _minMoleculeCnt(2), _maxMoleculeCnt(512), _moleculeCntIncreaseType(both),
	_cellProcessor(cellProcessor), _cutoffRadius(cutoffRadius), _LJCutoffRadius(LJCutoffRadius), _flopCounterBigRc(NULL), _flopCounterNormalRc(NULL){

}

VectorizationTuner::~VectorizationTuner() {
	// TODO Auto-generated destructor stub
}

void VectorizationTuner::readXML(XMLfileUnits& xmlconfig) {
	_outputPrefix = "mardyn";

	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << std::endl;

	xmlconfig.getNodeValue("minmoleculecnt", _minMoleculeCnt);
	global_log->info() << "Minimal molecule count: " << _minMoleculeCnt << std::endl;

	xmlconfig.getNodeValue("maxmoleculecnt", _maxMoleculeCnt);
	global_log->info() << "Maximal molecule count: " << _maxMoleculeCnt << std::endl;

	int type=2;
	xmlconfig.getNodeValue("moleculecntincreasetype", type);
	_moleculeCntIncreaseType = static_cast<MoleculeCntIncreaseTypeEnum>(type);
	global_log->info() << "Molecule count increase type: " << _moleculeCntIncreaseType << std::endl;

}

void VectorizationTuner::initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) {
	_flopCounterNormalRc = new FlopCounter(_cutoffRadius, _LJCutoffRadius);
	_flopCounterBigRc = new FlopCounter(_cutoffRadiusBig, _LJCutoffRadiusBig);
	tune(*(_simulation.getEnsemble()->components()));
}


using namespace std;

void VectorizationTuner::tune(std::vector<Component> ComponentList) {

	global_log->info() << "VT: begin VECTORIZATION TUNING "<< endl;

    double gflopsOwnBig, gflopsPairBig, gflopsOwnNormal, gflopsPairNormalFace, gflopsPairNormalEdge, gflopsPairNormalCorner;

    stringstream filenamestream;
	filenamestream << _outputPrefix;
	filenamestream << ".VT.csv";
    char const* value = filenamestream.str().c_str();

    global_log->info() << "VT: Writing to file " << value << endl;
    ofstream myfile;
    myfile.open(value, ofstream::out | ofstream::trunc);
    myfile << "Vectorization Tuner File" << endl
    		<< "The Cutoff Radii were: " << endl << "NormalRc=" << _cutoffRadius << " , LJCutoffRadiusNormal=" << _LJCutoffRadius << endl
			<< " , BigRC=" << _cutoffRadiusBig << " , BigLJCR=" << _LJCutoffRadiusBig << endl;

    if(_moleculeCntIncreaseType==linear or _moleculeCntIncreaseType==both){
		myfile << "Linearly distributed molecule counts" << endl;
		myfile << "Num. of Molecules, " << "Gflops for Own BigRc, " << "Gflops for Pair BigRc, " << "Gflops for Own NormalRc, " << "Gflops for Pair NormalRc Face, "
		    			<< "Gflops for Pair NormalRc Edge, "  << "Gflops for Pair NormalRc Corner"  << endl;
		for(unsigned int i = _minMoleculeCnt; i <= std::min(32u, _maxMoleculeCnt); i++){
			iterate(ComponentList, i,  gflopsOwnBig, gflopsPairBig, gflopsOwnNormal, gflopsPairNormalFace, gflopsPairNormalEdge, gflopsPairNormalCorner);
			myfile << i << ", " << gflopsOwnBig << ", " << gflopsPairBig << ", " << gflopsOwnNormal << ", " << gflopsPairNormalFace << ", " << gflopsPairNormalEdge << ", " << gflopsPairNormalCorner << endl;
		}
		myfile << endl;
    }
    if(_moleculeCntIncreaseType==exponential or _moleculeCntIncreaseType==both){
    	myfile << "Exponentially distributed molecule counts" << endl;
    	myfile << "Num. of Molecules," << "Gflops for Own BigRc, " << "Gflops for Pair BigRc, " << "Gflops for Own NormalRc, " << "Gflops for Pair NormalRc Face, "
    			<< "Gflops for Pair NormalRc Edge, "  << "Gflops for Pair NormalRc Corner"  << endl;
    	// logarithmically scaled axis -> exponentially increasing counts

    	for(unsigned int i = _minMoleculeCnt; i <= _maxMoleculeCnt; i*=2){
    		iterate(ComponentList, i, gflopsOwnBig, gflopsPairBig, gflopsOwnNormal, gflopsPairNormalFace, gflopsPairNormalEdge, gflopsPairNormalCorner);
    		myfile << i << ", " << gflopsOwnBig << ", " << gflopsPairBig << ", " << gflopsOwnNormal << ", " << gflopsPairNormalFace << ", " << gflopsPairNormalEdge << ", " << gflopsPairNormalCorner << endl;
    	}
    }
    myfile.close();

    _flopCounterBigRc->resetCounters();
    _flopCounterNormalRc->resetCounters();

	global_log->info() << "VECTORIZATION TUNING completed "<< endl;

}

void VectorizationTuner::iterateOwn(Timer timer, long long int numRepetitions,
		ParticleCell& cell, double& gflopsPair, FlopCounter& flopCounter) {
	runOwn(flopCounter, cell, 1);
	// run simulation for a pair of cells
	timer.start();
	runOwn(**_cellProcessor, cell, numRepetitions);
	timer.stop();
	// get Gflops for pair computations
	gflopsPair = flopCounter.getTotalFlopCount() * numRepetitions
			/ timer.get_etime() / (1024 * 1024 * 1024);
	global_log->info() << "FLOP-Count per Iteration: "
			<< flopCounter.getTotalFlopCount() << " FLOPs" << endl;
	global_log->info() << "FLOP-rate: " << gflopsPair << " GFLOPS" << endl << endl;
	flopCounter.resetCounters();
	timer.reset();
}

void VectorizationTuner::iteratePair(Timer timer, long long int numRepetitions,
		ParticleCell& firstCell, ParticleCell& secondCell, double& gflopsPair, FlopCounter& flopCounter) {
	//count/calculate the needed flops
	runPair(flopCounter, firstCell, secondCell, 1);
	// run simulation for a pair of cells
	timer.start();
	runPair(**_cellProcessor, firstCell, secondCell, numRepetitions);
	timer.stop();
	// get Gflops for pair computations
	gflopsPair = flopCounter.getTotalFlopCount() * numRepetitions
			/ timer.get_etime() / (1024 * 1024 * 1024);
	global_log->info() << "FLOP-Count per Iteration: "
			<< flopCounter.getTotalFlopCount() << " FLOPs" << endl;
	global_log->info() << "FLOP-rate: " << gflopsPair << " GFLOPS" << endl << endl;
	flopCounter.resetCounters();
	timer.reset();
}

void VectorizationTuner::iterate(std::vector<Component> ComponentList, unsigned int numMols, double& gflopsOwnBig, double& gflopsPairBig, double& gflopsOwnNormal, double& gflopsPairNormalFace,
		double& gflopsPairNormalEdge, double& gflopsPairNormalCorner){


	// get (first) component
	Component comp = ComponentList[0];

	// construct two cells
	ParticleCell firstCell;
	ParticleCell secondCell;
	firstCell.assignCellToInnerRegion();
	secondCell.assignCellToInnerRegion();

    #ifdef MASKING
    srand(time(NULL));
    #else
    srand(5);//much random, much wow :D
    #endif

	double BoxMin[3] = {0., 0., 0.};
	double BoxMax[3] = {1., 1., 1.};
	double dirxplus[3] = {1., 0., 0.};
	double diryplus[3] = {0., 1., 0.};
	double dirzplus[3] = {0., 0., 1.};

	Timer timer;
#ifdef ENABLE_MPI
	timer.set_sync(false);
#endif

	//initialize both cells with molecules between 0,0,0 and 1,1,1
    initUniformRandomMolecules(BoxMin, BoxMax, comp, firstCell, numMols);
    initUniformRandomMolecules(BoxMin, BoxMax, comp, secondCell, numMols);

	long long int numRepetitions = 10000;

	global_log->info() << "--------------------------Molecule count: " << numMols << "--------------------------" << endl;
    //1+2: bigRC
	(**_cellProcessor).setCutoffRadius(_cutoffRadiusBig);
	(**_cellProcessor).setLJCutoffRadius(_LJCutoffRadiusBig);
    //1. own, bigRC
		iterateOwn(timer, numRepetitions, firstCell, gflopsOwnBig, *_flopCounterBigRc);
	//2. pair, bigRC
		moveMolecules(dirxplus, secondCell);
		iteratePair(timer, numRepetitions, firstCell, secondCell, gflopsPairBig, *_flopCounterBigRc);
	//3,...: normalRC
	(**_cellProcessor).setCutoffRadius(_cutoffRadius);
	(**_cellProcessor).setLJCutoffRadius(_LJCutoffRadius);
	//3. own, normalRC
		iterateOwn(timer, numRepetitions, firstCell, gflopsOwnNormal, *_flopCounterNormalRc);
	//4. pair, normalRC face
		iteratePair(timer, numRepetitions, firstCell, secondCell, gflopsPairNormalFace, *_flopCounterNormalRc); //cell2s particles moved by 1,0,0 - common face
	//5. pair, normalRC edge
		moveMolecules(diryplus, secondCell);
		iteratePair(timer, numRepetitions, firstCell, secondCell, gflopsPairNormalEdge, *_flopCounterNormalRc); //cell2s particles moved by 1,1,0 - common edge
	//6. pair, normalRC corner
		moveMolecules(dirzplus, secondCell);
		iteratePair(timer, numRepetitions, firstCell, secondCell, gflopsPairNormalCorner, *_flopCounterNormalRc); //cell2s particles moved by 1,1,1 - common corner

	// clear cells
	clearMolecules(firstCell);
	clearMolecules(secondCell);

}

// returns a uniformly distributed random number between zero and one. (zero, one excluded)
double uniformRandom()
{
  return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 2. );
}

 // return a normally distributed random number (uses box-muller transform)
double normalRandom()
{
  double u1=uniformRandom();
  double u2=uniformRandom();
  return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1)); //box muller transform: cos(2*pi*u2)*sqrt(-2ln(u1))
}


// pass also second cell as argument, when it comes to that
void VectorizationTuner::runOwn(CellProcessor& cp, ParticleCell& cell1, int numRepetitions) {

	cp.initTraversal(1);

	cp.preprocessCell(cell1);

	for (int i = 0; i < numRepetitions; ++i) {
		cp.processCell(cell1);
	}

	cp.postprocessCell(cell1);
	cp.endTraversal();
}

void VectorizationTuner::runPair(CellProcessor& cp, ParticleCell& cell1, ParticleCell& cell2, int numRepetitions) {


	cp.initTraversal(2);

	cp.preprocessCell(cell1);
	cp.preprocessCell(cell2);

	for (int i = 0; i < numRepetitions; ++i) {
        cp.processCellPair(cell1, cell2);
	}

	cp.postprocessCell(cell1);
	cp.postprocessCell(cell2);
	cp.endTraversal();

}

// should work also if molecules are initialized via initMeshOfMolecules, initUniformRandomMolecules, initNormalRandomMolecules
void VectorizationTuner::clearMolecules(ParticleCell & cell) {
	std::vector<Molecule*>& cellMolecules = cell.getParticlePointers();

	int numMolecules = cellMolecules.size();
	for(int i = 0; i < numMolecules; ++i) {
		delete cellMolecules[i];
	}

	cell.removeAllParticles();
}


void VectorizationTuner::initMeshOfMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2) {
	//TODO: Bmax, Bmin
	int numMoleculesX = 3, numMoleculesY = 3, numMoleculesZ = 3;

	unsigned long id = 0;
	double vel[3] = { 0.0, 0.0, 0.0 };
	double orientation[4] = { 1.0, 0.0, 0.0, 0.0 }; // NOTE the 1.0 in the first coordinate
	double angularVelocity[3] = { 0.0, 0.0, 0.0 };

	double start_pos1[3] = { 0.0, 0.0, 0.0 }; // adapt for pair (Bmax Bmin)
	double start_pos2[3] = { 1.0, 0.0, 0.0 };
	double pos[3];


	double dx = (boxMax[0] - boxMin[0]) / numMoleculesX;

	for(int z = 0; z < numMoleculesZ; ++z) {
		for(int y = 0; y < numMoleculesY; ++y) {
			for(int x = 0; x < numMoleculesX; ++x) {

				pos[0] = start_pos1[0] + x*dx;
				pos[1] = start_pos1[1] + y*dx;
				pos[2] = start_pos1[2] + z*dx;

				Molecule* m = new Molecule(
						id, &comp,
						pos[0], pos[1], pos[2],
						vel[0], vel[1], vel[2],
						orientation[0], orientation[1], orientation[2], orientation[3],
						angularVelocity[0], angularVelocity[1], angularVelocity[2]
						);

				cell1.addParticle(m);
				id++; // id's need to be distinct
//	global_log->info() << "pos0 =  " << pos[0] << endl;
//	global_log->info() << "pos1 =  " << pos[1] << endl;
//	global_log->info() << "pos2 =  " << pos[2] << endl;
			}
		}
	}

    id = 0;
	for(int z = 0; z < numMoleculesZ; ++z) {
		for(int y = 0; y < numMoleculesY; ++y) {
			for(int x = 0; x < numMoleculesX; ++x) {

				pos[0] = start_pos2[0] + x*dx;
				pos[1] = start_pos2[1] + y*dx;
				pos[2] = start_pos2[2] + z*dx;

				Molecule* m = new Molecule(
						id, &comp,
						pos[0], pos[1], pos[2],
						vel[0], vel[1], vel[2],
						orientation[0], orientation[1], orientation[2], orientation[3],
						angularVelocity[0], angularVelocity[1], angularVelocity[2]
						);

				cell2.addParticle(m);
				id++; // id's need to be distinct
			}
		}
	}
}

void VectorizationTuner::initUniformRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell, unsigned int numMols) {
	unsigned long id = 0;
	double vel[3] = { 0.0, 0.0, 0.0 };
	double orientation[4] = { 1.0, 0.0, 0.0, 0.0 }; // NOTE the 1.0 in the first coordinate
	double angularVelocity[3] = { 0.0, 0.0, 0.0 };

	double pos[3];

	for(unsigned int i = 0; i < numMols; ++i) {
		pos[0] = boxMin[0] + ((double)rand()/(double)RAND_MAX)*(boxMax[0] - boxMin[0]);
		pos[1] = boxMin[1] + ((double)rand()/(double)RAND_MAX)*(boxMax[1] - boxMin[1]);
		pos[2] = boxMin[2] + ((double)rand()/(double)RAND_MAX)*(boxMax[2] - boxMin[2]);

		Molecule* m = new Molecule(
				id, &comp,
				pos[0], pos[1], pos[2],
				vel[0], vel[1], vel[2],
				orientation[0], orientation[1], orientation[2], orientation[3],
				angularVelocity[0], angularVelocity[1], angularVelocity[2]
				);
		cell.addParticle(m);
		id++; // id's need to be distinct
	}
}


void VectorizationTuner::initNormalRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, unsigned int numMols) {
//TODO: currently only cell 1
//TODO: does not really have/need/can_use Bmax, Bmin - is normal dist. proper???

	unsigned long id = 0;
	double vel[3] = { 0.0, 0.0, 0.0 };
	double orientation[4] = { 1.0, 0.0, 0.0, 0.0 }; // NOTE the 1.0 in the first coordinate
	double angularVelocity[3] = { 0.0, 0.0, 0.0 };

	double start_pos[3] = { 0.0, 0.0, 0.0 }; // adapt for pair (Bmax Bmin)
	double pos[3];

//	double dx = boxMax[0] - boxMin[0] / numMoleculesX;


	for(unsigned int i = 0; i < numMols; ++i) {


		pos[0] = normalRandom();
		pos[1] = start_pos[1] + normalRandom();
		pos[2] = start_pos[2] + normalRandom();

		Molecule* m = new Molecule(
				id, &comp,
				pos[0], pos[1], pos[2],
				vel[0], vel[1], vel[2],
				orientation[0], orientation[1], orientation[2], orientation[3],
				angularVelocity[0], angularVelocity[1], angularVelocity[2]
				);

		cell1.addParticle(m);
		id++; // id's need to be distinct
	}
}

void VectorizationTuner::moveMolecules(double direction[3], ParticleCell& cell){
	unsigned int cnt=cell.getMoleculeCount();
	for(unsigned int i=0; i < cnt; ++i){
		Molecule* mol = cell.getParticlePointers().at(i);
		mol->move(0, direction[0]);
		mol->move(0, direction[1]);
		mol->move(0, direction[2]);
	}
}
