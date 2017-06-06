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
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
//#include <string.h>

//preprocessing macros
//#define MASKING

VectorizationTuner::VectorizationTuner(std::string outputPrefix, unsigned int minMoleculeCnt, unsigned int maxMoleculeCnt,
		MoleculeCntIncreaseTypeEnum moleculeCntIncreaseType, double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor):
	_outputPrefix(outputPrefix), _minMoleculeCnt(minMoleculeCnt), _maxMoleculeCnt(maxMoleculeCnt), _moleculeCntIncreaseType(moleculeCntIncreaseType),
	_cellProcessor(cellProcessor), _cutoffRadius(cutoffRadius), _LJCutoffRadius(LJCutoffRadius), _flopCounterBigRc(NULL), _flopCounterNormalRc(NULL), _flopCounterZeroRc(NULL){

}

VectorizationTuner::VectorizationTuner(double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor):
	_outputPrefix("Mardyn"), _minMoleculeCnt(2), _maxMoleculeCnt(512), _moleculeCntIncreaseType(both),
	_cellProcessor(cellProcessor), _cutoffRadius(cutoffRadius), _LJCutoffRadius(LJCutoffRadius), _flopCounterBigRc(NULL), _flopCounterNormalRc(NULL), _flopCounterZeroRc(NULL){

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

void VectorizationTuner::initOutput(ParticleContainer* /*particleContainer*/,
			DomainDecompBase* /*domainDecomp*/, Domain* /*domain*/) {
	_flopCounterNormalRc = new FlopCounter(_cutoffRadius, _LJCutoffRadius);
	_flopCounterBigRc = new FlopCounter(_cutoffRadiusBig, _LJCutoffRadiusBig);
	_flopCounterZeroRc = new FlopCounter( 0., 0.);
	tune(*(_simulation.getEnsemble()->getComponents()));
}


using namespace std;

void VectorizationTuner::tune(std::vector<Component> ComponentList) {

	global_log->info() << "VT: begin VECTORIZATION TUNING "<< endl;

    double gflopsOwnBig=0., gflopsPairBig=0., gflopsOwnNormal=0., gflopsPairNormalFace=0., gflopsPairNormalEdge=0., gflopsPairNormalCorner=0.,gflopsOwnZero=0., gflopsPairZero=0.;

    string resultfile(_outputPrefix+".VT.csv");
    global_log->info() << "VT: Writing to file " << resultfile << endl;
    int rank = global_simulation->domainDecomposition().getRank();

    ofstream myfile;
    if (rank == 0) {
    	myfile.open(resultfile.c_str(), ofstream::out | ofstream::trunc);
    	myfile << "Vectorization Tuner File" << endl
			<< "The Cutoff Radii were: " << endl
			<< "NormalRc=" << _cutoffRadius << " , LJCutoffRadiusNormal=" << _LJCutoffRadius << endl
			<< "BigRC=" << _cutoffRadiusBig << " , BigLJCR=" << _LJCutoffRadiusBig << endl;
    }

    if(_moleculeCntIncreaseType==linear or _moleculeCntIncreaseType==both){
    	if (rank==0) {
			myfile << "Linearly distributed molecule counts" << endl;
			myfile << "Num. of Molecules, " << "Gflops for Own BigRc, " << "Gflops for Pair BigRc, " << "Gflops for Own NormalRc, " << "Gflops for Pair NormalRc Face, "
					<< "Gflops for Pair NormalRc Edge, "  << "Gflops for Pair NormalRc Corner, "  << "Gflops for Zero Rc (Own), " << "Gflops for Zero Rc (Pair)" << endl;
    	}
		for(unsigned int i = _minMoleculeCnt; i <= std::min(32u, _maxMoleculeCnt); i++){
			iterate(ComponentList, i,  gflopsOwnBig, gflopsPairBig, gflopsOwnNormal, gflopsPairNormalFace, gflopsPairNormalEdge, gflopsPairNormalCorner, gflopsOwnZero, gflopsPairZero);
			if (rank==0) {
				myfile << i << ", " << gflopsOwnBig << ", " << gflopsPairBig << ", " << gflopsOwnNormal << ", "
						<< gflopsPairNormalFace << ", " << gflopsPairNormalEdge << ", " << gflopsPairNormalCorner << ", "
						<< gflopsOwnZero << ", " << gflopsPairZero << endl;
			}
		}
		if (rank == 0) {
			myfile << endl;
		}
    }
    if(_moleculeCntIncreaseType==exponential or _moleculeCntIncreaseType==both){
    	if (rank == 0) {
			myfile << "Exponentially distributed molecule counts" << endl;
			myfile << "Num. of Molecules," << "Gflops for Own BigRc, " << "Gflops for Pair BigRc, " << "Gflops for Own NormalRc, " << "Gflops for Pair NormalRc Face, "
					<< "Gflops for Pair NormalRc Edge, "  << "Gflops for Pair NormalRc Corner, "  << "Gflops for Zero Rc (Own), " << "Gflops for Zero Rc (Pair)" << endl;
    	// logarithmically scaled axis -> exponentially increasing counts
    	}

    	for(unsigned int i = _minMoleculeCnt; i <= _maxMoleculeCnt; i*=2){
    		iterate(ComponentList, i, gflopsOwnBig, gflopsPairBig, gflopsOwnNormal, gflopsPairNormalFace, gflopsPairNormalEdge, gflopsPairNormalCorner, gflopsOwnZero, gflopsPairZero);
    		if (rank == 0) {
				myfile << i << ", " << gflopsOwnBig << ", " << gflopsPairBig << ", " << gflopsOwnNormal << ", "
									<< gflopsPairNormalFace << ", " << gflopsPairNormalEdge << ", " << gflopsPairNormalCorner << ", "
									<< gflopsOwnZero << ", " << gflopsPairZero << endl;
    		}
    	}
    }

    if (rank == 0) {
    	myfile.close();
    }

    _flopCounterZeroRc->resetCounters();
    _flopCounterBigRc->resetCounters();
    _flopCounterNormalRc->resetCounters();

	global_log->info() << "VECTORIZATION TUNING completed "<< endl;

}

void VectorizationTuner::iterateOwn(long long int numRepetitions,
		ParticleCell& cell, double& gflopsPair, FlopCounter& flopCounter) {
	runOwn(flopCounter, cell, 1);
	// run simulation for a pair of cells
	global_simulation->startTimer("VECTORIZATION_TUNER_TUNER");
	runOwn(**_cellProcessor, cell, numRepetitions);
	global_simulation->stopTimer("VECTORIZATION_TUNER_TUNER");
	// get Gflops for pair computations
	double tuningTime = global_simulation->getTime("VECTORIZATION_TUNER_TUNER");
	gflopsPair = flopCounter.getTotalFlopCount() * numRepetitions / tuningTime / (1024 * 1024 * 1024);
	global_log->info() << "FLOP-Count per Iteration: " << flopCounter.getTotalFlopCount() << " FLOPs" << endl;
	global_log->info() << "FLOP-rate: " << gflopsPair << " GFLOPS" << endl;
	global_log->info() << "number of iterations: " << numRepetitions << endl;
	global_log->info() << "total time: " << tuningTime << "s" << endl;
	global_log->info() << "time per iteration: " << tuningTime / numRepetitions << "s " << endl << endl;
	flopCounter.resetCounters();
	global_simulation->resetTimer("VECTORIZATION_TUNER_TUNER");
}

void VectorizationTuner::iteratePair(long long int numRepetitions, ParticleCell& firstCell,
		ParticleCell& secondCell, double& gflopsPair, FlopCounter& flopCounter) {
	//count/calculate the needed flops
	runPair(flopCounter, firstCell, secondCell, 1);
	// run simulation for a pair of cells
	global_simulation->startTimer("VECTORIZATION_TUNER_TUNER");
	runPair(**_cellProcessor, firstCell, secondCell, numRepetitions);
	global_simulation->stopTimer("VECTORIZATION_TUNER_TUNER");
	// get Gflops for pair computations
	double tuningTime = global_simulation->getTime("VECTORIZATION_TUNER_TUNER");
	gflopsPair = flopCounter.getTotalFlopCount() * numRepetitions / tuningTime / (1024 * 1024 * 1024);
	global_log->info() << "FLOP-Count per Iteration: " << flopCounter.getTotalFlopCount() << " FLOPs" << endl;
	global_log->info() << "FLOP-rate: " << gflopsPair << " GFLOPS" << endl;
	global_log->info() << "number of iterations: " << numRepetitions << endl;
	global_log->info() << "total time: " << tuningTime << "s" << endl;
	global_log->info() << "time per iteration: " << tuningTime / numRepetitions << "s " << endl << endl;
	flopCounter.resetCounters();
	global_simulation->resetTimer("VECTORIZATION_TUNER_TUNER");
}

void VectorizationTuner::iterate(std::vector<Component> ComponentList, unsigned int numMols, double& gflopsOwnBig,
		double& gflopsPairBig, double& /*gflopsOwnNormal*/, double& /*gflopsPairNormalFace*/, double& /*gflopsPairNormalEdge*/,
		double& /*gflopsPairNormalCorner*/, double& gflopsOwnZero, double& gflopsPairZero) {


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
	double dirxplus[3] = { 1., 0., 0.};
	double BoxMin2[3] = { 1., 0., 0. };
	double BoxMax2[3] = { 2., 1., 1. };

	firstCell.setBoxMin(BoxMin);
	secondCell.setBoxMin(BoxMin2);

	firstCell.setBoxMax(BoxMax);
	secondCell.setBoxMax(BoxMax2);

	//double diryplus[3] = {0., 1., 0.};
	//double dirzplus[3] = {0., 0., 1.};

	global_log->info() << "--------------------------Molecule count: " << numMols << "--------------------------" << endl;

	//initialize both cells with molecules between 0,0,0 and 1,1,1
    initUniformRandomMolecules(BoxMin, BoxMax, comp, firstCell, numMols);
    initUniformRandomMolecules(BoxMin2, BoxMax2, comp, secondCell, numMols);
    //moveMolecules(dirxplus, secondCell);

	firstCell.buildSoACaches();
	secondCell.buildSoACaches();

	long long int numRepetitions = std::max(20000000u / (numMols*numMols), 10u);


	//0a,0b: 0RC
		(**_cellProcessor).setCutoffRadius(0.);
		(**_cellProcessor).setLJCutoffRadius(0.);
		iterateOwn(numRepetitions, firstCell, gflopsOwnZero, *_flopCounterZeroRc);
		iteratePair(numRepetitions, firstCell, secondCell, gflopsPairZero, *_flopCounterZeroRc);
    //1+2: bigRC
	(**_cellProcessor).setCutoffRadius(_cutoffRadiusBig);
	(**_cellProcessor).setLJCutoffRadius(_LJCutoffRadiusBig);
    //1. own, bigRC
		iterateOwn(numRepetitions, firstCell, gflopsOwnBig, *_flopCounterBigRc);
	//2. pair, bigRC
		iteratePair(numRepetitions, firstCell, secondCell, gflopsPairBig, *_flopCounterBigRc);
#if 0
		TODO: redo these with mesh of molecules
	//3,...: normalRC
	(**_cellProcessor).setCutoffRadius(_cutoffRadius);
	(**_cellProcessor).setLJCutoffRadius(_LJCutoffRadius);
	//3. own, normalRC
		iterateOwn(numRepetitions, firstCell, gflopsOwnNormal, *_flopCounterNormalRc);
	//4. pair, normalRC face
		iteratePair(numRepetitions, firstCell, secondCell, gflopsPairNormalFace, *_flopCounterNormalRc); //cell2s particles moved by 1,0,0 - common face
	//5. pair, normalRC edge
		moveMolecules(diryplus, secondCell);
		iteratePair(numRepetitions, firstCell, secondCell, gflopsPairNormalEdge, *_flopCounterNormalRc); //cell2s particles moved by 1,1,0 - common edge
	//6. pair, normalRC corner
		moveMolecules(dirzplus, secondCell);
		iteratePair(numRepetitions, firstCell, secondCell, gflopsPairNormalCorner, *_flopCounterNormalRc); //cell2s particles moved by 1,1,1 - common corner
#endif

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

	cp.initTraversal();

	cp.preprocessCell(cell1);

	for (int i = 0; i < numRepetitions; ++i) {
		cp.processCell(cell1);
	}

	cp.postprocessCell(cell1);
	cp.endTraversal();
}

void VectorizationTuner::runPair(CellProcessor& cp, ParticleCell& cell1, ParticleCell& cell2, int numRepetitions) {


	cp.initTraversal();

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
	cell.deallocateAllParticles();
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

				Molecule m = Molecule(
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
	}

    id = 0;
	for(int z = 0; z < numMoleculesZ; ++z) {
		for(int y = 0; y < numMoleculesY; ++y) {
			for(int x = 0; x < numMoleculesX; ++x) {

				pos[0] = start_pos2[0] + x*dx;
				pos[1] = start_pos2[1] + y*dx;
				pos[2] = start_pos2[2] + z*dx;

				Molecule m = Molecule(
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

		Molecule m = Molecule(
				id, &comp,
				pos[0], pos[1], pos[2],
				vel[0], vel[1], vel[2],
				orientation[0], orientation[1], orientation[2], orientation[3],
				angularVelocity[0], angularVelocity[1], angularVelocity[2]
				);
		cell.addParticle(m);
		id++; // id's need to be distinct
		//global_log->info() << pos[0] << " " << pos[1] << " " << pos[2] << endl;
	}
}


void VectorizationTuner::initNormalRandomMolecules(double /*boxMin*/[3], double /*boxMax*/[3], Component& comp,
		ParticleCell& cell1, ParticleCell& /*cell2*/, unsigned int numMols) {
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

		Molecule m = Molecule(
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
		Molecule& mol = cell.moleculesAt(i);
		mol.move(0, direction[0]);
		mol.move(1, direction[1]);
		mol.move(2, direction[2]);
		//global_log->info() << mol->r(0) << " " << mol->r(1) << " " << mol->r(2) << endl;
	}
}
