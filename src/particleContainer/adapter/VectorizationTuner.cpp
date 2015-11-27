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

VectorizationTuner::VectorizationTuner(std::string outputPrefix, double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor):
_outputPrefix(outputPrefix), _cellProcessor(cellProcessor), _cutoffRadius(cutoffRadius), _LJCutoffRadius(LJCutoffRadius){

}

VectorizationTuner::VectorizationTuner(double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor):
		_outputPrefix("Mardyn"), _cellProcessor(cellProcessor), _cutoffRadius(cutoffRadius), _LJCutoffRadius(LJCutoffRadius){

}

VectorizationTuner::~VectorizationTuner() {
	// TODO Auto-generated destructor stub
}

void VectorizationTuner::initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) {
	FlopCounter flopCounter2 = FlopCounter(_cutoffRadius, _LJCutoffRadius);
	tune(*(_simulation.getEnsemble()->components()), **_cellProcessor , flopCounter2);
}


using namespace std;

void VectorizationTuner::tune(std::vector<Component> ComponentList, CellProcessor& vcp, FlopCounter& fc) {

	global_log->info() << "begin VECTORIZATION TUNING "<< endl;

    double gflopsOwn, gflopsPair;

    stringstream filenamestream;
	filenamestream << _outputPrefix;
	filenamestream << ".VT.csv";
    char const* value = filenamestream.str().c_str();
//    const char fileName = atoi(value);

    ofstream myfile;
//    myfile.open("./stats/Gflops.csv", std::ios::trunc);
    myfile.open(value, ios::trunc);
    myfile.close();
    myfile.open(value, ofstream::in | ofstream::out | ofstream::app);

    myfile << "Linearly" << endl;
    myfile << "Num. of Molecules," << "Gflops for Own," << "Gflops for Pair" << endl;
    // linearly
    for(int i = 2; i <= 32; i++){
        iterate(ComponentList, vcp, fc, i,  gflopsOwn, gflopsPair);
        myfile << i << ", " << gflopsOwn << ", " << gflopsPair << endl;
    }
    myfile << endl;
    myfile << "Logarithmically" << endl;
    myfile << "Num. of Molecules," << "Gflops for Own," << "Gflops for Pair" << endl;
    // logarithmically
    for(int i = 2; i <= 512; i*=2){
        iterate(ComponentList, vcp, fc, i, gflopsOwn, gflopsPair);
        myfile << i << ", " << gflopsOwn << ", " << gflopsPair << endl;
    }

    myfile.close();

    fc.resetCounters();

	global_log->info() << "VECTORIZATION TUNING completed "<< endl;

}

void VectorizationTuner::iterate(std::vector<Component> ComponentList, CellProcessor& vcp, FlopCounter& fc, int numMols, double& gflopsOwn, double& gflopsPair){
//, double& gflopsOwn, double& gflopsPair) {

	// get (first) component
	Component comp = ComponentList[0];

	// construct two cells
	ParticleCell first;
	ParticleCell second;
	first.assignCellToInnerRegion();
	second.assignCellToInnerRegion();

    #ifdef MASKING
    srand(time(NULL));
    #else
    srand(5);
    #endif

	double BoxMin[3] = {0.0, 0.0, 0.0};
	double BoxMax[3] = {1.0, 1.0, 1.0};

	Timer timer;
#ifdef ENABLE_MPI
	timer.set_sync(false);
#endif

//	 initialization

//  initSomeMolecules(comp, first, second);
//  initMeshOfMolecules(BoxMin, BoxMax, comp, first, second);
    initUniformRandomMolecules(BoxMin, BoxMax, comp, first, second, numMols);
//  initNormalRandomMolecules(BoxMin, BoxMax, comp, first, second, numMols);

	runOwn(fc, first, 1);

	// repeat many times and measure average time
	//	long long int numRepetitions = 10000000; // TODO: for realistic measurements (at least for the case with 8 molecules)
	long long int numRepetitions = 10000;

    // run simulation for a single cell
	timer.start();
	runOwn(vcp, first, numRepetitions);
	timer.stop();

    // get Gflops for own computations
	gflopsOwn = fc.getTotalFlopCount() * numRepetitions / timer.get_etime() / (1024*1024*1024);
	global_log->info() << "FLOP-Count per Iteration: " << fc.getTotalFlopCount() << " FLOPs" << endl;
	global_log->info() << "FLOP-rate: " << gflopsOwn << " GFLOPS" << endl;

    fc.resetCounters();
    timer.reset();


	runPair(fc, first, second, 1);

    // run simulation for a pair of cells
	timer.start();
	runPair(vcp, first, second, numRepetitions);
	timer.stop();

    // get Gflops for pair computations
	gflopsPair = fc.getTotalFlopCount() * numRepetitions / timer.get_etime() / (1024*1024*1024);
	global_log->info() << "FLOP-Count per Iteration: " << fc.getTotalFlopCount() << " FLOPs" << endl;
	global_log->info() << "FLOP-rate: " << gflopsPair << " GFLOPS" << endl;

    fc.resetCounters();
    timer.reset();

	// clear cells
	clearMolecules(first);
	clearMolecules(second);

}


double uniformRandom()
{
  return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}
 // return a normally distributed random number
double normalRandom()
{
  double u1=uniformRandom();
  double u2=uniformRandom();
  return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1)); 
}


// pass also second cell as argument, when it comes to that
void VectorizationTuner::runOwn(CellProcessor& cp, ParticleCell& cell1, int numRepetitions) {

	cp.initTraversal(1);

	cp.preprocessCell(cell1);

	for (int i = 0; i < numRepetitions; ++i) {
		cp.processCell(cell1);
//        cp.processCellPair(cell1, cell2);
	}

	cp.postprocessCell(cell1);
	cp.endTraversal();
}

void VectorizationTuner::runPair(CellProcessor& cp, ParticleCell& cell1, ParticleCell& cell2, int numRepetitions) {


	cp.initTraversal(2);

	cp.preprocessCell(cell1);
	cp.preprocessCell(cell2);

	for (int i = 0; i < numRepetitions; ++i) {
//		cp.processCell(cell1);
        cp.processCellPair(cell1, cell2);
	}

	cp.postprocessCell(cell1);
	cp.postprocessCell(cell2);
	cp.endTraversal();

}

// initialize eight molecules in the middle of the cell
void VectorizationTuner::initSomeMolecules(Component& comp, ParticleCell& cell1, ParticleCell& cell2) {

	unsigned long id = 0;
	double vel[3] = { 0.0, 0.0, 0.0 };
	double orientation[4] = { 1.0, 0.0, 0.0, 0.0 }; // NOTE the 1.0 in the first coordinate
	double angularVelocity[3] = { 0.0, 0.0, 0.0 };

	double start_pos1[3] = { 0.25, 0.25, 0.25 };
	double start_pos2[3] = { 1.25, 0.25, 0.25 };
	double pos[3];


	for (int z = 0; z < 2; ++z) {
		for (int y = 0; y < 2; ++y) {
			for (int x = 0; x < 2; ++x) {
				pos[0] = start_pos1[0] + x*0.5;
				pos[1] = start_pos1[1] + y*0.5;
				pos[2] = start_pos1[2] + z*0.5;

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
	}
	
    id = 0;
	for (int z = 0; z < 2; ++z) {
		for (int y = 0; y < 2; ++y) {
			for (int x = 0; x < 2; ++x) {
				pos[0] = start_pos2[0] + x*0.5;
				pos[1] = start_pos2[1] + y*0.5;
				pos[2] = start_pos2[2] + z*0.5;

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
//	global_log->info() << "pos0 =  " << pos[0] << endl;
//	global_log->info() << "pos1 =  " << pos[1] << endl;
//	global_log->info() << "pos2 =  " << pos[2] << endl;
			}
		}
	}
}

void VectorizationTuner::initUniformRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, int numMols) {


//    char* value;
//    value = getenv("VT_N");
//    int numMolecules = atoi(value);


//	int numMolecules = 2;
//  int numMoleculesY = 2, numMoleculesZ = 2;

	unsigned long id = 0;
	double vel[3] = { 0.0, 0.0, 0.0 };
	double orientation[4] = { 1.0, 0.0, 0.0, 0.0 }; // NOTE the 1.0 in the first coordinate
	double angularVelocity[3] = { 0.0, 0.0, 0.0 };

	double start_pos1[3] = { 0.0, 0.0, 0.0 }; // adapt for pair (Bmax Bmin)
	double start_pos2[3] = { 1.0, 0.0, 0.0 };
	double pos[3];

//	double dx = boxMax[0] - boxMin[0] / numMoleculesX;
//  double pos_randx = ((double)rand()/(double)RAND_MAX);
//  double pos_randy = ((double)rand()/(double)RAND_MAX);
//  double pos_randz = ((double)rand()/(double)RAND_MAX);

	for(int i = 0; i < numMols; ++i) {


				pos[0] = start_pos1[0] + ((double)rand()/(double)RAND_MAX);
				pos[1] = start_pos1[1] + ((double)rand()/(double)RAND_MAX);
				pos[2] = start_pos1[2] + ((double)rand()/(double)RAND_MAX);

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

    id = 0;
	for(int i = 0; i < numMols; ++i) {


				pos[0] = start_pos2[0] + ((double)rand()/(double)RAND_MAX);
				pos[1] = start_pos2[1] + ((double)rand()/(double)RAND_MAX);
				pos[2] = start_pos2[2] + ((double)rand()/(double)RAND_MAX);

				Molecule* m = new Molecule(
						id, &comp,
						pos[0], pos[1], pos[2],
						vel[0], vel[1], vel[2],
						orientation[0], orientation[1], orientation[2], orientation[3],
						angularVelocity[0], angularVelocity[1], angularVelocity[2]
						);

				cell2.addParticle(m);
				id++; // id's need to be distinct
//	global_log->info() << "pos0 =  " << pos[0] << endl;
//	global_log->info() << "pos1 =  " << pos[1] << endl;
//	global_log->info() << "pos2 =  " << pos[2] << endl;
	}
}


void VectorizationTuner::initNormalRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, int numMols) {

//	int numMoleculesX = 20;
//  int numMoleculesY = 2, numMoleculesZ = 2;

	unsigned long id = 0;
	double vel[3] = { 0.0, 0.0, 0.0 };
	double orientation[4] = { 1.0, 0.0, 0.0, 0.0 }; // NOTE the 1.0 in the first coordinate
	double angularVelocity[3] = { 0.0, 0.0, 0.0 };

	double start_pos[3] = { 0.0, 0.0, 0.0 }; // adapt for pair (Bmax Bmin)
	double pos[3];

//	double dx = boxMax[0] - boxMin[0] / numMoleculesX;
//  double pos_randx = ((double)rand()/(double)RAND_MAX);
//  double pos_randy = ((double)rand()/(double)RAND_MAX);
//  double pos_randz = ((double)rand()/(double)RAND_MAX);

	for(int i = 0; i < numMols; ++i) {


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
	global_log->info() << "pos0 =  " << pos[0] << endl;
	global_log->info() << "pos1 =  " << pos[1] << endl;
	global_log->info() << "pos2 =  " << pos[2] << endl;
	}
}
