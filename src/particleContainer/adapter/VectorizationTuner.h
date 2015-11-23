/*
 * VectorizationTuner.h
 *
 *  Created on: Jun 12, 2015
 *      Author: tchipevn
 */

#ifndef VECTORIZATIONTUNER_H_
#define VECTORIZATIONTUNER_H_

#include "CellProcessor.h"
#include "FlopCounter.h"
#include <vector>

class Component;
class ParticleCell;
//class VectorizedCellProcessor;
//class FlopCounter;

class VectorizationTuner {
public:
	VectorizationTuner();
	~VectorizationTuner();

	void tune(std::vector<Component> ComponentList, CellProcessor& vcp, FlopCounter& fc);
    void iterate(std::vector<Component> ComponentList,  CellProcessor& vcp, FlopCounter& fc, int numMols, double& gflopsOwn, double& gflopsPair);
//, double& gflopsOwn, double& gflopsPair);

private:
	void runOwn(CellProcessor& cp, ParticleCell& cell1, int numRepetitions);
	void runPair(CellProcessor& cp, ParticleCell& cell1, ParticleCell& cell2, int numRepetitions);

	void initSomeMolecules(Component& comp, ParticleCell& cell1, ParticleCell& cell2);
	void initMeshOfMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2);
    void initUniformRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, int numMols);
    void initNormalRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, int numMols);
	void clearMolecules(ParticleCell& cell1);
};

#endif /* VECTORIZATIONTUNER_H_ */
