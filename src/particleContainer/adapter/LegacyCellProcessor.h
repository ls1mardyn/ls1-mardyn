/*
 * LegacyCellProcessor.h
 *
 * @Date: 18.03.2012
 * @Author: eckhardw
 */

#ifndef LEGACYCELLPROCESSOR_H_
#define LEGACYCELLPROCESSOR_H_

#include "particleContainer/adapter/CellProcessor.h"

class ParticlePairsHandler;
#include "particleContainer/ParticleCellForwardDeclaration.h"

/**
 * This class is a simple extraction of the cell handling from the old
 * class Blocktraverse.
 * I expect it to be slightly refactored one day...
 */
class LegacyCellProcessor : public CellProcessor {

private:
	//const double _cutoffRadiusSquare;
	//const double _LJCutoffRadiusSquare;
	ParticlePairsHandler* const _particlePairsHandler;

public:
	LegacyCellProcessor& operator=(const LegacyCellProcessor&) = delete;

	LegacyCellProcessor(const double cutoffRadius, const double LJCutoffRadius,
			ParticlePairsHandler* particlePairsHandler);

	virtual ~LegacyCellProcessor();

	void initTraversal();

	void preprocessCell(ParticleCell& /*cell*/) {}

        void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false);

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2);

	void processCell(ParticleCell& cell);

	void postprocessCell(ParticleCell& /*cell*/) {}

	void endTraversal();
};

#endif /* LEGACYCELLPROCESSOR_H_ */
