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

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2);

	void processCell(ParticleCell& cell);

	void postprocessCell(ParticleCell& /*cell*/) {}

	void endTraversal();

protected:
	/**
	 * Implementation of processCellPair that only sums the macroscopic values of Molecule-Molecule pairs,
	 * not of Molecule-Halo pairs and does not calculate Halo-Halo pairs.
	 */
	virtual void processCellPairSumHalf(ParticleCell& cell1, ParticleCell& cell2) override;
	/**
	 * Implementation of processCellPair that sums all macroscopic values and calculates Halo-Halo pairs.
	 */
	virtual void processCellPairSumAll(ParticleCell& cell1, ParticleCell& cell2) override;
};

#endif /* LEGACYCELLPROCESSOR_H_ */
