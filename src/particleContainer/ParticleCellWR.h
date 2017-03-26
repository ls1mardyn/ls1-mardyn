/*
 * ParticleCellWR.h
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_PARTICLECELLWR_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLWR_H_

#include "ParticleCellBase.h"
#include "particleContainer/adapter/CellDataSoA_WR.h"

class ParticleCell_WR: public ParticleCellBase {
public:
	ParticleCell_WR();
	~ParticleCell_WR();

	void deallocateAllParticles();

	bool addParticle(Molecule& particle, bool checkWhetherDuplicate = false);

	Molecule& moleculesAt(size_t i);

	bool isEmpty() const;

	bool deleteMoleculeByIndex(size_t index);

	int getMoleculeCount() const;

	void preUpdateLeavingMolecules();

	void updateLeavingMoleculesBase(ParticleCellBase& otherCell);

	void postUpdateLeavingMolecules();

	void getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer = false);

	void buildSoACaches();

	void reserveMoleculeStorage(size_t numMols);

private:
	/**
	 * \brief object used for the moleculesAt() interface
	 */
	Molecule _dummy;

	/**
	 * \brief Structure of arrays for VectorizedCellProcessor.
	 * \author Johannes Heckl
	 */
	CellDataSoA_WR _cellDataSoA_WR;
};

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLWR_H_ */
