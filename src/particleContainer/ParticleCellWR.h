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

	ParticleCell_WR(const ParticleCell_WR& /*other*/):_cellDataSoA_WR(0){
		_dummy = Molecule();
	}

	~ParticleCell_WR();

	void deallocateAllParticles() override;

	bool addParticle(Molecule& particle, bool checkWhetherDuplicate = false) override;

	Molecule& moleculesAt(size_t i) override;

	const Molecule& moleculesAtConst(size_t i) const override;

	bool isEmpty() const override;

	bool deleteMoleculeByIndex(size_t index) override;

	int getMoleculeCount() const override;

	void preUpdateLeavingMolecules() override {}

	void updateLeavingMoleculesBase(ParticleCellBase& otherCell) override ;

	void postUpdateLeavingMolecules() override {}

	void getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer = false) override;

	void buildSoACaches() override {}

	void increaseMoleculeStorage(size_t numExtraMols) override;

	int countInRegion(double lowCorner[3], double highCorner[3]) const;

	void swapAndAppendToCell(ParticleCell_WR& other);

	void swapMolecules(int i, ParticleCell_WR& other, int j);

	CellDataSoA_WR & getCellDataSoA() {return _cellDataSoA_WR;}

	virtual size_t getMoleculeVectorDynamicSize() const override {return 0;}

	void prefetch() const;

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
