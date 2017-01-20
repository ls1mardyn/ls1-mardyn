/*
 * ParticleCellBase.h
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_PARTICLECELLBASE_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLBASE_H_

#include "Cell.h"
#include "molecules/Molecule.h"

class ParticleCellBase: public Cell {
public:
	ParticleCellBase();
	virtual ~ParticleCellBase();

	virtual void deallocateAllParticles() = 0;

	virtual bool addParticle(Molecule& particle, bool checkWhetherDuplicate = false) = 0;

	virtual Molecule& moleculesAt(size_t i) = 0;

	virtual bool isEmpty() const = 0;

	bool deleteMoleculeByID(unsigned long molid);

	virtual bool deleteMoleculeByIndex(size_t index) = 0;

	virtual int getMoleculeCount() const = 0;

	virtual void preUpdateLeavingMolecules() = 0;

	virtual void updateLeavingMoleculesBase(ParticleCellBase& otherCell) = 0;

	virtual void postUpdateLeavingMolecules() = 0;

	virtual void getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer = false) = 0;

    virtual void buildSoACaches() = 0;

	bool testInBox(const Molecule& particle) const {
		return particle.inBox(_boxMin, _boxMax);
	}
};

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLBASE_H_ */
