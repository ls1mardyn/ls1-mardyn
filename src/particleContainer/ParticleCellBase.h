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

	/**
	 * \param particle the particle to be added
	 * \param checkWhetherDuplicate if true, perform a check by molecule IDs,
	 * whether a particle with the same ID already exists
	 * \return true, if inserted
	 */
	virtual bool addParticle(Molecule& particle, bool checkWhetherDuplicate = false) = 0;

	virtual Molecule& moleculesAt(size_t i) = 0;

	virtual const Molecule& moleculesAtConst(size_t i) const = 0;

	virtual bool isEmpty() const = 0;

	bool isNotEmpty() const {return not isEmpty();}

	bool deleteMoleculeByID(unsigned long molid);

	virtual bool deleteMoleculeByIndex(size_t index) = 0;

	virtual int getMoleculeCount() const = 0;

	virtual void preUpdateLeavingMolecules() = 0;

	virtual void updateLeavingMoleculesBase(ParticleCellBase& otherCell) = 0;

	virtual void postUpdateLeavingMolecules() = 0;

	virtual void getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer = false) = 0;

	virtual void buildSoACaches() = 0;

	virtual void reserveMoleculeStorage(size_t numMols) = 0;

	bool testPointInCell(const double point[3]) const {
		return _boxMin[0] <= point[0] && _boxMin[1] <= point[1] && _boxMin[2] <= point[2] &&
				point[0] < _boxMax[0] && point[1] < _boxMax[1] && point[2] < _boxMax[2];
	}

	bool testInBox(const Molecule& particle) const {
		return particle.inBox(_boxMin, _boxMax);
	}

protected:
	void findMoleculeByID(bool& wasFound, size_t& index, unsigned long molid) const;
};

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLBASE_H_ */
