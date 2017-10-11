/*
 * ParticleCellBase.h
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_PARTICLECELLBASE_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLBASE_H_

#include "Cell.h"
#include "SingleCellIterator.h"
#include "molecules/Molecule.h"

#ifdef QUICKSCHED
#include <quicksched.h>
#endif

class Random;

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

	SingleCellIterator iteratorBegin();
	SingleCellIterator iteratorEnd();

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

	virtual void increaseMoleculeStorage(size_t numMols) = 0;

	virtual bool testPointInCell(const double point[3]) const {
		double boxMin[3] = {getBoxMin(0), getBoxMin(1), getBoxMin(2)};
		double boxMax[3] = {getBoxMax(0), getBoxMax(1), getBoxMax(2)};
		return boxMin[0] <= point[0] && boxMin[1] <= point[1] && boxMin[2] <= point[2] &&
				point[0] < boxMax[0] && point[1] < boxMax[1] && point[2] < boxMax[2];
	}

	virtual bool testInBox(const Molecule& particle) const {
		double boxMin[3] = {getBoxMin(0), getBoxMin(1), getBoxMin(2)};
		double boxMax[3] = {getBoxMax(0), getBoxMax(1), getBoxMax(2)};
		return particle.inBox(boxMin, boxMax);
	}

	virtual size_t getMoleculeVectorDynamicSize() const = 0;

	virtual void prefetchForForce() const {/*TODO*/}

	unsigned long initCubicGrid(int numMoleculesPerDimension, double spacing, Random & RNG);

//protected: Do not use! use SingleCellIterator instead!
	// multipurpose:
	// in FullParticleCell, this is set to point to one of the molecules in _molecules
	// in ParticleCellRMM, this points to an existing molecule, into which the correct data is written.
	virtual void moleculesAtNew(size_t i, Molecule *& multipurposePointer) = 0;
	virtual void moleculesAtConstNew(size_t i, Molecule *& multipurposePointer) const = 0;

	virtual void getLeavingMolecules(std::vector<Molecule> & appendBuffer) {
		// TODO: implement for FullParticleCell
	}

#ifdef QUICKSCHED
	qsched_res_t getRescourceId() const {
		return _resourceId;
	}

	void setResourceId(qsched_res_t resourceId){
		_resourceId = resourceId;
	}

	qsched_task_t getTaskId() const {
		return _taskId;
	}

	void setTaskId(qsched_task_t taskId){
		_taskId = taskId;
	}
#endif // QUICKSCHED

protected:
	bool findMoleculeByID(size_t& index, unsigned long molid) const;

#ifdef QUICKSCHED
	qsched_res_t  _resourceId;
	qsched_task_t _taskId;
#endif // QUICKSCHED
};

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLBASE_H_ */
