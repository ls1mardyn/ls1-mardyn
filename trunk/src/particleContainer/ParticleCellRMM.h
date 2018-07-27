#ifndef SRC_PARTICLECONTAINER_PARTICLECELLRMM_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLRMM_H_

#include "ParticleCellBase.h"
#include "particleContainer/adapter/CellDataSoARMM.h"
#include "CellBorderAndFlagManager.h"

class ParticleCellRMM: public ParticleCellBase {
public:
	ParticleCellRMM();

	ParticleCellRMM(const ParticleCellRMM& other) :
			ParticleCellBase(other), _cellDataSoARMM(other._cellDataSoARMM) {
	}

	~ParticleCellRMM();

	void deallocateAllParticles() override;

	bool addParticle(Molecule& particle, bool checkWhetherDuplicate = false) override;

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

	void swapAndAppendToCell(ParticleCellRMM& other);

	void swapMolecules(int i, ParticleCellRMM& other, int j);

	CellDataSoARMM & getCellDataSoA() {return _cellDataSoARMM;}

	size_t getMoleculeVectorDynamicSize() const override {return 0;}

	void prefetchForForce() const;

	void getLeavingMolecules(std::vector<Molecule> & appendBuffer);

	static CellBorderAndFlagManager _cellBorderAndFlagManager;
	bool isHaloCell() const {
		return _cellBorderAndFlagManager.isHaloCell(_cellIndex);
	}
	bool isBoundaryCell() const {
		return _cellBorderAndFlagManager.isBoundaryCell(_cellIndex);
	}
	bool isInnerCell() const {
		return _cellBorderAndFlagManager.isInnerCell(_cellIndex);
	}
	bool isInnerMostCell() const {
		return _cellBorderAndFlagManager.isInnerMostCell(_cellIndex);
	}

	void assignCellToHaloRegion() { mardyn_assert(isHaloCell()); }
	void assignCellToBoundaryRegion() { mardyn_assert(isBoundaryCell()); }
	void assignCellToInnerRegion() { mardyn_assert(isInnerCell()); }
	void assignCellToInnerMostAndInnerRegion() { mardyn_assert(isInnerMostCell() and isInnerCell()); }

	void skipCellFromHaloRegion() { mardyn_assert(not isHaloCell()); }
	void skipCellFromBoundaryRegion() { mardyn_assert(not isBoundaryCell()); }
	void skipCellFromInnerRegion() { mardyn_assert(not isInnerCell()); }
	void skipCellFromInnerMostRegion() { mardyn_assert(not isInnerMostCell()); }

	double getBoxMin(int d) const {
		return _cellBorderAndFlagManager.getBoundingBoxMin(_cellIndex, d);
	}
	double getBoxMax(int d) const {
		return _cellBorderAndFlagManager.getBoundingBoxMax(_cellIndex, d);
	}
	void setBoxMin(const double b[3]) {
		for (int d = 0; d < 3; ++d) {
			mardyn_assert(getBoxMin(d) == b[d]);
		}
	}
	void setBoxMax(const double b[3]) {
		for (int d = 0; d < 3; ++d) {
			mardyn_assert(getBoxMax(d) == b[d]);
		}
	}

//protected: do not use!
	void moleculesAtNew(size_t i, Molecule *& multipurposePointer) override;
	void moleculesAtConstNew(size_t i, Molecule *& multipurposePointer) const override;

private:
	Molecule buildAoSMolecule(size_t i) const;

	/**
	 * \brief Structure of arrays for VectorizedCellProcessor.
	 * \author Johannes Heckl
	 */
	CellDataSoARMM _cellDataSoARMM;
};

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLRMM_H_ */
