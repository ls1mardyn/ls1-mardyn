#ifndef SRC_PARTICLECONTAINER_PARTICLECELLRMM_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLRMM_H_

#include "ParticleCellBase.h"
#include "SingleCellIterator.h"
#include "particleContainer/adapter/CellDataSoARMM.h"

class ParticleCellRMM: public ParticleCellBase {
public:
	ParticleCellRMM();

	ParticleCellRMM(const ParticleCellRMM& other) = default;

	~ParticleCellRMM() override;

	SingleCellIterator<ParticleCellRMM> iterator() {
		return SingleCellIterator<ParticleCellRMM>(this);
	}

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

	void prefetchForForce() const override;

	void getLeavingMolecules(std::vector<Molecule> & appendBuffer) override;

	bool findMoleculeByID(size_t& index, unsigned long molid) const override;

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
