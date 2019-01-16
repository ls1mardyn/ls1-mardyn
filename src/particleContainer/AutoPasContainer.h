/**
 * @file AutoPasContainer.h
 * @author seckler
 * @date 19.09.18
 */

#pragma once

#include "ParticleContainer.h"

#include <autopas/AutoPas.h>

/**
 * A wrapper for the AutoPas library.
 */
class AutoPasContainer : public ParticleContainer {
public:
	AutoPasContainer();

	~AutoPasContainer() override = default;

	// from ParticleContainer
	void readXML(XMLfileUnits &xmlconfig) override;

	bool rebuild(double bBoxMin[3], double bBoxMax[3]) override;

	void update() override;

	bool addParticle(Molecule &particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false,
					 const bool &rebuildCaches = false) override;

	bool addHaloParticle(Molecule &particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false,
						 const bool &rebuildCaches = false) override;

	void addParticles(std::vector<Molecule> &particles, bool checkWhetherDuplicate = false) override;

	void traverseCells(CellProcessor &cellProcessor) override;

	void traverseNonInnermostCells(CellProcessor &cellProcessor) override;

	void traversePartialInnermostCells(CellProcessor &cellProcessor, unsigned int stage, int stageCount) override;

	ParticleIterator iterator(ParticleIterator::Type t = ParticleIterator::ALL_CELLS) override;

	RegionParticleIterator regionIterator(const double startCorner[3], const double endCorner[3],
										  ParticleIterator::Type t = ParticleIterator::ALL_CELLS) override;

	unsigned long getNumberOfParticles() override;

	void clear() override;

	void deleteOuterParticles() override;

	double get_halo_L(int index) const override;

	double getCutoff() override;

	void deleteMolecule(Molecule &molecule, const bool &rebuildCaches) override;

	double getEnergy(ParticlePairsHandler *particlePairsHandler, Molecule *m1, CellProcessor &cellProcessor) override;

	void updateInnerMoleculeCaches() override;

	void updateBoundaryAndHaloMoleculeCaches() override;

	void updateMoleculeCaches() override;

	bool getMoleculeAtPosition(const double pos[3], Molecule **result) override;

	unsigned long initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
								std::array<double, 3> simBoxLength) override;

	double *getCellLength() override;

	// from MemoryProfilable
	size_t getTotalSize() override { return 0; }

	void printSubInfo(int offset) override {}

	std::string getName() override { return "AutoPasContainer"; }

	void setCutoff(double cutoff) override { _cutoff = cutoff; }

private:
	double _cutoff;
	double _verletSkin;
	unsigned int _verletRebuildFrequency;
	unsigned int _tuningFrequency;
	unsigned int _tuningSamples;
	typedef autopas::FullParticleCell<Molecule> CellType;
	autopas::AutoPas<Molecule, CellType> _autopasContainer;

	std::vector<autopas::TraversalOptions> _traversalChoices;
	std::vector<autopas::ContainerOptions> _containerChoices;
	autopas::SelectorStrategy _traversalSelectorStrategy;
	autopas::SelectorStrategy _containerSelectorStrategy;
	autopas::DataLayoutOption _dataLayout;
};
