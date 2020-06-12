/**
 * @file AutoPasContainer.h
 * @author seckler
 * @date 19.09.18
 */

#pragma once

#include "ParticleContainer.h"

#include <autopas/AutoPas.h>
#include <autopas/molecularDynamics/autopasmd.h>

/**
 * A wrapper for the AutoPas library.
 */
class AutoPasContainer : public ParticleContainer {
public:
	AutoPasContainer(double cutoff);

	~AutoPasContainer() override {
#ifdef ENABLE_MPI
		_logFile.close();
#endif
	};

	// from ParticleContainer
	void readXML(XMLfileUnits &xmlconfig) override;

	bool rebuild(double bBoxMin[3], double bBoxMax[3]) override;

	void update() override;

	void forcedUpdate() override;

	bool addParticle(Molecule &particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false,
					 const bool &rebuildCaches = false) override;

	bool addHaloParticle(Molecule &particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false,
						 const bool &rebuildCaches = false) override;

	void addParticles(std::vector<Molecule> &particles, bool checkWhetherDuplicate = false) override;

	void traverseCells(CellProcessor &cellProcessor) override;

	void traverseNonInnermostCells(CellProcessor &cellProcessor) override;

	void traversePartialInnermostCells(CellProcessor &cellProcessor, unsigned int stage, int stageCount) override;

	ParticleIterator iterator(ParticleIterator::Type t) override;

	RegionParticleIterator regionIterator(const double startCorner[3], const double endCorner[3],
										  ParticleIterator::Type t) override;

	unsigned long getNumberOfParticles() override;

	void clear() override;

	void deleteOuterParticles() override;

	double get_halo_L(int index) const override;

	double getCutoff() const override;

	double getInteractionLength() const override;

	double getSkin() const override;

	void deleteMolecule(ParticleIterator &moleculeIter, const bool &rebuildCaches) override;

	double getEnergy(ParticlePairsHandler *particlePairsHandler, Molecule *m1, CellProcessor &cellProcessor) override;

	void updateInnerMoleculeCaches() override;

	void updateBoundaryAndHaloMoleculeCaches() override;

	void updateMoleculeCaches() override;

	bool getMoleculeAtPosition(const double pos[3], Molecule **result) override;

	unsigned long initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
								std::array<double, 3> simBoxLength) override;

	double *getCellLength() override;

	double *getHaloSize() override;

	// from MemoryProfilable
	size_t getTotalSize() override { return 0; }

	void printSubInfo(int offset) override {}

	std::string getName() override { return "AutoPasContainer"; }

	void setCutoff(double cutoff) override { _cutoff = cutoff; }

	std::vector<Molecule> getInvalidParticles() override {
		_hasInvalidParticles = false;
		return std::move(_invalidParticles);
	}

	bool hasInvalidParticles() override { return _hasInvalidParticles; }

	bool isInvalidParticleReturner() override { return true; }

private:
	/**
	 * Helper to get static value of shifting bool.
	 * @tparam shifting
	 */
	template <bool shifting>
	void traverseTemplateHelper();

	double _cutoff{0.};
	double _verletSkin;
	unsigned int _verletRebuildFrequency;
	unsigned int _verletClusterSize;
	unsigned int _tuningFrequency;
	unsigned int _tuningSamples;
	unsigned int _maxEvidence;
	using CellType = autopas::FullParticleCell<Molecule>;
	autopas::AutoPas<Molecule, CellType> _autopasContainer;

	std::set<autopas::TraversalOption> _traversalChoices;
	std::set<autopas::ContainerOption> _containerChoices;
	autopas::SelectorStrategyOption _selectorStrategy;
	autopas::TuningStrategyOption _tuningStrategyOption;
	autopas::AcquisitionFunctionOption _tuningAcquisitionFunction;
	std::set<autopas::DataLayoutOption> _dataLayoutChoices;
	std::set<autopas::Newton3Option> _newton3Choices;

	std::vector<Molecule> _invalidParticles;
	bool _hasInvalidParticles{false};

	ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;

#ifdef ENABLE_MPI
	std::ofstream _logFile;
#endif
};
