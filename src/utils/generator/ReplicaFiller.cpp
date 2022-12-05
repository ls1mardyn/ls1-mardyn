
#include "utils/generator/ReplicaFiller.h"

#include <string>
#include <utility>
#include <limits>
#include <vector>

#include "molecules/Molecule.h"
#include "io/InputBase.h"
// #include "io/InputPluginFactory.h"
#include "io/BinaryReader.h"
#include "utils/Logger.h"
#include "utils/Coordinate3D.h"
#include "Simulation.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/ParticleCellBase.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#else

#include "parallel/DomainDecompBase.h"

#endif

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "ensemble/CanonicalEnsemble.h"

using Log::global_log;
using std::endl;

/**
 * The ParticleContainerToBasisWrapper class is there to read in any phase space input and save it into a Basis object
 * instead of a regular particle container, so it can be used for the grid filler.
 * @warning While providing all interface methods the ParticleContainerToBasisWrapper is not a fully working particle
 * container!
 * @todo Find a better way to reuse all the different I/O classes and their readPhase space methods. Maybe bring back
 * the Dummy domain decomposition?
 */
class ParticleContainerToBasisWrapper : public ParticleContainer {
 public:
	ParticleContainerToBasisWrapper() = default;

	~ParticleContainerToBasisWrapper() override = default;

	void readXML(XMLfileUnits& xmlconfig) override {};

	void setBoundingBox(std::shared_ptr<Object> object) { _object = std::move(object); }

	bool addParticle(Molecule& particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false,
					 const bool& rebuildCaches = false) override {
		double r[3] = {particle.r(0), particle.r(1), particle.r(2)};
		if (_object && !_object->isInside(r)) {
			return false;
		}
		_basis.addMolecule(particle);
		return true;
	}

	/** @brief return reference to internal basis object. */
	Basis& getBasis() { return _basis; }

	void clear() override { _basis = Basis(); }

	unsigned long getNumberOfParticles(ParticleIterator::Type /* t */ = ParticleIterator::ALL_CELLS) override { return _basis.numMolecules(); }

	double getBoundingBoxMin(int dimension) const override {
		double min[3] = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};
		return min[dimension];
	}

	double getBoundingBoxMax(int dimension) const override {
		double max[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
		return max[dimension];
	}

	bool isInBoundingBox(double r[3]) const override {
		return true;
	}

	void update() override {}

	void addParticles(std::vector<Molecule>& particles, bool checkWhetherDuplicate = false) override {}

	void traverseCells(CellProcessor& cellProcessor) override {}

	void traverseNonInnermostCells(CellProcessor& cellProcessor) override {}

	void traversePartialInnermostCells(CellProcessor& cellProcessor, unsigned int stage, int stageCount) override {}

	ParticleIterator iterator(ParticleIterator::Type t) override { return ParticleIterator(); }

	RegionParticleIterator regionIterator(const double startCorner[3], const double endCorner[3],
										  ParticleIterator::Type t) override { return RegionParticleIterator(); }

	void deleteOuterParticles() override {}

	double get_halo_L(int index) const override { return 0.0; }

	double getCutoff() const override { return 0.0; }

	void deleteMolecule(ParticleIterator &moleculeIter, const bool& rebuildCaches) override {}

	double getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1, CellProcessor& cellProcessor) override {
		return 0.0;
	}

	void updateInnerMoleculeCaches() override {}

	void updateBoundaryAndHaloMoleculeCaches() override {}

	void updateMoleculeCaches() override {}

	// Pure virtual in ParticleContainer.h
	// Returns invalid iterator.
	std::variant<ParticleIterator, SingleCellIterator<ParticleCell>> getMoleculeAtPosition(const double pos[3]) override { return {}; }

	unsigned long initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
								std::array<double, 3> simBoxLength, size_t seed_offset) override { return 0; }

	size_t getTotalSize() override { return _basis.numMolecules() * sizeof(Molecule); }

	void printSubInfo(int offset) override { }

	std::string getName() override { return std::string("ParticleContainerToBasisWrapper"); }

	double* getCellLength() override { return nullptr; }

	string getConfigurationAsString() override {
		// give some dummy value
		return "{ParticleContainerToBasisWrapper: dummy}";
	}

 private:
	Basis _basis;
	std::shared_ptr<Object> _object;
};


void ReplicaFiller::setObject(std::shared_ptr<Object> object) { _object = object; }

std::shared_ptr<Object> ReplicaFiller::getObject() { return _object; }

void ReplicaFiller::readXML(XMLfileUnits& xmlconfig) {
	if (xmlconfig.changecurrentnode("input")) {
		std::string inputPluginName;
		xmlconfig.getNodeValue("@type", inputPluginName);
		if (inputPluginName != "BinaryReader") {
			global_log->error() << "[ReplicaFiller] ReplicaFiller only works with inputPlugins: BinaryReader at the moment" << endl;
			Simulation::exit(1);
		}
		setInputReader(std::make_shared<BinaryReader>());
		_inputReader->readXML(xmlconfig);
		if (_inputReader == nullptr) {
			global_log->error() << "[ReplicaFiller] Could not create input reader " << inputPluginName << endl;
			Simulation::exit(1);
		}
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->error() << "[ReplicaFiller] Input reader for original not specified." << endl;
		Simulation::exit(1);
	}
	if (xmlconfig.changecurrentnode("origin")) {
		Coordinate3D origin;
		origin.readXML(xmlconfig);
		origin.get(_origin);
		xmlconfig.changecurrentnode("..");
	}
	global_log->info() << "[ReplicaFiller] Base point for the replication: ["
					   << _origin[0] << "," << _origin[1] << "," << _origin[2]
					   << "]" << endl;

	unsigned int componentid = 0;
    // If the XML defines a new component ID for the replication, use it.
    // Otherwise, keep the value that is defined in the phase space file.
	if (xmlconfig.getNodeValue("componentid", componentid)) {
		const size_t numComps = global_simulation->getEnsemble()->getComponents()->size();
		if ((componentid < 1) || (componentid > numComps)) {
			global_log->error() << "[ReplicaFiller] Specified componentid is invalid. Valid range: 1 <= componentid <= " << numComps << endl;
			Simulation::exit(1);
		}
		_componentid = componentid - 1;  // Internally stored in array starting at index 0
		_keepComponent = false;
	}
	if (_keepComponent) {
		global_log->info() << "[ReplicaFiller] Keeping componentid of input" << endl;
	} else {
		global_log->info() << "[ReplicaFiller] Changing componentid of input to " << _componentid + 1 << endl;
	}
}


void ReplicaFiller::init() {
	ParticleContainerToBasisWrapper basisContainer;
	std::shared_ptr<Object> object = std::make_shared<ObjectShifter>(_object, _origin);

	// the following line is commented out, because the selection should not yet happen at this stage.
	// see https://github.com/ls1mardyn/ls1-mardyn/pull/64
	// basisContainer.setBoundingBox(object);

	Domain domain(0);
	_inputReader->readPhaseSpaceHeader(&domain, 0.0);
	_inputReader->readPhaseSpace(&basisContainer, &domain, &global_simulation->domainDecomposition());
	unsigned long numberOfParticles = basisContainer.getNumberOfParticles();
	global_log->info() << "[ReplicaFiller] Number of molecules in the replica: " << numberOfParticles << endl;

	if (numberOfParticles == 0) {
		global_log->error_always_output() << "[ReplicaFiller] No molecules in replica, aborting! " << endl;
		Simulation::exit(1);
	}

	global_log->info() << "[ReplicaFiller] Setting simulation time to 0.0" << endl;
	_simulation.setSimulationTime(0);

	Lattice lattice;
	double a[3] = {domain.getGlobalLength(0), 0.0, 0.0};
	double b[3] = {0.0, domain.getGlobalLength(1), 0.0};
	double c[3] = {0.0, 0.0, domain.getGlobalLength(2)};
	lattice.init(triclinic, primitive, a, b, c);
	_gridFiller.setObject(getObject());
	_gridFiller.init(lattice, basisContainer.getBasis(), _origin);
}

int ReplicaFiller::getMolecule(Molecule* molecule) {
	int ret = _gridFiller.getMolecule(molecule);
	if (ret != 0) {
		// change component if specified
		if (not _keepComponent && molecule->componentid() != _componentid) {
			molecule->setComponent(global_simulation->getEnsemble()->getComponent(_componentid));
		}
	}
	return ret;
}
