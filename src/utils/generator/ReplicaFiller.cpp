
#include "utils/generator/ReplicaFiller.h"

#include <string>

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


/** The ParticleContainerToBasisWrapper class is there to read in any phase space input and save it into a Basis object instead of a regular particle container, so it can be used for the grid filler.
 * @warning While providing all interface methods the ParticleContainerToBasisWrapper is not a fully working particle container!
 * @todo Find a better way to reuse all the different I/O classes and their readPhase space methods. Maybe bring back the Dummy domain decomposition?
 */
class ParticleContainerToBasisWrapper : public ParticleContainer {
public:
	ParticleContainerToBasisWrapper() {}

	~ParticleContainerToBasisWrapper() {}

	void readXML(XMLfileUnits& xmlconfig) {};

	void setBoundingBox(std::shared_ptr<Object> object) { _object = object; }

	bool addParticle(Molecule& particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false,
					 const bool& rebuildCaches = false) {
		double r[3] = {particle.r(0), particle.r(1), particle.r(2)};
		if(_object && !_object->isInside(r)) {
			return false;
		}
		_basis.addMolecule(particle);
		return true;
	}

	/** @brief return reference to internal basis object. */
	Basis& getBasis() { return _basis; }

	void clear() { _basis = Basis(); }

	unsigned long getNumberOfParticles() { return _basis.numMolecules(); }

	double getBoundingBoxMin(int dimension) const;

	bool isInBoundingBox(double r[3]) const;

	void update() {}

	void addParticles(std::vector<Molecule>& particles, bool checkWhetherDuplicate = false) {}

	void traverseCells(CellProcessor& cellProcessor) {}

	void traverseNonInnermostCells(CellProcessor& cellProcessor) {}

	void traversePartialInnermostCells(CellProcessor& cellProcessor, unsigned int stage, int stageCount) {}

	ParticleIterator iterator(ParticleIterator::Type t = ParticleIterator::ALL_CELLS) { return ParticleIterator(); }

	RegionParticleIterator regionIterator(const double startCorner[3], const double endCorner[3],
										  ParticleIterator::Type t = ParticleIterator::ALL_CELLS) { return RegionParticleIterator(); }

	double getBoundingBoxMax(int dimension) const;

	void deleteOuterParticles() {}

	double get_halo_L(int index) const { return 0.0; }

	double getCutoff() { return 0.0; }

	void deleteMolecule(Molecule& molecule, const bool& rebuildCaches) {}

	double
	getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1, CellProcessor& cellProcessor) { return 0.0; }

	void updateInnerMoleculeCaches() {}

	void updateBoundaryAndHaloMoleculeCaches() {}

	void updateMoleculeCaches() {}

	ParticleCellBase* getCell(unsigned cellIndex) { return nullptr; }

	const ParticleCellBase* getCell(unsigned cellIndex) const { return nullptr; }

	bool getMoleculeAtPosition(const double pos[3],
							   Molecule** result) { return false; } // pure virtual in particleContainer.h

	Molecule* getMoleculeCloseToPosition(const double pos[3], unsigned long id) override{};

	unsigned long initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension,
								std::array<double, 3> simBoxLength) { return 0; }

	size_t getTotalSize() { return _basis.numMolecules() * sizeof(Molecule); }

	void printSubInfo(int offset) { return; }

	std::string getName() { return std::string("ParticleContainerToBasisWrapper"); }

	double* getCellLength() override { return nullptr; }

private:
	Basis _basis;
	std::shared_ptr<Object> _object;
};


void ReplicaFiller::setObject(std::shared_ptr<Object> object) { _object = object; }

std::shared_ptr<Object> ReplicaFiller::getObject() { return _object; }

void ReplicaFiller::readXML(XMLfileUnits& xmlconfig) {
	if(xmlconfig.changecurrentnode("input")) {
		std::string inputPluginName;
		xmlconfig.getNodeValue("@type", inputPluginName);
		if(inputPluginName != "BinaryReader") {
			global_log->error() << "ReplicaFiller only works with inputPlugins: BinaryReader at the moment" << endl;
			Simulation::exit(1);
		}
// 	InputPluginFactory inputPluginFactory;
		InputBase* inputReader = new BinaryReader();
		if(inputReader == nullptr) {
			global_log->error() << "Could not create input reader " << inputPluginName << endl;
			Simulation::exit(1);
		}
		setInputReader(std::shared_ptr<InputBase>(inputReader));
		_inputReader->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->error() << "Input reader for original not specified." << endl;
		Simulation::exit(1);
	}
	if(xmlconfig.changecurrentnode("origin")) {
		Coordinate3D origin;
		origin.readXML(xmlconfig);
		origin.get(_origin);
		xmlconfig.changecurrentnode("..");
	}
	global_log->info() << "Base point for the replication: [" << _origin[0] << "," << _origin[1] << "," << _origin[2]
					   << "]" << endl;

	unsigned int componentid = 0;
	_componentid = 0;
	if(xmlconfig.getNodeValue("componentid", componentid))
		_componentid = componentid;
}


void ReplicaFiller::init() {
	ParticleContainerToBasisWrapper basisContainer;
	std::shared_ptr<Object> object = std::make_shared<ObjectShifter>(_object, _origin);
	basisContainer.setBoundingBox(object);


#ifdef ENABLE_MPI
	DomainDecomposition domainDecomp;
#else
	DomainDecompBase domainDecomp;
#endif
	Domain domain(0);
	_inputReader->readPhaseSpaceHeader(&domain, 0.0);
	_inputReader->readPhaseSpace(&basisContainer, &domain, &domainDecomp);
	global_log->info() << "Number of molecules in the replica: " << basisContainer.getNumberOfParticles() << endl;

	global_log->info() << "Setting simulation time to 0." << endl;
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

	// change component if specified
	if(molecule->componentid() != _componentid) {
		cout << "Set componentid: " << _componentid << endl;
		molecule->setComponent(global_simulation->getEnsemble()->getComponent(_componentid));
	}

	return ret;
}

