#include <cmath>
#include "particleContainer/LinkedCells.h"


#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "ParticleCell.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"
#include "utils/Random.h"
#include <array>
#include <algorithm>

#include "particleContainer/TraversalTuner.h"

// traversaal for each zonal method
#include "particleContainer/LinkedCellTraversals/C08CellPairTraversal.h"
#include "particleContainer/LinkedCellTraversals/OriginalCellPairTraversal.h"
#include "particleContainer/LinkedCellTraversals/HalfShellTraversal.h"
#include "particleContainer/LinkedCellTraversals/QuickschedTraversal.h"
#include "particleContainer/LinkedCellTraversals/SlicedCellPairTraversal.h"

#include "ResortCellProcessorSliced.h"


using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################

LinkedCells::LinkedCells() : ParticleContainer(), _traversalTuner(new TraversalTuner<ParticleCell>()), _resortCellProcessorSliced(nullptr) {
}

LinkedCells::LinkedCells(double bBoxMin[3], double bBoxMax[3],
		double cutoffRadius) :
		ParticleContainer(bBoxMin, bBoxMax), _traversalTuner(new TraversalTuner<ParticleCell>()) {
	int numberOfCells = 1;
	_cutoffRadius = cutoffRadius;

	global_log->debug() << "cutoff: " << cutoffRadius << endl;
	global_log->debug() << "# cells in cutoff hardcoded to 1 " << endl;

	for (int d = 0; d < 3; d++) {
		/* first calculate the cell length for this dimension */
		_boxWidthInNumCells[d] = floor(
				(_boundingBoxMax[d] - _boundingBoxMin[d]) / cutoffRadius);
		// in each dimension at least one layer of (inner+boundary) cells is necessary
		if (_boxWidthInNumCells[d] == 0) {
			_boxWidthInNumCells[d] = 1;
		}

		double diff = _boundingBoxMax[d] - _boundingBoxMin[d];
		_cellLength[d] = diff / _boxWidthInNumCells[d];
		_cellLengthReciprocal[d] = _boxWidthInNumCells[d] / diff;

		_haloWidthInNumCells[d] = 1;
		_haloLength[d] = _haloWidthInNumCells[d] * _cellLength[d];
		_haloBoundingBoxMin[d] = _boundingBoxMin[d] - _haloLength[d];
		_haloBoundingBoxMax[d] = _boundingBoxMax[d] + _haloLength[d];

		_cellsPerDimension[d] = _boxWidthInNumCells[d]
				+ 2 * _haloWidthInNumCells[d];

		numberOfCells *= _cellsPerDimension[d];
		mardyn_assert(numberOfCells > 0);
	}
	global_log->debug() << "Cell size (" << _cellLength[0] << ", "
			<< _cellLength[1] << ", " << _cellLength[2] << ")" << endl;
	global_log->info() << "Cells per dimension (incl. halo): "
			<< _cellsPerDimension[0] << " x " << _cellsPerDimension[1] << " x "
			<< _cellsPerDimension[2] << endl;

	_cells.resize(numberOfCells);

	// If the width of the inner region is less than the width of the halo
	// region a parallelization is not possible (with the used algorithms).
	// If a particle leaves this box, it would need to be communicated to the two next neighbors.
	if (_boxWidthInNumCells[0] < 2 * _haloWidthInNumCells[0]
			|| _boxWidthInNumCells[1] < 2 * _haloWidthInNumCells[1]
			|| _boxWidthInNumCells[2] < 2 * _haloWidthInNumCells[2]) {
		global_log->error_always_output()
				<< "LinkedCells (constructor): bounding box too small for calculated cell length"
				<< endl;
		global_log->error_always_output() << "_cellsPerDimension: " << _cellsPerDimension[0]
				<< " / " << _cellsPerDimension[1] << " / "
				<< _cellsPerDimension[2] << endl;
		global_log->error_always_output() << "_haloWidthInNumCells: "
				<< _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1]
				<< " / " << _haloWidthInNumCells[2] << endl;
		global_log->error_always_output() << "_boxWidthInNumCells: " << _boxWidthInNumCells[0]
				<< " / " << _boxWidthInNumCells[1] << " / "
				<< _boxWidthInNumCells[2] << endl;
		Simulation::exit(5);
	}

	initializeCells();
	initializeTraversal();

	_cellsValid = false;
	_resortCellProcessorSliced = nullptr;
}

LinkedCells::~LinkedCells() {
	delete _resortCellProcessorSliced;

	std::vector<ParticleCell>::iterator it;
	for (it = _cells.begin(); it != _cells.end(); ++it) {
		it->deallocateAllParticles();
	}
}

void LinkedCells::initializeTraversal() {
	std::array<long unsigned, 3> dims;
	for (int d = 0; d < 3; ++d) {
		dims[d] = _cellsPerDimension[d];
	}
	_traversalTuner->rebuild(_cells, dims);
}

void LinkedCells::readXML(XMLfileUnits& xmlconfig) {
	_cellsInCutoff = xmlconfig.getNodeValue_int("cellsInCutoffRadius", 1); // new
	mardyn_assert(_cellsInCutoff>=1); // new

	_traversalTuner = std::unique_ptr<TraversalTuner<ParticleCell>>(new TraversalTuner<ParticleCell>()); // new way to assign _traversalTuner
	_traversalTuner->readXML(xmlconfig);
}

bool LinkedCells::rebuild(double bBoxMin[3], double bBoxMax[3]) {
	global_log->info() << "REBUILD OF LinkedCells" << endl;

	for (int i = 0; i < 3; i++) {
		this->_boundingBoxMin[i] = bBoxMin[i];
		this->_boundingBoxMax[i] = bBoxMax[i];
//		_haloWidthInNumCells[i] = ::ceil(_cellsInCutoff);
		_haloWidthInNumCells[i] = _cellsInCutoff;
	}
	global_log->info() << "Bounding box: " << "[" << bBoxMin[0] << ", " << bBoxMax[0] << "]" << " x " << "["
			<< bBoxMin[1] << ", " << bBoxMax[1] << "]" << " x " << "[" << bBoxMin[2] << ", " << bBoxMax[2] << "]"
			<< std::endl;

	int numberOfCells = 1;

	global_log->info() << "Using " << _cellsInCutoff << " cells in cutoff." << endl;
	float rc = (_cutoffRadius / _cellsInCutoff);

	for (int dim = 0; dim < 3; dim++) {
		_boxWidthInNumCells[dim] = floor((_boundingBoxMax[dim] - _boundingBoxMin[dim]) / rc);

		_cellsPerDimension[dim] = _boxWidthInNumCells[dim] + 2 * _haloWidthInNumCells[dim];

		// in each dimension at least one layer of (inner+boundary) cells necessary
		if (_cellsPerDimension[dim] == 2 * _haloWidthInNumCells[dim]) {
			global_log->error_always_output() << "LinkedCells::rebuild: region too small" << endl;
			Simulation::exit(1);
		}

		numberOfCells *= _cellsPerDimension[dim];

		double diff = _boundingBoxMax[dim] - _boundingBoxMin[dim];
		_cellLength[dim] = diff / _boxWidthInNumCells[dim];
		_cellLengthReciprocal[dim] = _boxWidthInNumCells[dim] / diff;

		_haloLength[dim] = _haloWidthInNumCells[dim] * _cellLength[dim];

		_haloBoundingBoxMin[dim] = _boundingBoxMin[dim] - _haloLength[dim];
		_haloBoundingBoxMax[dim] = _boundingBoxMax[dim] + _haloLength[dim];
	}

	global_log->info() << "Cells per dimension (incl. halo): " << _cellsPerDimension[0] << " x "
			<< _cellsPerDimension[1] << " x " << _cellsPerDimension[2] << endl;


	_cells.resize(numberOfCells);

	bool sendParticlesTogether = true;
	// If the with of the inner region is less than the width of the halo region
	// leaving particles and halo copy must be sent separately.
	if (_boxWidthInNumCells[0] < 2 * _haloWidthInNumCells[0]
			|| _boxWidthInNumCells[1] < 2 * _haloWidthInNumCells[1]
			|| _boxWidthInNumCells[2] < 2 * _haloWidthInNumCells[2]) {
		sendParticlesTogether = false;
	}

	initializeCells();

	// TODO: We loose particles here as they are not communicated to the new owner
	// delete all Particles which are outside of the halo region
	deleteParticlesOutsideBox(_haloBoundingBoxMin, _haloBoundingBoxMax);

	initializeTraversal();

	_cellsValid = false;

	return sendParticlesTogether;

}

void LinkedCells::check_molecules_in_box() {
	std::vector<Molecule> badMolecules;
	unsigned numBadMolecules = 0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+ : numBadMolecules)
	#endif
	{
		for (ParticleIterator tM = iterator(); tM.isValid(); ++tM) {
			if (not tM->inBox(_haloBoundingBoxMin, _haloBoundingBoxMax)) {
				numBadMolecules++;

				#if defined(_OPENMP)
				#pragma omp critical
				#endif
				{
					badMolecules.push_back(*tM);
				}

			}
		}
	}

	if (numBadMolecules > 0) {
		global_log->error() << "Found " << numBadMolecules << " outside of bounding box:" << std::endl;
		for (auto & m : badMolecules) {
			global_log->error() << "Particle (id=" << m.getID() << "), (current position: x="
					<< m.r(0) << ", y=" << m.r(1) << ", z=" << m.r(2) << ")" << std::endl;
		}
		global_log->error() << "The bounding box is: [" << _haloBoundingBoxMin[0] << ", " << _haloBoundingBoxMax[0]
				<< ") x [" << _haloBoundingBoxMin[1] << ", " << _haloBoundingBoxMax[1] << ") x [" << _haloBoundingBoxMin[2]
				<< ", " << _haloBoundingBoxMax[2] << ")" << std::endl;
		global_log->error() << "Particles will be lost. Aborting simulation." << std::endl;
		Simulation::exit(311);
	}
}

void LinkedCells::update() {
#ifndef NDEBUG
	check_molecules_in_box();
#endif

	// TODO: replace via a cellProcessor and a traverseCells call ?
#ifndef ENABLE_REDUCED_MEMORY_MODE
	update_via_copies();
#else
//	update_via_coloring();
	std::array<long unsigned, 3> dims = {
		static_cast<long unsigned>(_cellsPerDimension[0]),
		static_cast<long unsigned>(_cellsPerDimension[1]),
		static_cast<long unsigned>(_cellsPerDimension[2])
	};
	if (_traversalTuner->isTraversalApplicable(TraversalTuner<ParticleCell>::traversalNames::SLICED, dims)) {
		if (_resortCellProcessorSliced == nullptr) {
			_resortCellProcessorSliced = new ResortCellProcessorSliced(this);
		}
		_traversalTuner->traverseCellPairs(TraversalTuner<ParticleCell>::traversalNames::SLICED, *_resortCellProcessorSliced);
	} else {
		update_via_traversal();
	}
#endif

	_cellsValid = true;



#ifndef NDEBUG
	unsigned numBadMolecules = 0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+: numBadMolecules)
	#endif
	{
		for (ParticleIterator tM = iterator(); tM.isValid(); ++tM) {
			if (not _cells[tM.getCellIndex()].testInBox(*tM)) {
				numBadMolecules++;
				global_log->error_always_output() << "particle " << tM->getID() << " in cell " << tM.getCellIndex()
						<< ", which is" << (_cells[tM.getCellIndex()].isBoundaryCell() ? "" : " NOT")
						<< " a boundarycell is outside of its cell after LinkedCells::update()." << std::endl;
				global_log->error_always_output() << "particle at (" << tM->r(0) << ", " << tM->r(1) << ", " << tM->r(2) << ")"
						<< std::endl << "cell: [" << _cells[tM.getCellIndex()].getBoxMin(0) << ", "
						<< _cells[tM.getCellIndex()].getBoxMax(0) << "] x [" << _cells[tM.getCellIndex()].getBoxMin(1)
						<< ", " << _cells[tM.getCellIndex()].getBoxMax(1) << "] x ["
						<< _cells[tM.getCellIndex()].getBoxMin(2) << ", " << _cells[tM.getCellIndex()].getBoxMax(2) << "]"
						<< std::endl;

			}
		}
	}


	if (numBadMolecules > 0) {
		global_log->error() << "Found " << numBadMolecules << " outside of their correct cells. Aborting." << std::endl;
		Simulation::exit(311);
	}
#endif
}

void LinkedCells::update_via_copies() {
	const vector<ParticleCell>::size_type numCells = _cells.size();
	vector<long> forwardNeighbourOffsets; // now vector
	vector<long> backwardNeighbourOffsets; // now vector
	calculateNeighbourIndices(forwardNeighbourOffsets, backwardNeighbourOffsets);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		#if defined(_OPENMP)
		#pragma omp for schedule(static)
		#endif
		for (vector<ParticleCell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++) {
			_cells[cellIndex].preUpdateLeavingMolecules();
		}

		#if defined(_OPENMP)
		#pragma omp for schedule(static)
		#endif
		for (vector<ParticleCell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++) {
			ParticleCell& cell = _cells[cellIndex];

			for (unsigned long j = 0; j < backwardNeighbourOffsets.size(); j++) {
				const unsigned long neighbourIndex = cellIndex - backwardNeighbourOffsets.at(j); // now vector
				if (neighbourIndex >= _cells.size()) {
					// handles cell_index < 0 (indices are unsigned!)
					mardyn_assert(cell.isHaloCell());
					continue;
				}
				cell.updateLeavingMoleculesBase(_cells[neighbourIndex]);
			}

			for (unsigned long j = 0; j < forwardNeighbourOffsets.size(); j++) {
				const unsigned long neighbourIndex = cellIndex + forwardNeighbourOffsets.at(j); // now vector
				if (neighbourIndex >= numCells) {
					mardyn_assert(cell.isHaloCell());
					continue;
				}
				cell.updateLeavingMoleculesBase(_cells[neighbourIndex]);
			}
		}

		#if defined(_OPENMP)
		#pragma omp for schedule(static)
		#endif
		for (vector<ParticleCell>::size_type cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
			_cells[cellIndex].postUpdateLeavingMolecules();
		}
	} // end pragma omp parallel
}
void LinkedCells::update_via_coloring() {
	std::array<std::pair<unsigned long, unsigned long>, 14> cellPairOffsets = calculateCellPairOffsets();

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int strides[3] = {2, 2, 2};
		for (int col = 0; col < 8; ++col) {
			int startIndices[3];
			threeDIndexOfCellIndex(col, startIndices, strides);

			#if defined (_OPENMP)
			#pragma omp for schedule(dynamic, 1) collapse(3)
			#endif
			for (int z = startIndices[2]; z < _cellsPerDimension[2]-1 ; z+= strides[2]) {
				for (int y = startIndices[1]; y < _cellsPerDimension[1]-1; y += strides[1]) {
					for (int x = startIndices[0]; x < _cellsPerDimension[0]-1; x += strides[0]) {
						long int baseIndex = cellIndexOf3DIndex(x, y, z);

						const int num_pairs = cellPairOffsets.size();
						for(int j = 0; j < num_pairs; ++j) {
							pair<long int, long int> current_pair = cellPairOffsets[j];

							long int offset1 = current_pair.first;
							long int cellIndex1 = baseIndex + offset1;
							if ((cellIndex1 < 0) || (cellIndex1 >= (int) (_cells.size())))
								continue;

							long int offset2 = current_pair.second;
							long int cellIndex2 = baseIndex + offset2;
							if ((cellIndex2 < 0) || (cellIndex2 >= (int) (_cells.size())))
								continue;

							ParticleCell& cell1 = _cells[cellIndex1];
							ParticleCell& cell2 = _cells[cellIndex2];

							if(cell1.isHaloCell() and cell2.isHaloCell()) {
								continue;
							}

							if(cellIndex1 == cellIndex2) {
								continue;
							}
							else {
								cell1.updateLeavingMoleculesBase(cell2);
							}
						}
					}
				}
			}
		}
	} // end pragma omp parallel
}

void LinkedCells::update_via_sliced_traversal() {
	if (_resortCellProcessorSliced == nullptr) {
		_resortCellProcessorSliced = new ResortCellProcessorSliced(this);
	}
	_traversalTuner->traverseCellPairs(*_resortCellProcessorSliced);
}

void LinkedCells::update_via_traversal() {
	class ResortCellProcessor : public CellProcessor {
	public:
		ResortCellProcessor() : CellProcessor(0.0, 0.0) {}
		void initTraversal() {}
		void preprocessCell(ParticleCell& ) {}

		void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) { // does this need a bool?
				cell1.updateLeavingMoleculesBase(cell2);
		}

		void processCell(ParticleCell& cell) {}
		double processSingleMolecule(Molecule*, ParticleCell& ) { return 0.0;}
		void postprocessCell(ParticleCell& ) {}
		void endTraversal() {}

	} resortCellProcessor;
	_traversalTuner->traverseCellPairs(resortCellProcessor);
}

bool LinkedCells::addParticle(Molecule& particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate, const bool& rebuildCaches) {
	bool wasInserted = false;
	const bool inBox = inBoxCheckedAlready or particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax);
	if (inBox) {
		int cellIndex = getCellIndexOfMolecule(&particle);
		wasInserted = _cells[cellIndex].addParticle(particle, checkWhetherDuplicate);
		if (rebuildCaches) {
			_cells[cellIndex].buildSoACaches();
		}
	}
	return wasInserted;
}

void LinkedCells::addParticles(vector<Molecule>& particles, bool checkWhetherDuplicate) {
	typedef vector<Molecule>::size_type mol_index_t;
	typedef vector<ParticleCell>::size_type cell_index_t;

#ifndef NDEBUG
	int oldNumberOfParticles = getNumberOfParticles();
#endif

	const mol_index_t N = particles.size();

	map<cell_index_t, vector<mol_index_t>> newPartsPerCell;

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		map<cell_index_t, vector<mol_index_t>> local_newPartsPerCell;

		#if defined(_OPENMP)
		#pragma omp for schedule(static)
		#endif
		for (mol_index_t i = 0; i < N; ++i) {
			Molecule & particle = particles[i];

			#ifndef NDEBUG
				if(!particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax)){
					global_log->error()<<"At particle with ID "<< particle.getID()<<" assertion failed..."<<endl;
				}
				mardyn_assert(particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax));
			#endif

			const unsigned long cellIndex = getCellIndexOfMolecule(&particle);
			mardyn_assert(cellIndex < _cells.size());
			local_newPartsPerCell[cellIndex].push_back(i);
		}

		#if defined(_OPENMP)
		#pragma omp critical(add_particles_reduce_maps)
		#endif
		{
			for (auto it = local_newPartsPerCell.begin(); it != local_newPartsPerCell.end(); ++it) {
				cell_index_t cellIndex = it->first;
				vector<mol_index_t> & global_vector = newPartsPerCell[cellIndex];
				vector<mol_index_t> & local_vector = it->second;
				global_vector.insert(global_vector.end(), local_vector.begin(), local_vector.end());
			}
		}

		#if defined(_OPENMP)
		#pragma omp barrier
		#endif

		const cell_index_t numCells = newPartsPerCell.size();
		const int thread_id = mardyn_get_thread_num();
		const int num_threads = mardyn_get_num_threads();
		const cell_index_t my_cells_start = (thread_id    ) * numCells / num_threads;
		const cell_index_t my_cells_end   = (thread_id + 1) * numCells / num_threads;

		map<cell_index_t, vector<mol_index_t>>::const_iterator thread_begin = newPartsPerCell.begin();
		advance(thread_begin, my_cells_start);
		map<cell_index_t, vector<mol_index_t>>::const_iterator thread_end = newPartsPerCell.begin();
		advance(thread_end, my_cells_end);

		for(auto it = thread_begin; it != thread_end; ++it) {
			cell_index_t cellIndex = it->first;
			const vector<mol_index_t> & global_vector = it->second;

			const size_t numMolsInCell = global_vector.size();
			_cells[cellIndex].increaseMoleculeStorage(numMolsInCell);

			for(size_t j = 0; j < numMolsInCell; ++j) {
				const mol_index_t molIndex = global_vector[j];
				Molecule & mol = particles[molIndex];
				_cells[cellIndex].addParticle(mol, checkWhetherDuplicate);
			}
		}

	} // end pragma omp parallel

#ifndef NDEBUG
	int numberOfAddedParticles = getNumberOfParticles() - oldNumberOfParticles;
	global_log->debug()<<"In LinkedCells::addParticles :"<<endl;
	global_log->debug()<<"\t#Particles to be added = "<<particles.size()<<endl;
	global_log->debug()<<"\t#Particles actually added = "<<numberOfAddedParticles<<endl;
#endif

	return;
}

void LinkedCells::traverseNonInnermostCells(CellProcessor& cellProcessor) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (traverseNonInnermostCells) invalid, call update first" << endl;
		Simulation::exit(1);
	}

	_traversalTuner->traverseCellPairsOuter(cellProcessor);
}

void LinkedCells::traversePartialInnermostCells(CellProcessor& cellProcessor, unsigned int stage, int stageCount) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (traversePartialInnermostCells) invalid, call update first" << endl;
		Simulation::exit(1);
	}

	_traversalTuner->traverseCellPairsInner(cellProcessor, stage, stageCount);
}

void LinkedCells::traverseCells(CellProcessor& cellProcessor) {
	if (_cellsValid == false) {
		global_log->error()
				<< "Cell structure in LinkedCells (traversePairs) invalid, call update first"
				<< endl;
		Simulation::exit(1);
	}

	cellProcessor.initTraversal();
	_traversalTuner->traverseCellPairs(cellProcessor);
	cellProcessor.endTraversal();
}

unsigned long LinkedCells::getNumberOfParticles() {
	unsigned long N = 0;
	unsigned long numCells = _cells.size();

	#if defined(_OPENMP)
	#pragma omp parallel for reduction(+:N)
	#endif
	for (unsigned long i = 0; i < numCells; ++i) {
		N += _cells.at(i).getMoleculeCount();
	}
	return N;
}

void LinkedCells::clear() {
	vector<ParticleCell>::iterator cellIter;
	for (cellIter = _cells.begin(); cellIter != _cells.end(); cellIter++) {
		cellIter->deallocateAllParticles();
	}
}

void LinkedCells::deleteParticlesOutsideBox(double boxMin[3], double boxMax[3]) {
	// This should be unimportant

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for (auto it = iterator(); it.isValid(); ++it) {
		bool outside = not it->inBox(boxMin, boxMax);
		if (outside) {
			it.deleteCurrentParticle();
		}
	}
}

void LinkedCells::deleteOuterParticles() {
	/*if (_cellsValid == false) {
		global_log->error()
				<< "Cell structure in LinkedCells (deleteOuterParticles) invalid, call update first"
				<< endl;
		Simulation::exit(1);
	}*/

	const int numHaloCells = _haloCellIndices.size();

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < numHaloCells; i++) {
		ParticleCell& currentCell = _cells[_haloCellIndices[i]];
		currentCell.deallocateAllParticles();
	}
}

double LinkedCells::get_halo_L(int index) const {
	return _haloLength[index];
}

RegionParticleIterator LinkedCells::regionIterator(const double startRegion[3], const double endRegion[3], ParticleIterator::Type type) {
	// parameter "type" not yet used
	// add functionality in a future version...
	unsigned int startRegionCellIndex;
	unsigned int endRegionCellIndex;

	getCellIndicesOfRegion(startRegion, endRegion, startRegionCellIndex, endRegionCellIndex);

	return getRegionParticleIterator(startRegion, endRegion, startRegionCellIndex, endRegionCellIndex, type);
}

//################################################
//############ PRIVATE METHODS ###################
//################################################

void LinkedCells::initializeCells() {
	//TODO: resize _haloCellIndices
	_haloCellIndices.clear();

	long int cellIndex;

	ParticleCell::_cellBorderAndFlagManager.init(_cellsPerDimension,
			_haloBoundingBoxMin, _haloBoundingBoxMax,
			_boundingBoxMin, _boundingBoxMax,
			_cellLength);

	for (int iz = 0; iz < _cellsPerDimension[2]; ++iz) {
		for (int iy = 0; iy < _cellsPerDimension[1]; ++iy) {
			for (int ix = 0; ix < _cellsPerDimension[0]; ++ix) {
				cellIndex = cellIndexOf3DIndex(ix, iy, iz);
				ParticleCell & cell = _cells[cellIndex];
				cell.setCellIndex(cellIndex); //set the index of the cell to the index of it...

				if (ix < _haloWidthInNumCells[0] ||
					iy < _haloWidthInNumCells[1] ||
					iz < _haloWidthInNumCells[2] ||
					ix >= _cellsPerDimension[0] - _haloWidthInNumCells[0] ||
					iy >= _cellsPerDimension[1] - _haloWidthInNumCells[1] ||
					iz >= _cellsPerDimension[2] - _haloWidthInNumCells[2]) {

					cell.assignCellToHaloRegion();
					_haloCellIndices.push_back(cellIndex);
				}
			}
		}
	}
}

void LinkedCells::calculateNeighbourIndices(std::vector<long>& forwardNeighbourOffsets, std::vector<long>& backwardNeighbourOffsets) const {
	global_log->debug() << "Setting up cell neighbour indice lists." << endl;

	// 13 neighbors for _haloWidthInNumCells = 1 or 64 for =2
	int nNeighbours = ( (2*_haloWidthInNumCells[0]+1) * (2*_haloWidthInNumCells[1]+1) * (2*_haloWidthInNumCells[2]+1) - 1) / 2;

	// Resize offset vector to number of neighbors and fill with 0
	forwardNeighbourOffsets.resize(nNeighbours, 0);
	backwardNeighbourOffsets.resize(nNeighbours, 0);

	int forwardNeighbourIndex = 0, backwardNeighbourIndex = 0;

	double xDistanceSquare;
	double yDistanceSquare;
	double zDistanceSquare;
	double cutoffRadiusSquare = pow(_cutoffRadius, 2);
	for (int zIndex = -_haloWidthInNumCells[2];
			zIndex <= _haloWidthInNumCells[2]; zIndex++) {
		// The distance in one dimension is the width of a cell multiplied with the number
		// of cells between the two cells (this is received by subtracting one of the
		// absolute difference of the cells, if this difference is not zero)
		if (zIndex != 0) {
			zDistanceSquare = pow((abs(zIndex) - 1) * _cellLength[2], 2);
		} else {
			zDistanceSquare = 0;
		}
		for (int yIndex = -_haloWidthInNumCells[1];
				yIndex <= _haloWidthInNumCells[1]; yIndex++) {
			if (yIndex != 0) {
				yDistanceSquare = pow((abs(yIndex) - 1) * _cellLength[1], 2);
			} else {
				yDistanceSquare = 0;
			}
			for (int xIndex = -_haloWidthInNumCells[0];
					xIndex <= _haloWidthInNumCells[0]; xIndex++) {
				if (xIndex != 0) {
					xDistanceSquare = pow((abs(xIndex) - 1) * _cellLength[0], 2);
				} else {
					xDistanceSquare = 0;
				}
				if (xDistanceSquare + yDistanceSquare + zDistanceSquare
						<= cutoffRadiusSquare) {
					long int offset = cellIndexOf3DIndex(xIndex, yIndex,
							zIndex);
					if (offset > 0) {
						forwardNeighbourOffsets.at(forwardNeighbourIndex) = offset; // now vector
						++forwardNeighbourIndex;
					}
					if (offset < 0) {
						backwardNeighbourOffsets.at(backwardNeighbourIndex) = abs(offset); // now vector
						++backwardNeighbourIndex;
					}
				}
			}
		}
	}

	mardyn_assert(forwardNeighbourIndex == nNeighbours);
	mardyn_assert(backwardNeighbourIndex == nNeighbours);
}
std::array<std::pair<unsigned long, unsigned long>, 14> LinkedCells::calculateCellPairOffsets() const {
	long int o   = cellIndexOf3DIndex(0,0,0); // origin
	long int x   = cellIndexOf3DIndex(1,0,0); // displacement to the right
	long int y   = cellIndexOf3DIndex(0,1,0); // displacement ...
	long int z   = cellIndexOf3DIndex(0,0,1);
	long int xy  = cellIndexOf3DIndex(1,1,0);
	long int yz  = cellIndexOf3DIndex(0,1,1);
	long int xz  = cellIndexOf3DIndex(1,0,1);
	long int xyz = cellIndexOf3DIndex(1,1,1);

	// minimize number of cells simultaneously in memory:
	std::array<std::pair<unsigned long, unsigned long>, 14> cellPairOffsets;

	cellPairOffsets[ 0] = make_pair(o, xyz);
	// evict xyz

	cellPairOffsets[ 1] = make_pair(o, yz );
	cellPairOffsets[ 2] = make_pair(x, yz );
	// evict yz

	cellPairOffsets[ 3] = make_pair(o, x  );

	cellPairOffsets[ 4] = make_pair(o, xy );
	cellPairOffsets[ 5] = make_pair(xy, z );
	// evict xy

	cellPairOffsets[ 6] = make_pair(o, z  );
	cellPairOffsets[ 7] = make_pair(x, z  );
	cellPairOffsets[ 8] = make_pair(y, z  );
	// evict z

	cellPairOffsets[ 9] = make_pair(o, y  );
	cellPairOffsets[10] = make_pair(x, y  );
	// evict x

	cellPairOffsets[11] = make_pair(o, xz );
	cellPairOffsets[12] = make_pair(y, xz );
	// evict xz

	cellPairOffsets[13] = make_pair(o, o  );

	return cellPairOffsets;
}

unsigned long int LinkedCells::getCellIndexOfMolecule(Molecule* molecule) const {
	double r[3] = {molecule->r(0), molecule->r(1), molecule->r(2)};
	return getCellIndexOfPoint(r);

	/*
	for (int dim = 0; dim < 3; dim++) {
		#ifndef NDEBUG
		if (molecule->r(dim) < _haloBoundingBoxMin[dim] || molecule->r(dim) >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "Molecule is outside of bounding box" << endl;
			global_log->error() << "Molecule:\n" << *molecule << endl;
			global_log->error() << "_haloBoundingBoxMin = (" << _haloBoundingBoxMin[0] << ", " << _haloBoundingBoxMin[1] << ", " << _haloBoundingBoxMin[2] << ")" << endl;
			global_log->error() << "_haloBoundingBoxMax = (" << _haloBoundingBoxMax[0] << ", " << _haloBoundingBoxMax[1] << ", " << _haloBoundingBoxMax[2] << ")" << endl;
			Simulation::exit(1);
		}
		#endif
		//this version is sensitive to roundoffs, if we have molecules (initialized) precisely at position 0.0:
		//cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);
		cellIndex[dim] = min(max(((int) floor((molecule->r(dim) - _boundingBoxMin[dim]) / _cellLength[dim])) + _haloWidthInNumCells[dim],0),_cellsPerDimension[dim]-1);

	}
	int cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
	// in very rare cases rounding is stupid, thus we need a check...
	//TODO: check if this can in any way be done better...
	if (_cells[cellIndex1d].testInBox(*molecule)) {
		return cellIndex1d;
	} else {
		for (int dim = 0; dim < 3; dim++) {
			if(molecule->r(dim)<_cells[cellIndex1d].getBoxMin(dim)){
				cellIndex[dim]--;
			}else if(molecule->r(dim)>=_cells[cellIndex1d].getBoxMax(dim)){
				cellIndex[dim]++;
			}
		}
		cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
		mardyn_assert(_cells[cellIndex1d].testInBox(*molecule));
		return cellIndex1d;
	}*/
}

unsigned long int LinkedCells::getCellIndexOfPoint(const double point[3]) const {
	int cellIndex[3]; // 3D Cell index
	double localPoint[] = {point[0], point[1], point[2]};

	for (int dim = 0; dim < 3; dim++) {
		// different than getCellIndexOfMolecule!!!

		// ignore a bit of rounding, if the point is outside of the box.
		if (localPoint[dim] <= _haloBoundingBoxMin[dim]){
			localPoint[dim] += _cellLength[dim] * 0.5;
		}
		else if(localPoint[dim] >= _haloBoundingBoxMax[dim]){
			localPoint[dim] -= _cellLength[dim] * 0.5;
		}

		#ifndef NDEBUG
		//this should never ever happen!
		if (localPoint[dim] < _haloBoundingBoxMin[dim] || localPoint[dim] >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "Point is outside of halo bounding box" << endl;
			global_log->error() << "Point p = (" << localPoint[0] << ", " << localPoint[1] << ", " << localPoint[2] << ")" << endl;
			global_log->error() << "_haloBoundingBoxMin = (" << _haloBoundingBoxMin[0] << ", " << _haloBoundingBoxMin[1] << ", " << _haloBoundingBoxMin[2] << ")" << endl;
			global_log->error() << "_haloBoundingBoxMax = (" << _haloBoundingBoxMax[0] << ", " << _haloBoundingBoxMax[1] << ", " << _haloBoundingBoxMax[2] << ")" << endl;
			Simulation::exit(1);
		}
		#endif

		//this version is sensitive to roundoffs, if we have molecules (initialized) precisely at position 0.0:
		//cellIndex[dim] = (int) floor(point[dim] - _haloBoundingBoxMin[dim]) * _cellLengthReciprocal[dim]);
		cellIndex[dim] = min(max(((int) floor((localPoint[dim] - _boundingBoxMin[dim]) * _cellLengthReciprocal[dim])) + _haloWidthInNumCells[dim],0),_cellsPerDimension[dim]-1);
	}

	int cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
	// in very rare cases rounding is stupid, thus we need a check...
	//TODO: check if this can in any way be done better...
	if (_cells[cellIndex1d].testPointInCell(localPoint)) {
		return cellIndex1d;
	}
	else {
		for (int dim = 0; dim < 3; dim++) {
			if(localPoint[dim]<_cells[cellIndex1d].getBoxMin(dim)){
				cellIndex[dim]--;
			}
			else{
				if(localPoint[dim]>=_cells[cellIndex1d].getBoxMax(dim)){
					cellIndex[dim]++;
				}
			}
		}
		cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
		mardyn_assert(_cells[cellIndex1d].testPointInCell(localPoint));
		return cellIndex1d;
	}
}

long int LinkedCells::cellIndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const {
	return (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
}

void LinkedCells::threeDIndexOfCellIndex(int ind, int r[3], const int dim[3]) const {
	r[2] = ind / (dim[0] * dim[1]);
	r[1] = (ind - r[2] * dim[0] * dim[1]) / dim[0];
	r[0] = ind - dim[0] * (r[1] + dim[1] * r[2]);
}

void LinkedCells::getCellIndicesOfRegion(const double startRegion[3], const double endRegion[3], unsigned int &startIndex, unsigned int &endIndex) {
	//find the cellIndices of the cells which contain the start and end corners of the region
	startIndex = getCellIndexOfPoint(startRegion);
	endIndex = getCellIndexOfPoint(endRegion);
}

unsigned long LinkedCells::initCubicGrid(std::array<unsigned long, 3> numMoleculesPerDimension, std::array<double, 3> simBoxLength) {
	const unsigned long numCells = _cells.size();

	std::vector<unsigned long> numMoleculesPerThread;
	const int numThreads = mardyn_get_max_threads();
	numMoleculesPerThread.resize(numThreads);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int myID = mardyn_get_thread_num();
		const unsigned long myStart = numCells * myID / numThreads;
		const unsigned long myEnd = numCells * (myID + 1) / numThreads;
		const int seed = myID;

		unsigned long numMoleculesByThisThread = 0;
		Random threadPrivateRNG = Random(seed);

		// manual "static" scheduling important, because later this thread needs to traverse the same cells
		for (unsigned long cellIndex = myStart; cellIndex < myEnd; ++cellIndex) {
			ParticleCell & cell = _cells[cellIndex];
			numMoleculesByThisThread += cell.initCubicGrid(numMoleculesPerDimension, simBoxLength, threadPrivateRNG);
		}

		// prefix sum of numMoleculesByThisThread


		#if defined(_OPENMP)
		#pragma omp critical(prefixInitCubicGridLC)
		#endif
		for (int i = myID; i < numThreads; ++i) {
			numMoleculesPerThread[i] += numMoleculesByThisThread;
		}

		/* wait for all threads to accumulate */
		#if defined(_OPENMP)
		#pragma omp barrier
		#endif

		// assign local IDs

		unsigned long threadIDsAssignedByThisThread = 0;
		if(myID > 0) {
			threadIDsAssignedByThisThread = numMoleculesPerThread[myID-1];
		}

		// manual "static" scheduling important, because later this thread needs to traverse the same cells
		for (unsigned long cellIndex = myStart; cellIndex < myEnd; ++cellIndex) {
			ParticleCell & cell = _cells[cellIndex];
			const int numMolecules = cell.getMoleculeCount();

			SingleCellIterator<ParticleCell> begin = cell.iterator();

			for (SingleCellIterator<ParticleCell> it = begin; it.isValid(); ++it) {
				it->setid(threadIDsAssignedByThisThread);
				++threadIDsAssignedByThisThread;
			}
		}
	} /* end of parallel */

	unsigned long totalNumberOfMolecules = numMoleculesPerThread.back();
	return totalNumberOfMolecules;
}

RegionParticleIterator LinkedCells::getRegionParticleIterator(
		const double startRegion[3], const double endRegion[3],
		const unsigned int startRegionCellIndex,
		const unsigned int endRegionCellIndex, ParticleIterator::Type type) {
	int start3DIndices[3], end3DIndices[3];
	int regionDimensions[3];
	threeDIndexOfCellIndex(startRegionCellIndex, start3DIndices, _cellsPerDimension);
	threeDIndexOfCellIndex(endRegionCellIndex, end3DIndices, _cellsPerDimension);
	for(int d = 0; d < 3; d++){
		regionDimensions[d] = end3DIndices[d] - start3DIndices[d] + 1;
	}
	ParticleIterator::CellIndex_T offset = mardyn_get_thread_num(); // starting position
	ParticleIterator::CellIndex_T stride = mardyn_get_num_threads(); // stride

	return RegionParticleIterator(type, &_cells, offset, stride, startRegionCellIndex, regionDimensions, _cellsPerDimension, startRegion, endRegion);
}

void LinkedCells::deleteMolecule(Molecule &molecule, const bool& rebuildCaches) {
	auto cellid = getCellIndexOfMolecule(&molecule);

	if (cellid >= _cells.size()) {
		global_log->error_always_output()
				<< "coordinates for atom deletion lie outside bounding box."
				<< endl;
		Simulation::exit(1);
	}

	bool found = this->_cells[cellid].deleteMoleculeByID(molecule.getID());

	if (!found) {
		global_log->error_always_output() << "could not delete molecule " << molecule.getID() << "."
				<< endl;
		Simulation::exit(1);
	}
	else if (rebuildCaches) {
		_cells[cellid].buildSoACaches();
	}
}

double LinkedCells::getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1, CellProcessor& cellProcessorI) {
	CellProcessor* cellProcessor;
	if (dynamic_cast<LegacyCellProcessor*>(&cellProcessorI)) {
		cellProcessor = &cellProcessorI;
	} else {
		cellProcessor = new LegacyCellProcessor(cellProcessorI.getCutoffRadius(), cellProcessorI.getLJCutoffRadius(),
				particlePairsHandler);
	}

	double u = 0.0;

	unsigned long cellIndex = getCellIndexOfMolecule(m1);

	ParticleCell& currentCell = _cells[cellIndex];

	mardyn_assert(not currentCell.isHaloCell());

	Molecule molWithSoA = *m1;
	molWithSoA.buildOwnSoA();

	cellProcessor->initTraversal();

	u += cellProcessor->processSingleMolecule(&molWithSoA, currentCell);

	vector<long> forwardNeighbourOffsets; // now vector
	vector<long> backwardNeighbourOffsets; // now vector
	calculateNeighbourIndices(forwardNeighbourOffsets, backwardNeighbourOffsets);

	// forward neighbours
	for (auto neighbourOffsetsIter = forwardNeighbourOffsets.begin();
			neighbourOffsetsIter != forwardNeighbourOffsets.end();
			neighbourOffsetsIter++) {
		ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
		u += cellProcessor->processSingleMolecule(&molWithSoA, neighbourCell);
	}
	// backward neighbours
	for (auto neighbourOffsetsIter = backwardNeighbourOffsets.begin();
			neighbourOffsetsIter != backwardNeighbourOffsets.end();
			neighbourOffsetsIter++) {
		ParticleCell& neighbourCell = _cells[cellIndex - *neighbourOffsetsIter];
		u += cellProcessor->processSingleMolecule(&molWithSoA, neighbourCell);
	}

	cellProcessor->endTraversal();

	molWithSoA.releaseOwnSoA();

	if (!dynamic_cast<LegacyCellProcessor*>(&cellProcessorI)) {
		delete cellProcessor;
	}

    mardyn_assert(not std::isnan(u)); // catches NaN

	return u;
}

void LinkedCells::updateInnerMoleculeCaches() {
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		if(_cells[cellIndex].isInnerCell()){
			_cells[cellIndex].buildSoACaches();
		}
	}
}

void LinkedCells::updateBoundaryAndHaloMoleculeCaches() {
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		if (_cells[cellIndex].isHaloCell() or _cells[cellIndex].isBoundaryCell()) {
			_cells[cellIndex].buildSoACaches();
		}
	}
}

void LinkedCells::updateMoleculeCaches() {
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		_cells[cellIndex].buildSoACaches();
	}
}

size_t LinkedCells::getTotalSize() {
	size_t totalSize = sizeof(LinkedCells);
	for (auto& cell : _cells) {
		totalSize += sizeof(ParticleCell);
		totalSize += cell.getCellDataSoA().getDynamicSize();
		totalSize += cell.getMoleculeVectorDynamicSize();
	}
	totalSize += _haloCellIndices.capacity() * sizeof(unsigned long);
	return totalSize;
}
void LinkedCells::printSubInfo(int offset) {
	size_t ownSize = sizeof(LinkedCells), cellTotal = 0, cellSoA = 0, cellMoleculeVectors = 0;
	for (auto& cell : _cells) {
		cellTotal += sizeof(ParticleCell);
		cellSoA += cell.getCellDataSoA().getDynamicSize();
		cellMoleculeVectors += cell.getMoleculeVectorDynamicSize();
	}
	cellTotal += cellSoA + cellMoleculeVectors;
	size_t indexVectors = 0;
	indexVectors += _haloCellIndices.capacity() * sizeof(unsigned long);
	std::stringstream offsetstream;
	for (int i = 0; i < offset; i++) {
		offsetstream << "\t";
	}
	global_log->info() << offsetstream.str() << "own datastructures:\t" << ownSize / 1.e6 << " MB" << std::endl;
	global_log->info() << offsetstream.str() << "cells total:\t\t" << cellTotal / 1.e6 << " MB" << std::endl;
	global_log->info() << offsetstream.str() << "cells SoAs:\t\t" << cellSoA / 1.e6 << " MB" << std::endl;
	global_log->info() << offsetstream.str() << "cells molecule vectors:\t" << cellMoleculeVectors / 1.e6 << " MB"
			<< std::endl;
	global_log->info() << offsetstream.str() << "indexVectors:\t\t" << indexVectors / 1.e6 << " MB" << std::endl;
}
std::string LinkedCells::getName() {
	return "LinkedCells";
}

bool LinkedCells::getMoleculeAtPosition(const double pos[3], Molecule** result) {
	const double epsi = this->_cutoffRadius * 1e-6;
	auto index = getCellIndexOfPoint(pos);
	auto& cell = _cells.at(index);

	// iterate through cell and compare position of molecules with given position

	SingleCellIterator<ParticleCell> begin1 = cell.iterator();

	for (SingleCellIterator<ParticleCell> it1 = begin1; it1.isValid(); ++it1) {
		auto& mol = *it1;

		if (fabs(mol.r(0) - pos[0]) <= epsi && fabs(mol.r(1) - pos[1]) <= epsi && fabs(mol.r(2) - pos[2]) <= epsi) {
			// found
			*result = &mol;
			return true;
		}
	}
	// not found
	return false;
}

Molecule* LinkedCells::getHaloMoleculeCloseToPosition(const double *pos, unsigned long id) {
	double lowCorner[3] = {pos[0] - _cutoffRadius, pos[1] - _cutoffRadius, pos[2] - _cutoffRadius};
	double highCorner[3] = {pos[0] + _cutoffRadius, pos[1] + _cutoffRadius, pos[2] + _cutoffRadius};
	for (auto iter = regionIterator(lowCorner, highCorner);
	     iter.isValid(); ++iter) {
		if (iter->getID() == id) {
			return &*iter;
		}
	}
	throw std::runtime_error("no fitting halo molecule found");
}

bool LinkedCells::requiresForceExchange() const {return _traversalTuner->getCurrentOptimalTraversal()->requiresForceExchange();}

std::vector<unsigned long> LinkedCells::getParticleCellStatistics() {
	int maxParticles = 0;
	for (auto& cell : _cells) {
		maxParticles = std::max(maxParticles, cell.getMoleculeCount());
	}

	std::vector<unsigned long> statistics(maxParticles + 1, 0ul);
	for (auto& cell : _cells) {
		statistics[cell.getMoleculeCount()]++;
	}
	return statistics;
}
