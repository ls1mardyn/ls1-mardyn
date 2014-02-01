#include "AdaptiveSubCells.h"

#include <cmath>
#include <iostream>

#include "Domain.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "ParticleCell.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################


AdaptiveSubCells::AdaptiveSubCells(
		double bBoxMin[3], double bBoxMax[3],
		double cutoffRadius, double LJCutoffRadius)
		: ParticleContainer(bBoxMin, bBoxMax)
{
	int numberOfCells = 1;
	_cutoffRadius = cutoffRadius;
	_LJCutoffRadius = LJCutoffRadius;

	for (int d = 0; d < 3; d++) {
		_haloWidthInNumCells[d] = 1;
		_cellsPerDimension[d] = (int) floor((_boundingBoxMax[d] - _boundingBoxMin[d]) / cutoffRadius) + 2*_haloWidthInNumCells[d];
		// in each dension at least one layer of (inner+boundary) cells necessary
		if (_cellsPerDimension[d] == 2 * _haloWidthInNumCells[d]) {
			_cellsPerDimension[d]++;
		}
		numberOfCells *= _cellsPerDimension[d];
		_cellLength[d] = (_boundingBoxMax[d] - _boundingBoxMin[d]) / (_cellsPerDimension[d] - 2*_haloWidthInNumCells[d]);
		_haloBoundingBoxMin[d] = _boundingBoxMin[d] - _haloWidthInNumCells[d] * _cellLength[d];
		_haloBoundingBoxMax[d] = _boundingBoxMax[d] + _haloWidthInNumCells[d] * _cellLength[d];
		_haloLength[d] = _haloWidthInNumCells[d] * _cellLength[d];
	}

	// The resize method initializes the _cells array with the number of all cells within the bounding box
	_cells.resize(numberOfCells);
	// Initialize the vector containing the _localRho value of each coarse cell
	_localRho.resize(numberOfCells);
	// Initialize the vector containing the _metaCellIndex of each coarse cell
	// The size of _metaCellIndex is by one larger than _cells.
	// To find out, whether a cell is refined, the _subCell-index difference of
	// the cell and the next coarse cell is calculated. If the difference is 8, it
	// is a refined cell, if the difference is 1, it is a coarse cell. To be able
	// to do this calculation for the last cell, the index of the "next cell" has
	// to be calculated. To be able to store this index, the vector _metaCellIndex
	// has to be one element larger.
	_metaCellIndex.resize(numberOfCells + 1);

	// If the with of the inner region is less than the width of the halo region
	// a parallelisation isn't possible (with the used algorithms).
	// In this case, print an error message
	// _cellsPerDimension is 2 times the halo width + the inner width
	// so it has to be at least 3 times the halo width
	if (_cellsPerDimension[0] < 3*_haloWidthInNumCells[0] ||
	    _cellsPerDimension[1] < 3*_haloWidthInNumCells[1] ||
	    _cellsPerDimension[2] < 3*_haloWidthInNumCells[2])
	{
		global_log->error() << "Error in AdaptiveSubCells (Constructor): bounding box too small for calculated cell Length" << endl;
		global_log->error() << "cellsPerDimension" << _cellsPerDimension[0] << " / " << _cellsPerDimension[1] << " / " << _cellsPerDimension[2] << endl;
		global_log->error() << "_haloWidthInNumCells" << _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1] << " / " << _haloWidthInNumCells[2] << endl;
		exit(5);
	}
	// initially, zero updates were performed:
	_numberOfUpdates = 0;
	_cellsValid = false;
}


AdaptiveSubCells::~AdaptiveSubCells() {
	// empty
}

void AdaptiveSubCells::readXML(XMLfileUnits& xmlconfig) {
	/* no parameters */
}

void AdaptiveSubCells::rebuild(double bBoxMin[3], double bBoxMax[3]) {
	for (int d = 0; d < 3; d++) {
		_boundingBoxMin[d] = bBoxMin[d];
		_boundingBoxMax[d] = bBoxMax[d];
	}

	int numberOfCells = 1;

	for (int d = 0; d < 3; d++) {
		_cellsPerDimension[d] = (int) floor((_boundingBoxMax[d] - _boundingBoxMin[d]) / (_cutoffRadius / _haloWidthInNumCells[d]))
		    + 2 * _haloWidthInNumCells[d];
		// in each dimension at least one layer of (inner+boundary) cells necessary
		if (_cellsPerDimension[d] == 2 * _haloWidthInNumCells[d]) {
			_cellsPerDimension[d]++;
		}
		numberOfCells *= _cellsPerDimension[d];
		_cellLength[d] = (_boundingBoxMax[d] - _boundingBoxMin[d]) / (_cellsPerDimension[d] - 2 * _haloWidthInNumCells[d]);
		_haloBoundingBoxMin[d] = _boundingBoxMin[d] - _haloWidthInNumCells[d] * _cellLength[d];
		_haloBoundingBoxMax[d] = _boundingBoxMax[d] + _haloWidthInNumCells[d] * _cellLength[d];
		_haloLength[d] = _haloWidthInNumCells[d] * _cellLength[d];
	}

	// The resize method initializes the _cells array with the number of all cells within the bounding box
	_cells.resize(numberOfCells);

	// Initialize the vector containing the _localRho value of each coarse cell
	_localRho.resize(numberOfCells);
	// Initialize the vector containing the _metaCellIndex of each coarse cell
	// The size of _metaCellIndex is by one larger than _cells.
	// To find out, whether a cell is refined, the _subCell-index difference of
	// the cell and the next coarse cell is calculated. If the difference is 8, it
	// is a refined cell, if the difference is 1, it is a coarse cell. To be able
	// to do this calculation for the last cell, the index of the "next cell" has
	// to be calculated. To be able to store this index, the vector _metaCellIndex
	// has to be one element larger.
	_metaCellIndex.resize(numberOfCells+1);

	// If the with of the inner region is less than the width of the halo region
	// a parallelisation isn't possible (with the used algorithms).
	// In this case, print an error message
	// _cellsPerDimension is 2 times the halo width + the inner width
	// so it has to be at least 3 times the halo width
	if (_cellsPerDimension[0] < 3*_haloWidthInNumCells[0] ||
	    _cellsPerDimension[1] < 3*_haloWidthInNumCells[1] ||
	    _cellsPerDimension[2] < 3*_haloWidthInNumCells[2]) {
		global_log->error() << "Error in AdaptiveSubCells::rebuild: bounding box too small for calculated cell Length" << endl;
		global_log->error() << "cellsPerDimension" << _cellsPerDimension[0] << " / " << _cellsPerDimension[1] << " / " << _cellsPerDimension[2] << endl;
		global_log->error() << "_haloWidthInNumCells" << _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1] << " / " << _haloWidthInNumCells[2] << endl;
		exit(5);
	}
	// initially, zero updates were perfored:
	_numberOfUpdates = 0;

	// delete all Particles which are outside of the halo region
	std::list<Molecule>::iterator particleIterator = _particles.begin();
	bool erase_mol;
	while (particleIterator != _particles.end()) {
		erase_mol = false;
		for (int d = 0; d < 3; d++) {
			const double& rd = particleIterator->r(d);
			// The molecules has to be within the domain of the process
			// If it is outside in at least one dimension, it has to be
			// erased /
			if (rd < _haloBoundingBoxMin[d] || rd >= _haloBoundingBoxMax[d])
				erase_mol = true;
		}
		if (erase_mol) {
			particleIterator = _particles.erase(particleIterator);
		}
		else {
			particleIterator++;
		}
	}
	_cellsValid = false;
}


void AdaptiveSubCells::update() {

	// the _cells vector is only needed in order to compute _localRho[] for each coarse grid cell. So we need to initialize _cells in 1 of 1000 timesteps only
	// FIXME: introduce global variabel instead of fixed numer '1000'
	if (_numberOfUpdates % 1000 == 0) {
		// clear all Cells
		std::vector<ParticleCell>::iterator celliter;
		for (celliter = (_cells).begin(); celliter != (_cells).end(); ++celliter) {
			(*celliter).removeAllParticles();
		}
		unsigned long index; // index of the cell into which the pointer has to be inserted
		std::list<Molecule>::iterator pos;
		for (pos = _particles.begin(); pos != _particles.end(); pos++) {
			index = getCellIndexOfMolecule(&(*pos));
			if (index >= _cells.size()) { // FIXME: how likely/critical is this, can this go to an assert?
				global_log->error() << "Found wrong index in AdaptiveSubCells::update()" << endl;
				exit(1);
			}
			(_cells[index]).addParticle(&(*pos));
		}
		// computes the _metaCellIndex for each Cell in the coarse grid
		calculateMetaCellIndex();
		int numberOfSubCells = _metaCellIndex[_cells.size()];
		// The .clear() method clears the whole _subCells vector
		_subCells.clear();
		// The resize method initializes the _subCells array with the number of all subCells within the bounding box
		_subCells.resize(numberOfSubCells);
		_forwardNeighbourSubOffsets.resize(numberOfSubCells);
		_backwardNeighbourSubOffsets.resize(numberOfSubCells);

		initializeSubCells();
		calculateSubNeighbourIndices();
	}
	// Here the adaptive datastructure is going to be created
	// subcelliter is the new iterator to pass through all subCells
	std::vector<ParticleCell>::iterator subcelliter;
	for (subcelliter = (_subCells).begin(); subcelliter != (_subCells).end(); subcelliter++) {
		(*subcelliter).removeAllParticles();
	}
	// subIndex is the index of the cell into which the pointer has to be inserted
	int subIndex;
	std::list<Molecule>::iterator pos;
	for (pos = _particles.begin(); pos != _particles.end(); pos++) {
		// getSubCellIndexOfMolecule computes the index of the subCell where the molecule has to be inserted
		subIndex = getSubCellIndexOfMolecule(&(*pos));

		if (subIndex < (int) 0 || subIndex >= (int) _subCells.size()) { // FIXME: should this go to an asssert?
			global_log->error() << "Found invalid index in AdaptiveSubCells::update()." << endl;
			exit(1);
		}
		_subCells[subIndex].addParticle(&(*pos));
	}
	// increas number of performed updates
	_numberOfUpdates++;
	_cellsValid = true;
}

void AdaptiveSubCells::addParticle(Molecule& particle) {

	double x = particle.r(0);
	double y = particle.r(1);
	double z = particle.r(2);

	if (x >= _haloBoundingBoxMin[0] && x < _haloBoundingBoxMax[0] &&
	    y >= _haloBoundingBoxMin[1] && y < _haloBoundingBoxMax[1] &&
	    z >= _haloBoundingBoxMin[2] && z < _haloBoundingBoxMax[2]) {
		_particles.push_front(particle);
		if (_cellsValid) {
			int subIndex = getSubCellIndexOfMolecule(&particle);
			if (subIndex < (int) 0 || subIndex >= (int) _subCells.size()) // FIXME: Should this go to an assert?
			{
				global_log->error() << "Found invalid index in AdaptiveSubCells::addParticle()" << endl;
				exit(1);
			}
			_subCells[subIndex].addParticle(&(_particles.front()));
		}
	}
}


unsigned AdaptiveSubCells::countParticles(unsigned int cid) {
	unsigned N = 0;
	std::vector<Molecule*>::iterator molIter1;

	for (unsigned long i = 0; i < _cells.size(); i++) {
		ParticleCell& currentCell = _cells[i];
		if (currentCell.isHaloCell())
			continue;
		for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
			if ((*molIter1)->componentid() == cid)
				N++;
		}
	}
	return N;
}

/**
 * @todo move this method to the ChemicalPotential, using a call to ParticleContainer::getRegion() !?
 */
unsigned AdaptiveSubCells::countParticles(unsigned int cid, double* cbottom, double* ctop) {
	int minIndex[3];
	int maxIndex[3];

	for (int d = 0; d < 3; d++) {
		if (cbottom[d] < _haloBoundingBoxMin[d])
			minIndex[d] = 0;
		else
			minIndex[d] = (int) floor((cbottom[d] - _haloBoundingBoxMin[d]) / _cellLength[d]);

		if (ctop[d] > _haloBoundingBoxMax[d])
			maxIndex[d] = (int) floor((_haloBoundingBoxMax[d] - _haloBoundingBoxMin[d]) / _cellLength[d]);
		else
			maxIndex[d] = (int) floor((ctop[d] - _haloBoundingBoxMin[d]) / _cellLength[d]);
	}

	unsigned N = 0;
	int cix[3];
	std::vector<Molecule*>::iterator molIter1;
	bool individualCheck;
	int cellid;

	for (cix[0] = minIndex[0]; maxIndex[0] >= cix[0]; (cix[0])++) {
		for (cix[1] = minIndex[1]; maxIndex[1] >= cix[1]; (cix[1])++) {
			for (cix[2] = minIndex[2]; maxIndex[2] >= cix[2]; (cix[2])++) {
				individualCheck = (cix[0] == minIndex[0]) || (cix[0] == maxIndex[0]) || (cix[1] == maxIndex[1]) || (cix[1] == maxIndex[1]) || (cix[2] == maxIndex[2]) || (cix[2] == maxIndex[2]);
				cellid = this->subCellIndexOf3DIndex(cix[0], cix[1], cix[2]);
				ParticleCell& currentCell = _subCells[cellid];

				if (currentCell.isHaloCell())
					continue;

				if (individualCheck) {
					for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
						if (((*molIter1)->r(0) > cbottom[0]) && ((*molIter1)->r(1) > cbottom[1]) && ((*molIter1)->r(2) > cbottom[2]) && ((*molIter1)->r(0) < ctop[0]) && ((*molIter1)->r(1) < ctop[1]) && ((*molIter1)->r(2) < ctop[2]) && ((*molIter1)->componentid() == cid))
							N++;
					}
				}
				else {
					for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
						if ((*molIter1)->componentid() == cid)
							N++;
					}
				}
			}
		}
	}

	return N;
}

void AdaptiveSubCells::deleteMolecule(unsigned long molid, double x, double y, double z) {
	int ix = (int) floor((x - _haloBoundingBoxMin[0]) / _cellLength[0]);
	int iy = (int) floor((y - _haloBoundingBoxMin[1]) / _cellLength[1]);
	int iz = (int) floor((z - _haloBoundingBoxMin[2]) / _cellLength[2]);
	unsigned long hash = this->subCellIndexOf3DIndex(ix, iy, iz);
	if (hash >= _subCells.size()) {
		global_log->error() << "SEVERE ERROR: coordinates for atom deletion lie outside bounding box.\n";
		exit(1);
	}
	bool found = _subCells[hash].deleteMolecule(molid);
	if (!found) {
		global_log->error() << "SEVERE ERROR: could not delete molecule " << molid << ".\n";
		exit(1);
	}
}


void AdaptiveSubCells::traverseCells(CellProcessor& cellProcessor) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (traversePairs) invalid, call update first" << endl;
		exit(1);
	}

	vector<unsigned long>::iterator neighbourOffsetsIter;

#ifndef NDEBUG
	global_log->debug() << "LinkedCells::traverseCells: Processing pairs and preprocessing Tersoff pairs." << endl;
	global_log->debug() << "_minNeighbourOffset=" << _minNeighbourOffset << "; _maxNeighbourOffset=" << _maxNeighbourOffset<< endl;
#endif

	cellProcessor.initTraversal(_maxNeighbourOffset + _minNeighbourOffset +1);
	// open the window of cells activated
	for (unsigned int cellIndex = 0; cellIndex < _maxNeighbourOffset; cellIndex++) {
		cellProcessor.preprocessCell(_subCells[cellIndex]);
	}

	// loop over all inner cells and calculate forces to forward neighbours
	for (unsigned int cellIndex = 0; cellIndex < _subCells.size(); cellIndex++) {
		ParticleCell& currentCell = _subCells[cellIndex];

		// extend the window of cells with cache activated
		if (cellIndex + _maxNeighbourOffset < _subCells.size()) {
			#ifndef NDEBUG
			global_log->debug() << "Opening cached cells window for cell index=" << (cellIndex + _maxNeighbourOffset)
					<< " with numMolecules()="<< _subCells[cellIndex + _maxNeighbourOffset].getMoleculeCount()
					<< " currentCell " << cellIndex << endl;
			#endif
			cellProcessor.preprocessCell(_subCells[cellIndex + _maxNeighbourOffset]);
		}

		if (currentCell.isInnerCell()) {
			cellProcessor.processCell(currentCell);
			// loop over all neighbours
			for (neighbourOffsetsIter = _forwardNeighbourSubOffsets[cellIndex].begin(); neighbourOffsetsIter != _forwardNeighbourSubOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _subCells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}

		if (currentCell.isHaloCell()) {
			cellProcessor.processCell(currentCell);
			for (neighbourOffsetsIter = _forwardNeighbourSubOffsets[cellIndex].begin(); neighbourOffsetsIter != _forwardNeighbourSubOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
				if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) (_subCells.size())))
					continue;
				ParticleCell& neighbourCell = _subCells[neighbourCellIndex];
				if (!neighbourCell.isHaloCell())
					continue;

				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}

		// loop over all boundary cells and calculate forces to forward and backward neighbours
		if (currentCell.isBoundaryCell()) {
			cellProcessor.processCell(currentCell);

			// loop over all forward neighbours
			for (neighbourOffsetsIter = _forwardNeighbourSubOffsets[cellIndex].begin(); neighbourOffsetsIter != _forwardNeighbourSubOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _subCells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}

			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for (neighbourOffsetsIter = _backwardNeighbourSubOffsets[cellIndex].begin(); neighbourOffsetsIter != _backwardNeighbourSubOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _subCells[cellIndex + *neighbourOffsetsIter];  // minus oder plus?
				if (neighbourCell.isHaloCell()) {
					cellProcessor.processCellPair(currentCell, neighbourCell);
				}
			}
		} // if ( isBoundaryCell() )

		// narrow the window of cells activated
		if (cellIndex >= _minNeighbourOffset) {
#ifndef NDEBUG
			global_log->debug() << "Narrowing cells window for cell index=" << (cellIndex - _minNeighbourOffset)
									<< " with size()="<<_subCells[cellIndex - _minNeighbourOffset].getMoleculeCount()
									<< " currentCell " << cellIndex << endl;
#endif
			cellProcessor.postprocessCell(_subCells[cellIndex - _minNeighbourOffset]);
		}
	} // loop over all cells

	// close the window of cells with cache activated
	for (unsigned int cellIndex = _subCells.size() - _minNeighbourOffset; cellIndex < _subCells.size(); cellIndex++) {
#ifndef NDEBUG
			global_log->debug() << "Narrowing cached cells window for cell index=" << cellIndex
					<< " size()="<<_subCells[cellIndex].getMoleculeCount() << endl;
#endif
			cellProcessor.postprocessCell(_subCells[cellIndex]);
	}
	cellProcessor.endTraversal();
}

double AdaptiveSubCells::getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1, CellProcessor& cellProcessor)
{
	double u = 0.0;

	unsigned long subCellIndex = getSubCellIndexOfMolecule(m1);
	ParticleCell& currentSubCell = _subCells[subCellIndex];
	vector<unsigned long>::iterator neighbourOffsetsIter;
	cellProcessor.initTraversal(_maxNeighbourOffset + _minNeighbourOffset + 1);

	// extend the window of cells with cache activated
	for (unsigned int windowCellIndex = subCellIndex - _minNeighbourOffset; windowCellIndex < subCellIndex + _maxNeighbourOffset+1 ; windowCellIndex++) {
		cellProcessor.preprocessCell(_cells[windowCellIndex]);
	}

	if (m1->numTersoff() > 0) {
		global_log->error() << "The grand canonical ensemble is not implemented for solids.\n";
		exit(848);
	}
	u += cellProcessor.processSingleMolecule(m1, currentSubCell);

	// forward neighbours
	for (neighbourOffsetsIter = _forwardNeighbourSubOffsets[subCellIndex].begin(); neighbourOffsetsIter != _forwardNeighbourSubOffsets[subCellIndex].end(); neighbourOffsetsIter++)
	{
		ParticleCell& neighbourCell = _subCells[subCellIndex + *neighbourOffsetsIter];
		u += cellProcessor.processSingleMolecule(m1, neighbourCell);
	}
	// backward neighbours
	for (neighbourOffsetsIter = _backwardNeighbourSubOffsets[subCellIndex].begin(); neighbourOffsetsIter != _backwardNeighbourSubOffsets[subCellIndex].end(); neighbourOffsetsIter++)
	{
		ParticleCell& neighbourCell = _subCells[subCellIndex + *neighbourOffsetsIter];  // minus oder plus?
		u += cellProcessor.processSingleMolecule(m1, neighbourCell);
	}
	
	// close the window of cells activated
	for (unsigned int windowCellIndex = subCellIndex - _minNeighbourOffset; windowCellIndex < subCellIndex + _maxNeighbourOffset+1; windowCellIndex++) {
		cellProcessor.postprocessCell(_cells[windowCellIndex]);
	}

	cellProcessor.endTraversal();
	return u;
}

unsigned long AdaptiveSubCells::getNumberOfParticles() {
	return _particles.size();
}

// this method remains unchanged

Molecule* AdaptiveSubCells::begin() {
	_particleIter = _particles.begin();
	if (_particleIter != _particles.end()) {
		return &(*_particleIter);
	}
	else {
		return NULL;
	}
}

// this method remains unchanged

Molecule* AdaptiveSubCells::next() {
	_particleIter++;
	if (_particleIter != _particles.end()) {
		return &(*_particleIter);
	}
	else {
		return NULL;
	}
}

// this method remains unchanged

Molecule* AdaptiveSubCells::end() {
	return NULL;
}

// this method remains unchanged

Molecule* AdaptiveSubCells::deleteCurrent() {
	_particleIter = _particles.erase(_particleIter);
	if (_particleIter != _particles.end()) {
		return &(*_particleIter);
	}
	else {
		return NULL;
	}
}

// this method remains unchanged

void AdaptiveSubCells::deleteOuterParticles() {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in AdaptiveSubCells (deleteOuterParticles) invalid, call update first" << endl;
		exit(1);
	}

	vector<unsigned long>::iterator subCellIndexIter;
	//std::list<Molecule*>::iterator molIter1;
	for (subCellIndexIter = _haloSubCellIndices.begin(); subCellIndexIter != _haloSubCellIndices.end(); subCellIndexIter++) {
		ParticleCell& currentSubCell = _subCells[*subCellIndexIter];
		currentSubCell.removeAllParticles();
	}

	std::list<Molecule>::iterator particleIterator = _particles.begin();
	bool erase_mol;
	while (particleIterator != _particles.end()) {
		erase_mol = false;
		for (unsigned short d = 0; d < 3; ++d) {
			const double& rd = particleIterator->r(d);
			// The molecules has to be within the domain of the process
			// If it is outside in at least one dimension, it has to be
			// erased /
			if (rd < _boundingBoxMin[d] || rd >= _boundingBoxMax[d])
				erase_mol = true;
		}
		if (erase_mol) {
			particleIterator = _particles.erase(particleIterator);
		}
		else {
			particleIterator++;
		}
	}
}

// this method remains unchanged

double AdaptiveSubCells::get_halo_L(int index) const {
	return _haloLength[index];
}


void AdaptiveSubCells::getHaloParticles(list<Molecule*> &haloParticlePtrs) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in AdaptiveSubCells (getHaloParticles) invalid, call update first" << endl;
		exit(1);
	}

	std::vector<Molecule*>::iterator particlePtrIter;
	vector<unsigned long>::iterator subCellIndexIter;

	// loop over all halo cells
	for (subCellIndexIter = _haloSubCellIndices.begin(); subCellIndexIter != _haloSubCellIndices.end(); subCellIndexIter++) {
		ParticleCell& currentSubCell = _subCells[*subCellIndexIter];
		// loop over all molecules in the cell
		for (particlePtrIter = currentSubCell.getParticlePointers().begin(); particlePtrIter != currentSubCell.getParticlePointers().end(); particlePtrIter++) {
			haloParticlePtrs.push_back(*particlePtrIter);
		}
	}
}

void AdaptiveSubCells::getRegion(double lowCorner[3], double highCorner[3], list<Molecule*> &particlePtrs) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in AdaptiveSubCells (getRegion) invalid, call update first" << endl;
		exit(1);
	}
	int startIndex[3];
	int stopIndex[3];
	int globalCellIndex;
	std::vector<Molecule*>::iterator particleIter;

	for (int dim = 0; dim < 3; dim++) {
		if (lowCorner[dim] < _boundingBoxMax[dim] && highCorner[dim] > _boundingBoxMin[dim]) {
			startIndex[dim] = (int) floor((lowCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) - 1;
			stopIndex[dim] = (int) floor((highCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) + 1;
			if (startIndex[dim] < 0)
				startIndex[dim] = 0;
			if (stopIndex[dim] > _cellsPerDimension[dim] - 1)
				stopIndex[dim] = _cellsPerDimension[dim] - 1;
		}
		else {
			// No Part of the given region is owned by this process
			// --> chose some startIndex which is higher than the stopIndex
			startIndex[dim] = 1;
			stopIndex[dim] = 0;
		}
	}

	for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
		for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
			for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
				// globalCellIndex is the cellIndex of the molecule on the coarse Cell level.
				globalCellIndex = (iz * _cellsPerDimension[1] + iy) * _cellsPerDimension[0] + ix;
				// loop over all subcells (either 1 or 8)
				for (int sCIdx = _metaCellIndex[globalCellIndex]; sCIdx < _metaCellIndex[globalCellIndex + 1]; sCIdx++) {
					// traverse all molecules in the current cell
					for (particleIter = _subCells[sCIdx].getParticlePointers().begin(); particleIter != _subCells[sCIdx].getParticlePointers().end(); particleIter++) {
						if ((*particleIter)->r(0) >= lowCorner[0] && (*particleIter)->r(0) < highCorner[0] && (*particleIter)->r(1) >= lowCorner[1] && (*particleIter)->r(1) < highCorner[1] && (*particleIter)->r(2) >= lowCorner[2] && (*particleIter)->r(2) < highCorner[2]) {
							particlePtrs.push_back(*particleIter);
						}
					}
				}
			}
		}
	}
}

//################################################
//############ PRIVATE METHODS ###################
//################################################

// assigns each subCell to the inner, boundary or halo region

void AdaptiveSubCells::initializeSubCells() {
	_innerSubCellIndices.clear();
	_boundarySubCellIndices.clear();
	_haloSubCellIndices.clear();
	int subCellIndex, nextSubCellIndex;
	for (int iz = 0; iz < _cellsPerDimension[2]; ++iz) {
		for (int iy = 0; iy < _cellsPerDimension[1]; ++iy) {
			for (int ix = 0; ix < _cellsPerDimension[0]; ++ix) {
				// calls the method that computes the subCellIndex of the current coarse Cell
				// maybe there's an inline method possible!
				subCellIndex = subCellIndexOf3DIndex(ix, iy, iz);
				// nextSubCellIndex is the index of the subCell that corresponds with (cellIndex + 1)
				nextSubCellIndex = subCellIndexOf3DIndex(ix + 1, iy, iz);
				if (ix < _haloWidthInNumCells[0] || iy < _haloWidthInNumCells[1] || iz < _haloWidthInNumCells[2] || ix >= _cellsPerDimension[0] - _haloWidthInNumCells[0] || iy >= _cellsPerDimension[1] - _haloWidthInNumCells[1] || iz >= _cellsPerDimension[2] - _haloWidthInNumCells[2]) {
					// assign the (first) subCell to the halo region
					_subCells[subCellIndex].assignCellToHaloRegion();
					_haloSubCellIndices.push_back(subCellIndex);
					// if nextSubCellIndex-subCellIndex == 8 (or 7, at the end of the array) the current Cell contains subCells
					if (nextSubCellIndex - subCellIndex > 1) {
						for (int i = 1; i < 8; i++) {
							// assign the remaining 7 subCells to the halo region
							_subCells[subCellIndex + i].assignCellToHaloRegion();
							_haloSubCellIndices.push_back(subCellIndex);
						}
					}
				}
				else if (ix < 2 * _haloWidthInNumCells[0] || iy < 2 * _haloWidthInNumCells[1] || iz < 2 * _haloWidthInNumCells[2] || ix >= _cellsPerDimension[0] - 2 * _haloWidthInNumCells[0] || iy >= _cellsPerDimension[1] - 2 * _haloWidthInNumCells[1] || iz >= _cellsPerDimension[2] - 2 * _haloWidthInNumCells[2]) {
					// assign the (first) subCell to the boundary region
					_subCells[subCellIndex].assignCellToBoundaryRegion();
					_boundarySubCellIndices.push_back(subCellIndex);
					// if nextSubCellIndex-subCellIndex > 1 the current Cell contains subCells
					if (nextSubCellIndex - subCellIndex > 1) {
						for (int i = 1; i < 8; i++) {
							// assign the remaining 7 subCells to the boundary region
							_subCells[subCellIndex + i].assignCellToBoundaryRegion();
							_boundarySubCellIndices.push_back(subCellIndex + i);
						}
					}
				}
				else {
					// assign the (first) subCell to the inner region
					_subCells[subCellIndex].assignCellToInnerRegion();
					_innerSubCellIndices.push_back(subCellIndex);
					// if nextSubCellIndex-subCellIndex > 1 the current Cell contains subCells
					if (nextSubCellIndex - subCellIndex > 1) {
						for (int i = 1; i < 8; i++) {
							// assign the remaining 7 subCells to the inner region
							_subCells[subCellIndex + i].assignCellToInnerRegion();
							_innerSubCellIndices.push_back(subCellIndex + i);
						}
					}
				}
			}
		}
	}
}

void AdaptiveSubCells::calculateSubNeighbourIndices() {
	// The distance between two coarse cells (the centers of the cells) in each coordinate
	// direction can be expressed using the unit "coarse cells" (Distance to a direct
	// neighbour 1, e.g.) Now the distance between find cells (of fine and coarse cells)
	// shall be expressed using the same unit. The centre of a subcell is shifted by
	// a absolute value of 0.25 in each coordinate direction compared to the centre of
	// the corresponding coarse cell. The following table shows the shift for subcell
	// 0 to 7 and all three coordinate direction.
	// dir \ k|   0   |   1   |   2   |   3   |   4   |   5   |   6   |   7   |
	//   x    | -0.25 | +0.25 | -0.25 | +0.25 | -0.25 | +0.25 | -0.25 | +0.25 |
	//   y    | -0.25 | -0.25 | +0.25 | +0.25 | -0.25 | -0.25 | +0.25 | +0.25 |
	//   z    | -0.25 | -0.25 | -0.25 | -0.25 | +0.25 | +0.25 | +0.25 | +0.25 |
	// those values are stored in subCellShifts and are needed for the distance calculation
	double subCellShifts[3][8] = { { -0.25, +0.25, -0.25, +0.25, -0.25, +0.25, -0.25, +0.25 }, { -0.25, -0.25, +0.25, +0.25, -0.25, -0.25, +0.25, +0.25 }, { -0.25, -0.25, -0.25, -0.25, +0.25, +0.25, +0.25, +0.25 } };
	// The values of the following variables {x,y,z}CoarseCellDistanceNotZero are used for the
	// distance calculations of the subcells and have two possible values: 1.0 and 0.0
	// If the corresponding coarse cells don't have a distance in a given coordinate direction
	// (e.g. xIndex = 0), then the minimal distance between the subcells is also Zero,
	// and therefore {x,y,z}CoarseCellDistanceNotZero is 0.0 e.g.:
	// ========================
	// || a3 | a4 || b3 | b4 ||
	// ||----A----||----B----||
	// || a1 | a2 || b1 | b2 ||
	// ========================
	// In this example, The distance, yIndex=0, so the centers of A and B have the same y-value
	// The minimal distance in y-direction between two subcells is in this case always zero
	// in all other cases {x,y,z}CoarseCellDistanceNotZero is 1.0.
	double xCoarseCellDistanceNotZero, yCoarseCellDistanceNotZero, zCoarseCellDistanceNotZero;
	// the minimum squared distance between 2 cells for each spacial dimension is stored in
	// the variables  xDistanceSquare, yDistanceSquare and zDistanceSquare.
	// The spacial offsets for the local subCells (...K) and the neighbouring
	// subCells (...P) must be taken into regard when calculating the distance.
	// Consider the two coarse cells A and B in the following picture.
	// ========================
	// || a3 | a4 || b3 | b4 ||
	// ||----A----||----B----||
	// || a1 | a2 || b1 | b2 ||
	// ========================
	// Both cells have been refined.
	// xIndex can be seen as distance (unit: cells) between A and B, which
	// is in this case 1. Now let's consider the distance between a1 and b2.
	// The centre of SubCell a1 is shifted by -0.25 compared to the centre
	// of the coarse Cell A, this shift is stored in subCellShift[0][k].
	// For b2, a shift of 0.25 is stored in subCellShift[0][p]. The distance
	// between b2 and a1 is abs(position_of_b2 - position_of_a1)
	// Let position_of_b2 be (xIndex+subCellShift[0][p]) = (1+0.25) = 1.25.
	// Then position_of_a1 is subCellShift[0][k] = -0.25
	// ==> distance = abs(position_of_b2 - position_of_a1) = abs(1.25 - (- 0.25))
	//              = 1.5
	// From this value, 0.5 has to be substracted (but only if the distance is not
	// zero, so if e.g. xIndex is not null, see variable xCoarseCellDistanceNotZero)
	// as 1.5 is the distance between the to cell centers. The distance of the center
	// to the border is 0.25 for subCells. For the "real" distance, the resulting
	// value (here 1.0) has to be multiplied with the cellLength.
	double xDistanceSquare, yDistanceSquare, zDistanceSquare;
	int coarseCellOffset; // Each cell is either a coarse cell or is located in a coarse cell. For
	// two neighbouring cells, coarseCellOffset is the index offset of the
	// corresponding coarse cells in the vector _cells
	int metaCellOffset; // Each coarse cell (or "meta" cell if the cell is refined) has also an index in
	// the vector _subCells, which contains coarse and fine cells. metaCellOffset is
	// the index offset of two coarse of "meta" cellsin the vector _subCells
	int SubCellOffset; // The offset between two arbitrary (coarse or fine) cells in the vector _subCells
	double cutoffRadiusSquare = pow(_cutoffRadius, 2); // square of the _cutoffRadius

	// clear the old neighbour offsets of the last adaptive grid in order to create room for the new offsets
	for (unsigned int f = 0; f < _subCells.size(); f++) {
		_forwardNeighbourSubOffsets[f].clear();
		_backwardNeighbourSubOffsets[f].clear();
	}
	_maxNeighbourOffset = 0;
	_minNeighbourOffset = 0;
	// loop over all cellIndices of the coarse grid
	// @todo only necessary cells!
	for (unsigned int i = 0; i < _cells.size(); i++) {
		// Bugfix Johannes Weissl: Exclude halo, otherwise negative index access in _metaCellIndex
		if (_subCells[_metaCellIndex[i]].isHaloCell())
			continue;
		// If two coarse cells, which have a index difference of 1 (cell i and cell /i+1) in the
		// array _cells (usually direct neighbours in x-direction) have an index difference of
		// 8 in the array _subCells (corresponding indeces stored in _metaCellIndex, then the
		// cell i is refined and countains 8 subcells. Otherwise the index difference is 1 and
		// cell i is not refined.
		if (_metaCellIndex[i + 1] - _metaCellIndex[i] == 8) {
			// loop over all local subCells
			for (int k = 0; k < 8; k++) {
				// first find the neighbours whithin the cell itself (7 other subCells)
				// compute the forward index offsets of the local subCells
				for (int l = k + 1; l < 8; l++) {
					SubCellOffset = (int) (l - k);
					_forwardNeighbourSubOffsets[_metaCellIndex[i] + k].push_back(SubCellOffset);
					assert(SubCellOffset > 0);
					if ((unsigned) SubCellOffset > _maxNeighbourOffset) {
						_maxNeighbourOffset = SubCellOffset;
					}
				}
				// compute the backward index offsets of the local subCells
				for (int m = 0; m < k; m++) {
					SubCellOffset = (int) (m - k);
					_backwardNeighbourSubOffsets[_metaCellIndex[i] + k].push_back(SubCellOffset);
					if ((unsigned) abs(SubCellOffset) > _minNeighbourOffset) {
						_minNeighbourOffset = abs(SubCellOffset);
					}
				}
				// loop over all neighbouring ("continue" in case of the cell itself) coarse cells
				for (int zIndex = -_haloWidthInNumCells[2]; zIndex <= _haloWidthInNumCells[2]; zIndex++) {
					for (int yIndex = -_haloWidthInNumCells[1]; yIndex <= _haloWidthInNumCells[1]; yIndex++) {
						for (int xIndex = -_haloWidthInNumCells[0]; xIndex <= _haloWidthInNumCells[0]; xIndex++) {
							// The cell itself should not be searched for neighbours!
							if (xIndex == 0 && yIndex == 0 && zIndex == 0)
								continue;
							// compute the global coarseCellOffset for the current coarse grid neighbour cell
							coarseCellOffset = (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
							// _cells.size()-1 is the upper cellIndex boundary! So i + coarseCellOffset must be smaller than _cells.size()
							// @todo only calculate neighbour offets for necessary cells, so this  if  could then be removed
							if ((int) i + coarseCellOffset < (int) _cells.size()) {
								// check for each dimension if the cells have the same coordinate value in that dimension
								if (xIndex == 0)
									xCoarseCellDistanceNotZero = 0.0;
								else
									xCoarseCellDistanceNotZero = 1.0;
								if (yIndex == 0)
									yCoarseCellDistanceNotZero = 0.0;
								else
									yCoarseCellDistanceNotZero = 1.0;
								if (zIndex == 0)
									zCoarseCellDistanceNotZero = 0.0;
								else
									zCoarseCellDistanceNotZero = 1.0;

								// Index Offset of the two coarse cells in the array _subCells
								metaCellOffset = _metaCellIndex[i + coarseCellOffset] - _metaCellIndex[i];
								// if "delta"==8 the neighbour Cell with index i+coarseCellOffset contains subCells
								// ##################################
								// ## CASE 1: BOTH CELLS REFINED   ##
								// ##################################
								if (_metaCellIndex[i + coarseCellOffset + 1] - _metaCellIndex[i + coarseCellOffset] == 8) {
									// loop over all subCells that the neighbour Cell contains
									for (int p = 0; p < 8; p++) {
										// compute the minimum distance between 2 cells for each spacial dimension.
										// The spacial offsets for the local subCells (...K) and the neighbouring
										// subCells (...P) must be taken into regard.
										// Consider the two coarse cells A and B in the following picture.
										// ========================
										// || a3 | a4 || b3 | b4 ||
										// ||----A----||----B----||
										// || a1 | a2 || b1 | b2 ||
										// ========================
										// xIndex can be seen as distance (unit: cells) between A and B, which
										// is in this case 1. Now let's consider the distance between a1 and b2.
										// The centre of SubCell a1 is shifted by -0.25 compared to the centre
										// of the coarse Cell A, this shift is stored in subCellShift[0][k].
										// For b2, a shift of 0.25 is stored in subCellShift[0][p]. The distance
										// between b2 and a1 is abs(position_of_b2 - position_of_a1)
										// Let position_of_b2 be (xIndex+subCellShift[0][p]) = (1+0.25) = 1.25.
										// Then position_of_a1 is subCellShift[0][k] = -0.25
										// ==> distance = abs(position_of_b2 - position_of_a1) = abs(1.25 - (- 0.25))
										//              = 1.5
										// From this value, 0.5 has to be substracted (but only if the distance is not
										// zero, so if e.g. xIndex is not null) as 1.5 is the distance between
										// the to cell centers. The distance of the center to the border is 0.25 for
										// subCells. For the "real" distance, the resulting value (here 1.0) has to be
										// multiplied with the cellLength
										xDistanceSquare = pow((abs(xIndex - subCellShifts[0][k] + subCellShifts[0][p]) - 0.5 * xCoarseCellDistanceNotZero) * _cellLength[0], 2);
										yDistanceSquare = pow((abs(yIndex - subCellShifts[1][k] + subCellShifts[1][p]) - 0.5 * yCoarseCellDistanceNotZero) * _cellLength[1], 2);
										zDistanceSquare = pow((abs(zIndex - subCellShifts[2][k] + subCellShifts[2][p]) - 0.5 * zCoarseCellDistanceNotZero) * _cellLength[2], 2);
										// check whether the distance of the cells is smaller than the cutoff radius
										if (xDistanceSquare + yDistanceSquare + zDistanceSquare < cutoffRadiusSquare) {
											// You need to add the local index p of the neighbouring subCell and to
											// substract the local index k of the regarded subCell to compute the right
											// subCell Index Offset
											SubCellOffset = metaCellOffset + p - k;
											if (SubCellOffset > 0) {
												_forwardNeighbourSubOffsets[_metaCellIndex[i] + k].push_back(SubCellOffset);
												assert(SubCellOffset > 0);
												if ((unsigned) SubCellOffset > _maxNeighbourOffset) {
													_maxNeighbourOffset = SubCellOffset;
												}
											}
											else { // 0 can't happen as only real neighbours are considered
												_backwardNeighbourSubOffsets[_metaCellIndex[i] + k].push_back(SubCellOffset);
												if ((unsigned) abs(SubCellOffset) > _minNeighbourOffset) {
													_minNeighbourOffset = abs(SubCellOffset);
												}
											}
										} // if Distance < cutoffRadius
									} // loop over the subcells in the neighbouring cell
								} // CASE 1
								// ########################################
								// ## CASE 2: CELL A REFINED, CELL B NOT ##
								// ########################################
								else {
									// As in case 1, first the distance between the centers is calculated.
									// Only Cell A (loop index k) is shifted, Cell B is not as it is a coarse cell.
									// One of the cells is coarse (distance to boundary 0.5, one is fine (distance 0.25),
									// so for the minimal distance, 0.75 has to be substracted if the two cells don't have
									// the same x (or y or z) -value (See comment of for the variable xCoarseCellDistanceNotZero).
									xDistanceSquare = pow((abs(xIndex - subCellShifts[0][k]) - 0.75 * xCoarseCellDistanceNotZero) * _cellLength[0], 2);
									yDistanceSquare = pow((abs(yIndex - subCellShifts[1][k]) - 0.75 * yCoarseCellDistanceNotZero) * _cellLength[1], 2);
									zDistanceSquare = pow((abs(zIndex - subCellShifts[2][k]) - 0.75 * zCoarseCellDistanceNotZero) * _cellLength[2], 2);
									if (xDistanceSquare + yDistanceSquare + zDistanceSquare <= cutoffRadiusSquare) {
										// You need to substract the local index k of the regarded subCell in order to compute the subCellOffset
										SubCellOffset = metaCellOffset - k;
										if (SubCellOffset > 0) {
											_forwardNeighbourSubOffsets[_metaCellIndex[i] + k].push_back(SubCellOffset);
											assert(SubCellOffset > 0);
											if ((unsigned) SubCellOffset > _maxNeighbourOffset) {
												_maxNeighbourOffset = SubCellOffset;
											}
										}
										else { // 0 can't happen as only real neighbours are considered
											_backwardNeighbourSubOffsets[_metaCellIndex[i] + k].push_back(SubCellOffset);
											if ((unsigned) abs(SubCellOffset) > _minNeighbourOffset) {
												_minNeighbourOffset = abs(SubCellOffset);
											}
										}
									} // if Distance < cutoffRadius
								} // CASE 2
							} // end if for checking if the cell is a boundary cell (to be removed)
						}// loop over xIndex
					} // loop over yIndex
				} // loop over zIndex
			} // loop over the subcells of cell A
		} // if cell A is subcell
		// the easier case: There are no subCells in the current coarse Cell
		// (but may be there are subCells in the neighbouring Cells!)
		else {
			// loop over all neighbouring ("continue" in case of the cell itself) coarse cells
			for (int zIndex = -_haloWidthInNumCells[2]; zIndex <= _haloWidthInNumCells[2]; zIndex++) {
				for (int yIndex = -_haloWidthInNumCells[1]; yIndex <= _haloWidthInNumCells[1]; yIndex++) {
					for (int xIndex = -_haloWidthInNumCells[0]; xIndex <= _haloWidthInNumCells[0]; xIndex++) {
						// The cell itself should not be searched for neighbours!
						if (xIndex == 0 && yIndex == 0 && zIndex == 0)
							continue;
						// compute the global coarseCellOffset for the current coarse grid neighbour cell
						coarseCellOffset = (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
						// _cells.size()-1 is the upper cellIndex boundary! So i + coarseCellOffset must be smaller than _cells.size()
						// @todo only calculate neighbour offets for necessary cells, so this  if  could then be removed
						if ((int) i + coarseCellOffset < (int) _cells.size() && (int) i + coarseCellOffset >= 0) {
							// check for each dimension if the cells have the same coordinate value in that dimension
							if (xIndex == 0)
								xCoarseCellDistanceNotZero = 0.0;
							else
								xCoarseCellDistanceNotZero = 1.0;
							if (yIndex == 0)
								yCoarseCellDistanceNotZero = 0.0;
							else
								yCoarseCellDistanceNotZero = 1.0;
							if (zIndex == 0)
								zCoarseCellDistanceNotZero = 0.0;
							else
								zCoarseCellDistanceNotZero = 1.0;

							// Index Offset of the two coarse cells in the array _subCells
							metaCellOffset = _metaCellIndex[i + coarseCellOffset] - _metaCellIndex[i];
							// ###########################################
							// ## CASE 3: CELL A COARSE, CELL B REFINED ##
							// ###########################################
							if (_metaCellIndex[i + coarseCellOffset + 1] - _metaCellIndex[i + coarseCellOffset] == 8) {
								// loop over all subCells that the neighbour Cell contains
								for (int p = 0; p < 8; p++) {
									xDistanceSquare = pow((abs(xIndex + subCellShifts[0][p]) - 0.75 * xCoarseCellDistanceNotZero) * _cellLength[0], 2);
									yDistanceSquare = pow((abs(yIndex + subCellShifts[1][p]) - 0.75 * yCoarseCellDistanceNotZero) * _cellLength[1], 2);
									zDistanceSquare = pow((abs(zIndex + subCellShifts[2][p]) - 0.75 * zCoarseCellDistanceNotZero) * _cellLength[2], 2);
									if (xDistanceSquare + yDistanceSquare + zDistanceSquare <= cutoffRadiusSquare) {
										// You need to add the local index p of the neighbour subCell in order to compute the subCellOffset
										SubCellOffset = metaCellOffset + p;
										if (SubCellOffset > 0) {
											_forwardNeighbourSubOffsets[_metaCellIndex[i]].push_back(SubCellOffset);
											assert(SubCellOffset > 0);
											if ((unsigned) SubCellOffset > _maxNeighbourOffset) {
												_maxNeighbourOffset = SubCellOffset;
											}
										}
										else { // 0 can't happen as only real neighbours are considered
											_backwardNeighbourSubOffsets[_metaCellIndex[i]].push_back(SubCellOffset);
											if ((unsigned) abs(SubCellOffset) > _minNeighbourOffset) {
												_minNeighbourOffset = abs(SubCellOffset);
											}
										}
									} // if Distance < cutoffRadius
								} // loop over the subcells in the neighbouring cell
							} // CASE 3
							// #################################
							// ## CASE 4: BOTH CELLS COARSE   ##
							// #################################
							else {
								// in the "else else" case (neither the current Cell nor the neighbour Cell contain subCells)
								// the distance in one dimension is the width of a coarse Cell multiplied with the number
								// of coarse cells between the two cells (this is received by substracting one of the
								// absolute difference of the cells, if this difference is not zero)
								xDistanceSquare = pow((abs(xIndex) - 1 * xCoarseCellDistanceNotZero) * _cellLength[0], 2);
								yDistanceSquare = pow((abs(yIndex) - 1 * yCoarseCellDistanceNotZero) * _cellLength[1], 2);
								zDistanceSquare = pow((abs(zIndex) - 1 * zCoarseCellDistanceNotZero) * _cellLength[2], 2);
								if (xDistanceSquare + yDistanceSquare + zDistanceSquare <= cutoffRadiusSquare) {
									// subCellOffset==metaCellOffset, because there are neither local subCells nor subCells
									// within the current neighbour Cell. So we don't need any additional index offsets
									SubCellOffset = metaCellOffset;
									if (SubCellOffset > 0) {
										_forwardNeighbourSubOffsets[_metaCellIndex[i]].push_back(SubCellOffset);
										assert(SubCellOffset > 0);
										if ((unsigned) SubCellOffset > _maxNeighbourOffset) {
											_maxNeighbourOffset = SubCellOffset;
										}
									}
									else { // 0 can't happen as only real neighbours are considered
										_backwardNeighbourSubOffsets[_metaCellIndex[i]].push_back(SubCellOffset);
										if ((unsigned) abs(SubCellOffset) > _minNeighbourOffset) {
											_minNeighbourOffset = abs(SubCellOffset);
										}
									}
								}
							} // CASE 4
						} // end if for checking if the cell is a boundary cell (to be removed)
					}// loop over xIndex
				} // loop over yIndex
			} // loop over zIndex
		} // else case (cell A is not refined)
	} // loop over all cells
#ifndef NDEBUG
	global_log->info() << "Neighbour offsets are bounded by "
			<< _minNeighbourOffset << ", " << _maxNeighbourOffset << endl;
#endif
} // end method


unsigned long AdaptiveSubCells::getCellIndexOfMolecule(Molecule* molecule) {
	int cellIndex[3]; // 3D Cell index

	for (int dim = 0; dim < 3; dim++) {
		if (molecule->r(dim) < _haloBoundingBoxMin[dim] || molecule->r(dim) >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "AdaptiveSubCells::getCellIndexOfMolecule(Molecule* molecule): Molecule is outside of the bounding box" << endl;
		}
		cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);

	}
	return (cellIndex[2] * _cellsPerDimension[1] + cellIndex[1]) * _cellsPerDimension[0] + cellIndex[0];
}

unsigned long AdaptiveSubCells::getSubCellIndexOfMolecule(Molecule* molecule) {
	// 3D Cell index
	int cellIndex[3], globalCellIndex, metaSubCellIndex, subCell3DIndexOffset[3], subCellIndexOffset;

	for (int dim = 0; dim < 3; dim++) {
		if (molecule->r(dim) < _haloBoundingBoxMin[dim] || molecule->r(dim) >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "AdaptiveSubCells::getSubCellIndexOfMolecule(Molecule* molecule): Molecule is outside of the bounding box" << endl;
		}
		cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);

	}
	// globalCellIndex is the cellIndex of the molecule on the coarse Cell level.
	// It needs to be transformed to the subCell level
	globalCellIndex = (cellIndex[2] * _cellsPerDimension[1] + cellIndex[1]) * _cellsPerDimension[0] + cellIndex[0];
	// metaSubCellIndex is the index of the subCell at the lower, left, front end of the current coarse Cell.
	// For adressing the other local subCells an index offset needs to be added.
	metaSubCellIndex = _metaCellIndex[globalCellIndex];
	// if (_metaCellIndex[r+1]-_metaCellIndex[r]==8) the Cell with global Index r contains subCells
	if (_metaCellIndex[globalCellIndex + 1] - _metaCellIndex[globalCellIndex] == 8) {
		// computes an Index Offset for each spacial dimension (x,y,z).
		for (int d = 0; d < 3; d++) {
			// the case matching is needed to decide for each spacial dimension (x,y,z) whether the
			// molecule must be inserted into the first or into the last subCell.
			if (((molecule->r(d) - _haloBoundingBoxMin[d]) / _cellLength[d]) - floor((molecule->r(d) - _haloBoundingBoxMin[d]) / _cellLength[d]) >= 0.5) {

				// insert molecule into the last subCell of dimension d, d=0,1,2
				subCell3DIndexOffset[d] = 1;
			}
			else {
				// insert molecule into the first subCell of dimension d, d=0,1,2
				subCell3DIndexOffset[d] = 0;
			}
		}
		// transforms the 3D index offset of the subCell into the global subCellIndexOffset:
		subCellIndexOffset = 4 * subCell3DIndexOffset[2] + 2 * subCell3DIndexOffset[1] + subCell3DIndexOffset[0];
		// add the subCellIndexOffset and return the resulting index of the subCell where the molecule has to be inserted
		return (metaSubCellIndex + subCellIndexOffset);
	}
	else {
		// if the current Cell doesn't contain subCells the molecule's globalSubCellIndex is idem with its metaSubCellIndex
		return metaSubCellIndex;
	}
}

unsigned long AdaptiveSubCells::subCellIndexOf3DIndex(int xIndex, int yIndex, int zIndex) {

	int globalCellIndex, metaSubCellIndex;
	// globalCellIndex is the cellIndex of the molecule on the Cell level. It needs to be transformed to the subCell level
	globalCellIndex = (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
	// this case matching prevents the method from returning indices that are overstepping the upper index boundary of the array
	if (globalCellIndex >= (int) _cells.size()) {
		metaSubCellIndex = _subCells.size() - 1;
	}

	else {
		// metaSubCellIndex is the index of the subCell at the lower, left, front end of the current coarse Cell
		metaSubCellIndex = _metaCellIndex[globalCellIndex];
	}
	// returns the metaSubCellIndex. A suitable index offset needs to be added
	//in the calling function in order to adress all subCells
	return metaSubCellIndex;
}

void AdaptiveSubCells::calculateMetaCellIndex() {

	int metaCellIter = 0;
	// calls the computation of _localRho[]
	// may be you can integrate this method as inline method
	calculateLocalRho();
	// loop over all cellIndices in the coarse grid
	for (int r = 0; r < (int) _cells.size(); r++) {
		_metaCellIndex[r] = metaCellIter;
		// if localRho >= rhoLimit you need to increment the cellIter 8 times in order to create space for the subCell index offsets. For example rhoLimit=44 is the minimum density of molecules per Cell in a standard fluid
		if (_localRho[r] >= 16) {
			metaCellIter += 7;
		}
		// if _localRho < rhoLimit you need to increment the cellIter only once
		metaCellIter++;
	}
	//compute the numberOfSubCells (=the maximum subCellIndex + 1) in order to be able do determine the size of the subCell array
	int numberOfSubCells = metaCellIter;
	// insert an upper index boundary into the array in order to be able to handle the getSubCellIndexOfMolecule method
	_metaCellIndex[_cells.size()] = numberOfSubCells;
}

void AdaptiveSubCells::calculateLocalRho() {
	// may be you can get along without the integer variables cellIntIterator and minIterator
	int cellIntIterator = 0;
	// loop over all cells
	vector<ParticleCell>::iterator cellIterator;
	std::vector<Molecule*>::iterator molIterator;
	// traverse all coarse cells
	for (cellIterator = _cells.begin(); cellIterator != _cells.end(); cellIterator++) {
		int molIntIterator = 0;
		for (molIterator = cellIterator->getParticlePointers().begin();
		// traverse all molecules in the current coarse cell
		molIterator != cellIterator->getParticlePointers().end(); molIterator++) {
			// count the molecules that the current coarse cell contains
			molIntIterator++;
		}
		_localRho[cellIntIterator] = molIntIterator;
		// cellIntIterator is the integer value of cellIterator
		cellIntIterator++;
	}
}


int AdaptiveSubCells::grandcanonicalBalance(DomainDecompBase* comm) {
	comm->collCommInit(1);
	comm->collCommAppendInt(_localInsertionsMinusDeletions);
	comm->collCommAllreduceSum();
	int universalInsertionsMinusDeletions = comm->collCommGetInt();
	comm->collCommFinalize();
	return universalInsertionsMinusDeletions;
}

void AdaptiveSubCells::grandcanonicalStep(ChemicalPotential* mu, double T, Domain* domain, CellProcessor& cellProcessor) {
	bool accept = true;
	double DeltaUpot;
	Molecule* m;
	ParticlePairs2PotForceAdapter particlePairsHandler(*domain);

	_localInsertionsMinusDeletions = 0;

	mu->submitTemperature(T);
	double minco[3];
	double maxco[3];
	for (int d = 0; d < 3; d++) {
		minco[d] = this->getBoundingBoxMin(d);
		maxco[d] = this->getBoundingBoxMax(d);
	}

	bool hasDeletion = true;
	bool hasInsertion = true;
	double ins[3];
	unsigned nextid = 0;
	while (hasDeletion || hasInsertion) {
		if (hasDeletion)
			hasDeletion = mu->getDeletion(this, minco, maxco);
		if (hasDeletion) {
			m = &(*(_particleIter));
			DeltaUpot = -1.0 * getEnergy(&particlePairsHandler, m, cellProcessor);

			accept = mu->decideDeletion(DeltaUpot / T);
#ifndef NDEBUG
			if(accept) global_log->debug() << "r" << mu->rank() << "d" << m->id() << endl;
			else global_log->debug() << "   (r" << mu->rank() << "-d" << m->id() << ")" << endl;
#endif
			if (accept) {
				m->upd_cache();
				// reset forces and momenta to zero
				{
					double zeroVec[3] = {0.0, 0.0, 0.0};
					m->setF(zeroVec);
					m->setM(zeroVec);
				}
				mu->storeMolecule(*m);
				deleteMolecule(m->id(), m->r(0), m->r(1), m->r(2));
				_particles.erase(_particleIter);
				_particleIter = _particles.begin();
				_localInsertionsMinusDeletions--;
			}
		}

		if (mu->isWidom()){
			m = &(*(_particles.begin()));
			mu->storeMolecule(*m);
		}
		if (hasInsertion) {
			nextid = mu->getInsertion(ins);
			hasInsertion = (nextid > 0);
		}
		if (hasInsertion) {
			// for(int d = 0; d < 3; d++)
			//    ins[d] = ins[d]-coords[d]*proc_domain_L[d]-m_rmin[d];
			Molecule tmp = mu->loadMolecule();
			for (int d = 0; d < 3; d++)
				tmp.setr(d, ins[d]);
			tmp.setid(nextid);
			this->_particles.push_back(tmp);

			std::list<Molecule>::iterator mit = _particles.end();
			mit--;
			m = &(*mit);
			m->upd_cache();
			// reset forces and momenta to zero
			if(!mu->isWidom()) {
				double zeroVec[3] = {0.0, 0.0, 0.0};
				m->setF(zeroVec);
				m->setM(zeroVec);
			}
			m->check(nextid);
#ifndef NDEBUG
			global_log->debug() << "rank " << mu->rank() << ": insert " << m->id()
			<< " at the reduced position (" << ins[0] << "/" << ins[1] << "/" << ins[2] << ")? " << endl;
#endif

			unsigned long cellid = this->getCellIndexOfMolecule(m);
			_cells[cellid].addParticle(m);
			DeltaUpot = getEnergy(&particlePairsHandler, m, cellProcessor);
                        domain->submitDU(mu->getComponentID(), DeltaUpot, ins);
			accept = mu->decideInsertion(DeltaUpot / T);

#ifndef NDEBUG
			if(accept) global_log->debug() << "r" << mu->rank() << "i" << mit->id() << ")" << endl;
			else global_log->debug() << "   (r" << mu->rank() << "-i" << mit->id() << ")" << endl;
#endif
			if (accept) {
				_localInsertionsMinusDeletions++;
			}
			else {
			//	deleteMolecule(m->id(), m->r(0), m->r(1), m->r(2));
				this->_cells[cellid].deleteMolecule(m->id());

				mit->check(m->id());
				_particles.erase(mit);
			}
		}
	}
	for (m = this->begin(); m != this->end(); m = this->next()) {
#ifndef NDEBUG
		m->check(m->id());
#endif
	}
}
