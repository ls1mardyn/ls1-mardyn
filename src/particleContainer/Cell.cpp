#include "particleContainer/Cell.h"
#include "molecules/Molecule.h"

using namespace std;

Cell::Cell() {
	this->_haloCellState = false;
	this->_boundaryCellState = false;
	this->_innerCellState = false;
}

void Cell::removeAllParticles() {
	this->_particlePointers.clear();
}

void Cell::addParticle(Molecule* particle_ptr) {
	_particlePointers.push_back(particle_ptr);
}

vector<Molecule*>& Cell::getParticlePointers() {
	return this->_particlePointers;
}

void Cell::assignCellToHaloRegion() {
	this->_haloCellState = true;
}

void Cell::assignCellToBoundaryRegion() {
	this->_boundaryCellState = true;
}

void Cell::assignCellToInnerRegion() {
	this->_innerCellState = true;
}

bool Cell::isHaloCell() const {
	return _haloCellState;
}

bool Cell::isBoundaryCell() const {
	return _boundaryCellState;
}

bool Cell::isInnerCell() const {
	return _innerCellState;
}

int Cell::getMoleculeCount() const {
	return _particlePointers.size();
}

bool Cell::deleteMolecule(unsigned long molid) {
	bool found = false;
	vector<Molecule*>::iterator cellit;

	for (cellit = _particlePointers.begin(); cellit != _particlePointers.end(); cellit++) {
		if ((*cellit)->id() == molid) {
			found = true;
			_particlePointers.erase(cellit);
			break;
		}
	}
	return found;
}
