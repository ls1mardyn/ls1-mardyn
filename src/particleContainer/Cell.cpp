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

list<Molecule*>& Cell::getParticlePointers() {
	return this->_particlePointers;
}

void Cell::assingCellToHaloRegion() {
	this->_haloCellState = true;
}

void Cell::assignCellToBoundaryRegion() {
	this->_boundaryCellState = true;
}

void Cell::assignCellToInnerRegion() {
	this->_innerCellState = true;
}

bool Cell::isHaloCell() {
	return _haloCellState;
}

bool Cell::isBoundaryCell() {
	return _boundaryCellState;
}

bool Cell::isInnerCell() {
	return _innerCellState;
}

int Cell::getMoleculeCount() const {
	return _particlePointers.size();
}

bool Cell::deleteMolecule(unsigned long molid) {
	bool found = false;
	list<Molecule*>::iterator cellit;

	for (cellit = this->_particlePointers.begin(); cellit != this->_particlePointers.end(); cellit++) {
		if ((*cellit)->id() == molid) {
			found = true;
			this->_particlePointers.remove(*cellit);
			break;
		}
	}
	return found;
}
