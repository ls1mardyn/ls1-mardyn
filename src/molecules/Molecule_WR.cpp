/*
 * Molecule_WR.cpp
 *
 *  Created on: 21 Jan 2017
 *      Author: tchipevn
 */

#include "Molecule_WR.h"
#include "Simulation.h"
#include "molecules/Component.h"
#include "molecules/Quaternion.h"
#include "ensemble/EnsembleBase.h"
#include "particleContainer/adapter/CellDataSoA_WR.h"


bool 			Molecule_WR::_initCalled = false;
Component * 	Molecule_WR::_component;
Quaternion 		Molecule_WR::_quaternion;

void Molecule_WR::initStaticVars() {
	if (not _initCalled) {
		_component = global_simulation->getEnsemble()->getComponent(0);
		_quaternion = Quaternion(1.0, 0.0, 0.0, 0.0);
		_initCalled = true;
	}
}

void Molecule_WR::setSoA(CellDataSoABase * const s) {
	assert(_state == AOS);
	CellDataSoA_WR * derived;
#ifndef NDEBUG
	derived = nullptr;
	derived = dynamic_cast<CellDataSoA_WR *>(s);
	if(derived == nullptr and s != nullptr) {
		global_log->error() << "expected CellDataSoA_WR pointer for m" << _id << std::endl;
		assert(false);
	}
#else
	derived = static_cast<CellDataSoA_WR *>(s);
#endif
	_soa = derived;
}

double Molecule_WR::r(unsigned short d) const {
	assert(_state == SOA or _state == AOS);

	if (_state == AOS) {
		return _r[d];
	} else {
		size_t linOffset = _soa->_mol_r.dimensionToOffset(d);
		return _soa->_mol_r.linearCrossAccess(linOffset + _soa_index);
	}
}

double Molecule_WR::v(unsigned short d) const {
	assert(_state == SOA or _state == AOS);

	if (_state == AOS) {
		return _v[d];
	} else {
		size_t linOffset = _soa->_mol_v.dimensionToOffset(d);
		return _soa->_mol_v.linearCrossAccess(linOffset + _soa_index);
	}
}

unsigned long Molecule_WR::id() const {
	assert(_state == SOA or _state == AOS);

	if (_state == AOS) {
		return _id;
	} else {
		return _soa->_mol_uid[_soa_index];
	}
}

void Molecule_WR::setr(unsigned short d, double r) {
	assert(_state == SOA or _state == AOS);

	if (_state == AOS) {
		_r[d] = r;
	} else {
		size_t linOffset = _soa->_mol_r.dimensionToOffset(d);
		_soa->_mol_r.linearCrossAccess(linOffset + _soa_index) = r;
	}
}

void Molecule_WR::setv(unsigned short d, double v) {
	assert(_state == SOA or _state == AOS);

	if (_state == AOS) {
		_v[d] = v;
	} else {
		size_t linOffset = _soa->_mol_r.dimensionToOffset(d);
		_soa->_mol_v.linearCrossAccess(linOffset + _soa_index) = v;
	}
}

void Molecule_WR::setid(unsigned long id) {
	assert(_state == SOA or _state == AOS);

	if (_state == AOS) {
		_id = id;
	} else {
		_soa->_mol_uid[_soa_index] = id;
	}
}

std::string Molecule_WR::getWriteFormat(){
	// TODO?
	assert(false);
	return std::string("ICRV");
}

std::ostream& operator<<( std::ostream& os, const Molecule_WR& m ) {
	os << "ID: " << m.id() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	return os;
}
