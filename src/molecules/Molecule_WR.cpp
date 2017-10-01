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
		if (global_simulation == nullptr)
			return;
		if (global_simulation->getEnsemble() == nullptr)
			return;
		if (global_simulation->getEnsemble()->getComponents() == nullptr){
			return;
		}
		if (global_simulation->getEnsemble()->getComponents()->size() == 0) {
			return;
		}
		_component = global_simulation->getEnsemble()->getComponent(0);

		if (_component == nullptr) {
			// we are in some constructor and are not ready to initialise yet
			return;
		}
		_quaternion = Quaternion(1.0, 0.0, 0.0, 0.0);
		_initCalled = true;
	}
}

void Molecule_WR::setSoA(CellDataSoABase * const s) {
	CellDataSoA_WR * derived;
#ifndef NDEBUG
	derived = nullptr;
	derived = dynamic_cast<CellDataSoA_WR *>(s);
	if(derived == nullptr and s != nullptr) {
		global_log->error() << "expected CellDataSoA_WR pointer for m" << _id << std::endl;
		mardyn_assert(false);
	}
#else
	derived = static_cast<CellDataSoA_WR *>(s);
#endif
	_soa = derived;
}

double Molecule_WR::r(unsigned short d) const {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		return _r[d];
	} else {
		return _soa->getMolR(d, _soa_index);
	}
}

double Molecule_WR::v(unsigned short d) const {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		return _v[d];
	} else {
		return _soa->getMolV(d, _soa_index);
	}
}

unsigned long Molecule_WR::id() const {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		return _id;
	} else {
		return _soa->getMolUid(_soa_index);
	}
}

void Molecule_WR::setr(unsigned short d, double r) {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		_r[d] = r;
	} else {
		_soa->setMolR(d, _soa_index, r);
	}
}

void Molecule_WR::setv(unsigned short d, double v) {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		_v[d] = v;
	} else {
		_soa->setMolV(d, _soa_index, v);
	}
}

void Molecule_WR::setid(unsigned long id) {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		_id = id;
	} else {
		_soa->setMolUid(_soa_index, id);
	}
}

std::string Molecule_WR::getWriteFormat(){
	return std::string("IRV");
}

void Molecule_WR::write(std::ostream& ostrm) const {
	ostrm << id() << "\t"
		  << r(0) << " " << r(1) << " " << r(2) << "\t"
		  << v(0) << " " << v(1) << " " << v(2) << "\t"
		  << endl;
}

std::ostream& operator<<( std::ostream& os, const Molecule_WR& m ) {
	os << "ID: " << m.id() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	return os;
}
