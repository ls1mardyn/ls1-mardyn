/*
 * MoleculeRMM.cpp
 *
 *  Created on: 21 Jan 2017
 *      Author: tchipevn
 */

#include "MoleculeRMM.h"
#include "Simulation.h"
#include "molecules/Component.h"
#include "molecules/Quaternion.h"
#include "ensemble/EnsembleBase.h"
#include "particleContainer/adapter/CellDataSoARMM.h"


bool 			MoleculeRMM::_initCalled = false;
Component * 	MoleculeRMM::_component;
Quaternion 		MoleculeRMM::_quaternion;

void MoleculeRMM::initStaticVars() {
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

void MoleculeRMM::setSoA(CellDataSoABase * const s) {
	CellDataSoARMM * derived;
	derived = static_cast<CellDataSoARMM *>(s);
	_soa = derived;
}

double MoleculeRMM::r(unsigned short d) const {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		return _r[d];
	} else {
		return _soa->getMolR(d, _soa_index);
	}
}

double MoleculeRMM::v(unsigned short d) const {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		return _v[d];
	} else {
		return _soa->getMolV(d, _soa_index);
	}
}

unsigned long MoleculeRMM::getID() const {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		return _id;
	} else {
		return _soa->getMolUid(_soa_index);
	}
}

void MoleculeRMM::setr(unsigned short d, double r) {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		_r[d] = r;
	} else {
		_soa->setMolR(d, _soa_index, r);
	}
}

void MoleculeRMM::setv(unsigned short d, double v) {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		_v[d] = v;
	} else {
		_soa->setMolV(d, _soa_index, v);
	}
}

void MoleculeRMM::setid(unsigned long id) {
	mardyn_assert(_state == STORAGE_SOA or _state == STORAGE_AOS);

	if (_state == STORAGE_AOS) {
		_id = id;
	} else {
		_soa->setMolUid(_soa_index, id);
	}
}

void MoleculeRMM::setF(unsigned short /*d*/, double /*F*/) {

}

std::string MoleculeRMM::getWriteFormat(){
	return std::string("IRV");
}

void MoleculeRMM::write(std::ostream& ostrm) const {
	ostrm << getID() << "\t"
		  << r(0) << " " << r(1) << " " << r(2) << "\t"
		  << v(0) << " " << v(1) << " " << v(2) << "\t"
		  << "\n";
}

std::ostream& operator<<( std::ostream& os, const MoleculeRMM& m ) {
	os << "ID: " << m.getID() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	return os;
}
