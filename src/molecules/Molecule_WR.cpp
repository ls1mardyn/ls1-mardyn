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
