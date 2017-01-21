/*
 * Molecule_WR.cpp
 *
 *  Created on: 21 Jan 2017
 *      Author: tchipevn
 */

#include "Molecule_WR.h"
#include "Simulation.h"

Component Molecule_WR::_component = global_simulation->getEnsemble()->getComponent(0);
Quaternion Molecule_WR::_quaternion = Quaternion(1.0, 0.0, 0.0, 0.0);
