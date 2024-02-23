/*
 * MoleculeInterface.cpp
 *
 *  Created on: 21 Jan 2017
 *      Author: tchipevn
 */

#include "MoleculeInterface.h"

#include "utils/mardyn_assert.h"
#include <cmath>
#include <fstream>

#include "utils/Logger.h"
#include "Simulation.h"


bool MoleculeInterface::isLessThan(const MoleculeInterface& m2) const {
	if (r(2) < m2.r(2))
		return true;
	else if (r(2) > m2.r(2))
		return false;
	else {
		if (r(1) < m2.r(1))
			return true;
		else if (r(1) > m2.r(1))
			return false;
		else {
			if (r(0) < m2.r(0))
				return true;
			else if (r(0) > m2.r(0))
				return false;
			else {
				Log::global_log->error() << "LinkedCells::isFirstParticle: both Particles have the same position" << std::endl;
				Simulation::exit(1);
			}
		}
	}
	return false; /* Silence warnings about missing return statement */
}

