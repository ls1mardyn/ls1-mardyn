/*
 * LeapFrogWR.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: tchipevn
 */

#include "ExplicitEulerWR.h"

#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "ensemble/PressureGradient.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/ParticleIterator.h"

using namespace std;
using Log::global_log;

ExplicitEuler_WR::ExplicitEuler_WR(double timestepLength) :
		Integrator(timestepLength) {
}

void ExplicitEuler_WR::readXML(XMLfileUnits & xmlconfig) {
	_timestepLength = 0;
	xmlconfig.getNodeValueReduced("timestep", _timestepLength);
	global_log->info() << "Timestep: " << _timestepLength << endl;
	mardyn_assert(_timestepLength > 0);
}

void ExplicitEuler_WR::eventNewTimestep(ParticleContainer* moleculeContainer,
		Domain* domain) {
	computePositions(moleculeContainer, domain);
}

void ExplicitEuler_WR::computePositions(ParticleContainer* molCont, Domain* dom) {

}
