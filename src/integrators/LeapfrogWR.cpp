/*
 * LeapFrogWR.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: tchipevn
 */

#include "LeapfrogWR.h"

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

Leapfrog_WR::Leapfrog_WR(double timestepLength) :
		Integrator(timestepLength) {
}

void Leapfrog_WR::readXML(XMLfileUnits & xmlconfig) {
	_timestepLength = 0;
	xmlconfig.getNodeValueReduced("timestep", _timestepLength);
	global_log->info() << "Timestep: " << _timestepLength << endl;
	mardyn_assert(_timestepLength > 0);
}

void Leapfrog_WR::computePositions(ParticleContainer* molCont, Domain* dom) {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const ParticleIterator begin = molCont->iteratorBegin();
		const ParticleIterator end = molCont->iteratorEnd();
		for(auto i = begin; i != end; ++i) {
			i->ee_upd_preF(_timestepLength);
		}
	}
}

void Leapfrog_WR::computeVelocities(ParticleContainer* molCont, Domain* dom) {
	// TODO: Thermostat functionality is duplicated X times and needs to be rewritten!
	map<int, unsigned long> N;
	map<int, unsigned long> rotDOF;
	map<int, double> summv2;
	map<int, double> sumIw2;
	{
		unsigned long red_N = 0;
		unsigned long red_rotDOF = 0;
		double red_summv2 = 0.0;
		double red_sumIw2 = 0.0;
		#if defined(_OPENMP)
		#pragma omp parallel reduction(+: red_N, red_rotDOF, red_summv2, red_sumIw2)
		#endif
		{
			const ParticleIterator begin = molCont->iteratorBegin();
			const ParticleIterator end = molCont->iteratorEnd();

			for (ParticleIterator i = begin; i != end; ++i) {
				double dummy = 0.0;
				i->ee_upd_postF(_timestepLength, red_summv2);
				mardyn_assert(red_summv2 >= 0.0);
				red_N++;
				red_rotDOF += i->component()->getRotationalDegreesOfFreedom();
			}
		} // end pragma omp parallel
		N[0] += red_N;
		rotDOF[0] += red_rotDOF;
		summv2[0] += red_summv2;
		sumIw2[0] += red_sumIw2;
	}
	for (map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++) {
		dom->setLocalSummv2(thermit->second, thermit->first);
		dom->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
		dom->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
	}
}
