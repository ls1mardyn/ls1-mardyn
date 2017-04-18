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

void ExplicitEuler_WR::computePositions(ParticleContainer* molCont, Domain* dom) {
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

void ExplicitEuler_WR::computeVelocities(ParticleContainer* molCont, Domain* dom) {
	// TODO: I hate this piece of code. It will be rewritten from within the Leapfrog integrator ASAP.
	map<int, unsigned long> N;
	map<int, unsigned long> rotDOF;
	map<int, double> summv2;
	map<int, double> sumIw2;
	{
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			unsigned long Ngt_l = 0;
			unsigned long rotDOFgt_l = 0;
			double summv2gt_l = 0.0;
			double sumIw2gt_l = 0.0;

			const ParticleIterator begin = molCont->iteratorBegin();
			const ParticleIterator end = molCont->iteratorEnd();

			for (ParticleIterator i = begin; i != end; ++i) {
				double dummy = 0.0;
				i->ee_upd_postF(_timestepLength, summv2gt_l);
				mardyn_assert(summv2gt_l >= 0.0);
				Ngt_l++;
				rotDOFgt_l += i->component()->getRotationalDegreesOfFreedom();
			}

			#if defined(_OPENMP)
			#pragma omp critical (thermostat)
			#endif
			{
				N[0] += Ngt_l;
				rotDOF[0] += rotDOFgt_l;
				summv2[0] += summv2gt_l;
				sumIw2[0] += sumIw2gt_l;
			}
		} // end pragma omp parallel
	}
	for (map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++) {
		dom->setLocalSummv2(thermit->second, thermit->first);
		dom->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
		dom->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
	}
}
