#include "LeapfrogRMM.h"

#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/ParticleIterator.h"

#include "PositionCellProcessorRMM.h"
#include "VelocityCellProcessorRMM.h"


LeapfrogRMM::LeapfrogRMM() {
	_velocityCellProcessor = nullptr;
}

LeapfrogRMM::~LeapfrogRMM() {
	if (_velocityCellProcessor != nullptr) {
		delete _velocityCellProcessor;
	}
}

LeapfrogRMM::LeapfrogRMM(double timestepLength) :
		Integrator(timestepLength) {
	_velocityCellProcessor = nullptr;
}

void LeapfrogRMM::readXML(XMLfileUnits & xmlconfig) {
	_timestepLength = 0;
	xmlconfig.getNodeValueReduced("timestep", _timestepLength);
	global_log->info() << "Timestep: " << _timestepLength << std::endl;
	mardyn_assert(_timestepLength > 0);

	mardyn_assert(_velocityCellProcessor == nullptr);
	_velocityCellProcessor = new VelocityCellProcessorRMM();
}

void LeapfrogRMM::computePositions(ParticleContainer* molCont, Domain* dom) {
#ifndef ENABLE_REDUCED_MEMORY_MODE
	// leaving old functionality for debugging purposes
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const auto begin = molCont->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		for(auto i = begin; i.isValid(); ++i) {
			i->ee_upd_preF(_timestepLength);
		}
	}
#else
	// this is actually called in RMM
	PositionCellProcessorRMM cellProc(_timestepLength);
	molCont->traverseCells(cellProc);
#endif
}

void LeapfrogRMM::computeVelocities(ParticleContainer* molCont, Domain* dom) {
#ifndef ENABLE_REDUCED_MEMORY_MODE
	// leaving old functionality for debugging purposes

	// TODO: Thermostat functionality is duplicated X times and needs to be rewritten!
	std::map<int, unsigned long> N;
	std::map<int, unsigned long> rotDOF;
	std::map<int, double> summv2;
	std::map<int, double> sumIw2;
	{
		unsigned long red_N = 0;
		unsigned long red_rotDOF = 0;
		double red_summv2 = 0.0;
		double red_sumIw2 = 0.0;
		#if defined(_OPENMP)
		#pragma omp parallel reduction(+: red_N, red_rotDOF, red_summv2, red_sumIw2)
		#endif
		{

			for (auto i = molCont->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); i.isValid(); ++i) {
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
	for (std::map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++) {
		dom->setLocalSummv2(thermit->second, thermit->first);
		dom->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
		dom->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
	}
#else
	molCont->traverseCells(*_velocityCellProcessor);
	unsigned long N = _velocityCellProcessor->getN();
	double summv2 = _velocityCellProcessor->getSummv2();

	dom->setLocalSummv2(summv2, 0);
	dom->setLocalSumIw2(0.0, 0);
	dom->setLocalNrotDOF(0, N, 0);
#endif
}
