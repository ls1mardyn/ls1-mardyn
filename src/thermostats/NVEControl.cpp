#include "thermostats/NVEControl.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"


NVEControl::NVEControl() {
	this->_globalBetaTrans = 1.0;
	this->_Nglobal = _simulation.getDomain()->getglobalNumMolecules(false, nullptr, nullptr);
	
	double globalT = _simulation.getDomain()->getGlobalCurrentTemperature();
	double globalUpot = _simulation.getDomain()->getGlobalUpot();
	
	this->_Utarget = globalUpot + 3./2. * _Nglobal * globalT;
	Log::global_log->info() << "[NVEControl] initialized. _Nglobal = " << this->_Nglobal << ", _Utarget = "<< this->_Utarget;
}

NVEControl::~NVEControl() {
}

void NVEControl::readXML(XMLfileUnits& xmlconfig) {	



	xmlconfig.getNodeValue("control/start", this->_startTime);
	xmlconfig.getNodeValue("control/stop", this->_endTime);
	xmlconfig.getNodeValue("control/frequency", this->_frequency);
	xmlconfig.getNodeValue("control/writefreq", this->_writeFrequency);
	xmlconfig.getNodeValue("control/fileprefix", this->_fileprefix);
	Log::global_log->info() << "[NVEControl] xml read.";
}



void NVEControl::setBetaTrans(double beta) {
	_globalBetaTrans = beta;
}

void NVEControl::apply(ParticleContainer *moleculeContainer) {
	unsigned long simstep = _simulation.getSimulationStep();
	if ((simstep % _frequency) == 0)
	{
		this->_Nglobal =  _simulation.getDomain()->getglobalNumMolecules(false, nullptr, nullptr);
		double globalT = _simulation.getDomain()->getGlobalCurrentTemperature();
		double globalUpot = _simulation.getDomain()->getGlobalUpot();
		double Tscalingfactor = (globalT * 3/2 * _Nglobal)/(_Utarget - globalUpot);
		this->_globalBetaTrans = std::sqrt(Tscalingfactor);
		// this->betaRot = vscalingfactor;

		Log::global_log->info() << "[NVEControl] _Utarget = " << _Utarget;
		Log::global_log->info() << "[NVEControl] globalT = " << globalT;
		Log::global_log->info() << "[NVEControl] globalUpot = " << globalUpot;
		Log::global_log->info() << "[NVEControl] TScalingfactor = " << Tscalingfactor;
		Log::global_log->info() << "[NVEControl] _globalBetaTrans = " << _globalBetaTrans;

		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			double betaTrans = _globalBetaTrans;
			// double betaRot = _globalBetaRot;
			Log::global_log->debug() << "Beta trans: " << betaTrans << std::endl;
			// Log::global_log->debug() << "Beta rot: " << betaRot << std::endl;

			for (auto i = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); i.isValid(); ++i) {
				i->scale_v(betaTrans);
				// i->scale_D(betaRot);
			}
		} // end pragma omp parallel
	}
}
