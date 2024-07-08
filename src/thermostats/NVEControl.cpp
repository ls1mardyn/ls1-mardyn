#include "thermostats/NVEControl.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"


NVEControl::NVEControl() {
	Log::global_log->info() << "[NVEControl] initialized."<<std::endl;
}

NVEControl::~NVEControl() {
}

void NVEControl::readXML(XMLfileUnits& xmlconfig) {
	Log::global_log->info() << "[NVEControl] called." << std::endl;
	xmlconfig.getNodeValue("control/start", this->_startTime);
	xmlconfig.getNodeValue("control/stop", this->_endTime);
	xmlconfig.getNodeValue("control/frequency", this->_frequency);
	xmlconfig.getNodeValue("control/writefreq", this->_writeFrequency);
	xmlconfig.getNodeValue("control/fileprefix", this->_fileprefix);
	Log::global_log->info() << "[NVEControl] xml read." << std::endl;
}

void NVEControl::init(){
	this->_globalBetaTrans = 1.0;
	this->_Nglobal = _simulation.getDomain()->getglobalNumMolecules(false, nullptr, nullptr);
	
	double globalT = _simulation.getDomain()->getGlobalCurrentTemperature();
	double globalUpot = _simulation.getDomain()->getGlobalUpot();
	
	this->_Utarget = globalUpot + 3./2. * _Nglobal * globalT;

	Log::global_log->info() << "[NVEControl] _Nglobal = " << this->_Nglobal <<std::endl;
	Log::global_log->info() << "[NVEControl] globalT = " << globalT <<std::endl;
	Log::global_log->info() << "[NVEControl] globalUpot = " << globalUpot <<std::endl;
	Log::global_log->info() << "[NVEControl] _Utarget = " << this->_Utarget <<std::endl;

	// Log::global_log->info() << "[NVEControl] initialized. _Nglobal = " << this->_Nglobal << ", _Utarget = "<< this->_Utarget<<std::endl;
	return;
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
		double Tscalingfactor = (_Utarget - globalUpot)/(globalT * 3/2 * _Nglobal);
		this->_globalBetaTrans = std::sqrt(Tscalingfactor);
		// this->betaRot = vscalingfactor;

		// Log::global_log->info() << "[NVEControl] _Utarget = " << _Utarget << std::endl;
		// Log::global_log->info() << "[NVEControl] globalT = " << globalT << std::endl;
		// Log::global_log->info() << "[NVEControl] globalUpot = " << globalUpot << std::endl;
		// Log::global_log->info() << "[NVEControl] TScalingfactor = " << Tscalingfactor << std::endl;
		// Log::global_log->info() << "[NVEControl] _globalBetaTrans = " << _globalBetaTrans << std::endl;

		Log::global_log->info() << "[NVEControl] _Utarget : TScalingfactor : _globalBetaTrans = " << _Utarget << ":" << Tscalingfactor << ":" <<_globalBetaTrans << std::endl;

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
