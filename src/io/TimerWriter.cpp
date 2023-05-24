//
// Created by seckler on 10.10.19.
//

#include <chrono>
#include <iomanip>
#include <ctime>

#include "TimerWriter.h"

#include "Simulation.h"
#include "parallel/DomainDecompBase.h"

void TimerWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << std::endl;

	XMLfile::Query query = xmlconfig.query("timers/timer");
	std::string oldpath = xmlconfig.getcurrentnodepath();
	for (auto timerIter = query.begin(); timerIter; ++timerIter) {
		xmlconfig.changecurrentnode(timerIter);
		std::string timername;
		xmlconfig.getNodeValue("name", timername);
		_timerNames.push_back(timername);

		bool incrementalTimer = false;  // if the timer is increasing in every time step, this should be true
		xmlconfig.getNodeValue("incremental", incrementalTimer);
		_incremental.push_back(incrementalTimer);

		global_log->info() << "Added timer for LB monitoring: " << timername << ", incremental: " << incrementalTimer
						   << std::endl;
	}
	if (_timerNames.empty()) {
		global_log->error() << "TimerWriter: no timers given. make sure you specify them correctly." << std::endl;
		Simulation::exit(242367);
	}
	xmlconfig.changecurrentnode(oldpath);
}

void TimerWriter::init(ParticleContainer* /*particleContainer*/, DomainDecompBase* domainDecomp, Domain* /*domain*/) {
	auto rank = domainDecomp->getRank();
	std::stringstream filename;

	const auto time_tNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

	const auto maxRank = domainDecomp->getNumProcs();
	const auto numDigitsMaxRank = std::to_string(maxRank).length();
	tm unused{};
	filename << _outputPrefix << "-rank" << std::setfill('0') << std::setw(numDigitsMaxRank) << rank << "_"
			 << std::put_time(localtime_r(&time_tNow, &unused), "%Y-%m-%d_%H-%M-%S") << ".dat";

	_fileStream.open(filename.str(), std::ofstream::out);

	_fileStream << "timestep";

	for (const auto& timerName : _timerNames) {
		_fileStream << ";" << timerName;
	}
	_fileStream << std::endl;

	_times.resize(_timerNames.size(), 0.);  // initialize to 0.
}
void TimerWriter::endStep(ParticleContainer* /*particleContainer*/, DomainDecompBase* /*domainDecomp*/,
						  Domain* /*domain*/, unsigned long simstep) {
	for (size_t i = 0; i < _timerNames.size(); ++i) {
		if (not _incremental[i]) {
			// save times for each time step, so we can average it later.
			_times[i] += global_simulation->timers()->getTimer(_timerNames[i])->get_etime();
		}
		// nothing to do for incremental timers, as they track their overall time themselves.
	}
	// increase amount of written steps
	++_stepsSinceLastWrite;

	if (simstep % _writeFrequency == 0 and simstep != 0) {
		_fileStream << simstep;
		for (size_t i = 0; i < _timerNames.size(); ++i) {
			double timeOneStep;
			if (_incremental[i]) {
				double currentTime = global_simulation->timers()->getTimer(_timerNames[i])->get_etime();
				timeOneStep = (currentTime - _times[i]) / static_cast<double>(_stepsSinceLastWrite);
				// incremental timer, so set to current time!
				_times[i] = currentTime;
			} else {
				// we have added _stepsSinceLastWrite individual timing values, so divide by this number!
				timeOneStep = _times[i] / static_cast<double>(_stepsSinceLastWrite);
				// reset time to 0, so we can add them again.
				_times[i] = 0;
			}
			_fileStream << ";" << timeOneStep;
		}
		_fileStream << std::endl;

		// reset write state
		_stepsSinceLastWrite = 0ul;
	}
}
