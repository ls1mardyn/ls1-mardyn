/*
 * TimerProfiler.cpp
 *
 *  Created on: Apr 9, 2017
 *      Author: Andrei Costinescu
 */

#include <cmath>
#include <tuple>
#include <sstream>

#include "TimerProfiler.h"
#include "utils/Logger.h"
#include "utils/String_utils.h"
#include "utils/xmlfileUnits.h"
#include "utils/mardyn_assert.h"



const std::string TimerProfiler::_baseTimerName = "_baseTimer";

TimerProfiler::TimerProfiler(): _numElapsedIterations(0), _displayMode(Displaymode::ALL) {
	_timers[_baseTimerName] = _Timer(_baseTimerName);
	readInitialTimersFromFile("");
}

void TimerProfiler::readXML(XMLfileUnits& xmlconfig) {
	std::string displayMode;
	if(xmlconfig.getNodeValue("displaymode", displayMode)) {
		Log::global_log->info() << "Timer display mode: " << displayMode << std::endl;
		if(displayMode == "all") {
			setDisplayMode(Displaymode::ALL);
		} else if (displayMode == "active") {
			setDisplayMode(Displaymode::ACTIVE);
		} else if (displayMode == "non-zero") {
			setDisplayMode(Displaymode::NON_ZERO);
		} else if (displayMode == "none") {
			setDisplayMode(Displaymode::NONE);
		} else {
			std::ostringstream error_message;
			error_message << "Unknown display mode: " << displayMode << std::endl;
			MARDYN_EXIT(error_message.str());
		}
	}
}


Timer* TimerProfiler::getTimer(std::string timerName){
	auto timerProfiler = _timers.find(timerName);
	if(timerProfiler != _timers.end()) {
		return (timerProfiler->second)._timer.get();
	}
	return nullptr;
}

void TimerProfiler::registerTimer(std::string timerName, std::vector<std::string> parentTimerNames, Timer *timer, bool activate){
	Log::global_log->debug() << "Registering timer: " << timerName << "  [parents: " << string_utils::join(parentTimerNames, std::string(", ")) << "]" << std::endl;

	if (!activate && timer){
		timer->deactivateTimer();
	}
	_timers[timerName] = _Timer(timerName, timer);

	if (parentTimerNames.empty()) {
		parentTimerNames.push_back(_baseTimerName);
	}
	for (const auto& parentTimerName : parentTimerNames) {
		_timers[timerName]._parentTimerNames.push_back(parentTimerName);
		_timers[parentTimerName]._childTimerNames.push_back(timerName);
	}
}

void TimerProfiler::activateTimer(std::string timerName){
	if (!_timers.count(timerName)) return ;
	_timers[timerName]._timer->activateTimer();
}

void TimerProfiler::deactivateTimer(std::string timerName){
	if (!_timers.count(timerName)) return ;
	_timers[timerName]._timer->deactivateTimer();
}

void TimerProfiler::setSyncTimer(std::string timerName, bool sync){
	if (!_timers.count(timerName)) return ;
	_timers[timerName]._timer->set_sync(sync);
}

void TimerProfiler::print(std::string timerName, std::string outputPrefix){
	if ( ! _checkTimer(timerName, false)) {
		_debugMessage(timerName);
		return;
	}
	if( (getDisplayMode() == Displaymode::ALL) ||
		(getDisplayMode() == Displaymode::ACTIVE && getTimer(timerName)->isActive()) ||
		(getDisplayMode() == Displaymode::NON_ZERO && getTimer(timerName)->get_etime() > 0)
	) {
		Log::global_log->info() << outputPrefix << getOutputString(timerName) << getTime(timerName) << " sec" << std::endl;
	}
}

void TimerProfiler::printTimers(std::string timerName, std::string outputPrefix){
	if (!_timers.count(timerName)) return ;
	if (_checkTimer(timerName)){
		print(timerName, outputPrefix);
		outputPrefix += "\t";
	}
	for(const auto& childTimerName : _timers[timerName]._childTimerNames){
		printTimers(childTimerName, outputPrefix);
	}
}

void TimerProfiler::start(std::string timerName){
	#ifdef _OPENMP
	#pragma omp critical
	#endif
	{
		if (_checkTimer(timerName)) {
			getTimer(timerName)->start();
		} else {
			_debugMessage(timerName);
		}
	}
}

void TimerProfiler::stop(std::string timerName){
	#ifdef _OPENMP
	#pragma omp critical
	#endif
	{
		if (_checkTimer(timerName)) {
			getTimer(timerName)->stop();
		} else {
			_debugMessage(timerName);
		}
	}
}

void TimerProfiler::reset(std::string timerName){
	if (_checkTimer(timerName)){
		getTimer(timerName)->reset();
	}
	else{
		_debugMessage(timerName);
	}
}

void TimerProfiler::resetTimers(std::string startingTimerName){
	if (startingTimerName.compare(_baseTimerName)) {
		_numElapsedIterations = 0;
	}

	if (!_timers.count(startingTimerName)) return ;
	reset(startingTimerName);
	for(auto childTimerName : _timers[startingTimerName]._childTimerNames){
		resetTimers(childTimerName);
	}
}

void TimerProfiler::readInitialTimersFromFile(std::string fileName){
	//temporary until read from .xml file is implemented

	//timer classes for grouping timers -> easier printing and resetting
	std::vector<std::string> timerClasses = {
		"COMMUNICATION_PARTNER",
		"SIMULATION",
		"UNIFORM_PSEUDO_PARTICLE_CONTAINER",
		"GENERATORS",
		"CELL_PROCESSORS",
		"TUNERS",
		"IO"
	};

	/**************************************************************
	* There are 4 unused timers in UniformPseudoParticleContainer
	* 1) UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_CELLS -> there were calls to set_sync and to reset, but no start/stop
	* 2) UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_FAR_FIELD -> there were calls to set_sync and to reset, but no start/stop
	* 3) UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_M -> not used at all...
	* 4) UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_LM -> not used at all...
	**************************************************************/
	std::vector<std::tuple<std::string, std::vector<std::string>, bool>> timerAttrs = {
		std::make_tuple("VECTORIZATION_TUNER_TUNER", std::vector<std::string>{"TUNERS"}, true),
		std::make_tuple("AQUEOUS_NA_CL_GENERATOR_INPUT", std::vector<std::string>{"GENERATORS"}, true),
		std::make_tuple("CRYSTAL_LATTICE_GENERATOR_INPUT", std::vector<std::string>{"GENERATORS"}, true),
		std::make_tuple("CUBIC_GRID_GENERATOR_INPUT", std::vector<std::string>{"GENERATORS"}, true),
		std::make_tuple("DROPLET_GENERATOR_INPUT", std::vector<std::string>{"GENERATORS"}, true),
		std::make_tuple("MS2RST_GENERATOR_INPUT", std::vector<std::string>{"GENERATORS"}, true),
		std::make_tuple("REYLEIGH_TAYLOR_GENERATOR_INPUT", std::vector<std::string>{"GENERATORS"}, true),
		std::make_tuple("REPLICA_GENERATOR_VLE_INPUT", std::vector<std::string>{"GENERATORS"}, true),
		std::make_tuple("L2P_CELL_PROCESSOR_L2P", std::vector<std::string>{"CELL_PROCESSORS"}, true),
		std::make_tuple("P2M_CELL_PROCESSOR_P2M", std::vector<std::string>{"CELL_PROCESSORS"}, true),
		std::make_tuple("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P", std::vector<std::string>{"CELL_PROCESSORS"}, true),
		std::make_tuple("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P", std::vector<std::string>{"CELL_PROCESSORS"}, true),
		std::make_tuple("BINARY_READER_INPUT", std::vector<std::string>{"IO"}, true),
		std::make_tuple("INPUT_OLDSTYLE_INPUT", std::vector<std::string>{"IO"}, true),
		std::make_tuple("MPI_IO_READER_INPUT", std::vector<std::string>{"IO"}, true),
		std::make_tuple("MPI_CHECKPOINT_WRITER_INPUT", std::vector<std::string>{"IO"}, true),
		std::make_tuple("SIMULATION_LOOP", std::vector<std::string>{"SIMULATION"}, true),
		std::make_tuple("SIMULATION_DECOMPOSITION", std::vector<std::string>{"SIMULATION_LOOP"}, true),
		std::make_tuple("SIMULATION_COMPUTATION", std::vector<std::string>{"SIMULATION_LOOP"}, true),
#ifdef QUICKSCHED
		std::make_tuple("QUICKSCHED", std::vector<std::string>{"SIMULATION_LOOP"}, true),
#endif
		std::make_tuple("SIMULATION_PER_STEP_IO", std::vector<std::string>{"SIMULATION_LOOP"}, true),
		std::make_tuple("SIMULATION_IO", std::vector<std::string>{"SIMULATION"}, true),
		std::make_tuple("SIMULATION_UPDATE_CONTAINER", std::vector<std::string>{"SIMULATION_DECOMPOSITION"}, true),
		std::make_tuple("SIMULATION_MPI_OMP_COMMUNICATION", std::vector<std::string>{"SIMULATION_DECOMPOSITION"}, true),
		std::make_tuple("SIMULATION_UPDATE_CACHES", std::vector<std::string>{"SIMULATION_DECOMPOSITION"}, true),
		std::make_tuple("SIMULATION_FORCE_CALCULATION", std::vector<std::string>{"SIMULATION_COMPUTATION"}, true),
		std::make_tuple("COMMUNICATION_PARTNER_INIT_SEND", std::vector<std::string>{"COMMUNICATION_PARTNER", "SIMULATION_MPI_OMP_COMMUNICATION"}, true),
		std::make_tuple("COMMUNICATION_PARTNER_TEST_RECV", std::vector<std::string>{"COMMUNICATION_PARTNER", "SIMULATION_MPI_OMP_COMMUNICATION"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_CELLS", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_GLOBAL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_LOKAL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_LOKAL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_FAR_FIELD", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMMUNICATION_HALOS", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_HALO_GATHER", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_CALCULATION", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_INIT", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_FINALIZE", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_M", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_LM", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE_ME", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		std::make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL", std::vector<std::string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true)
	};

	for (auto timerClass : timerClasses){
		registerTimer(timerClass, std::vector<std::string>());
	}
	for (auto timerAttr : timerAttrs){
		registerTimer(std::get<0>(timerAttr), std::get<1>(timerAttr), new Timer(), std::get<2>(timerAttr));
	}
}

void TimerProfiler::setOutputString(std::string timerName, std::string outputString){
	if (!_timers.count(timerName)) return ;
	if (outputString[outputString.length() - 1] != ' ') {
		outputString += " ";
	}
	_timers[timerName]._outputString = outputString;
}

std::string TimerProfiler::getOutputString(std::string timerName){
	if (!_timers.count(timerName)) return std::string("");
	std::string output = _timers[timerName]._outputString;
	if (output.compare("") == 0){
		output = "Timer "+timerName+" took: ";
	}
	return output;
}

double TimerProfiler::getTime(std::string timerName){
	auto timer = getTimer(timerName);
	if(timer != nullptr) {
		return timer->get_etime();
	}
	return 0.0;
}

/* private Methods */

bool TimerProfiler::_checkTimer(std::string timerName, bool checkActive){
	return _timers.count(timerName) && _timers[timerName]._timer && (!checkActive || _timers[timerName]._timer->isActive());
}

void TimerProfiler::_debugMessage(std::string timerName){
	if(_timers.count(timerName)){
		Log::global_log->debug()<<"Timer "<<timerName<<" is not "<<(!_checkTimer(timerName, false) ? "a timer" : "active")<<".\n";
	}
	else{
		Log::global_log->debug()<<"Timer "<<timerName<<" is not registered."<< std::endl;
	}
}
