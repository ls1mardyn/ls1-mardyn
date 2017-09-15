/*
 * TimerProfiler.cpp
 *
 *  Created on: Apr 9, 2017
 *      Author: Andrei Costinescu
 */

#include <assert.h>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <thread>
#include <chrono>
#include <cmath>

#include "TimerProfiler.h"
#include "utils/Logger.h"
#include "utils/String_utils.h"
#include "mardyn_assert.h"

using namespace std;
using Log::global_log;

const string TimerProfiler::_baseTimerName = "_baseTimer";

TimerProfiler::TimerProfiler(): _numElapsedIterations(0) {
	_timers[_baseTimerName] = _Timer(_baseTimerName);
	readInitialTimersFromFile("");
}

TimerProfiler::~TimerProfiler() {
	_clearTimers();
	mardyn_assert(_timers.size() == 0);
}

Timer* TimerProfiler::getTimer(string timerName){
	if (!_timers.count(timerName)) throw invalid_argument("No timer with name "+timerName+".");
	if (!_checkTimer(timerName, false)) throw invalid_argument("Timer "+timerName+" is not a timer");
	return _timers[timerName]._timer;
}

void TimerProfiler::registerTimer(string timerName, vector<string> parentTimerNames, Timer *timer, bool activate){
	global_log->debug() << "Registering timer: " << timerName << "  [parents: " << string_utils::join(parentTimerNames, string(", ")) << "]" << endl;

	if (!activate && timer){
		timer->deactivateTimer();
	}
	_timers[timerName] = _Timer(timerName, timer);

	if (!parentTimerNames.size()) parentTimerNames.push_back(_baseTimerName);
	for (size_t i=0; i<parentTimerNames.size(); i++) {
		_timers[parentTimerNames[i]]._childTimerNames.push_back(timerName);
		_timers[timerName]._parentTimerNames.push_back(parentTimerNames[i]);
	}
}

void TimerProfiler::activateTimer(string timerName){
	if (!_timers.count(timerName)) return ;
	_timers[timerName]._timer->activateTimer();
}

void TimerProfiler::deactivateTimer(string timerName){
	if (!_timers.count(timerName)) return ;
	_timers[timerName]._timer->deactivateTimer();
}

void TimerProfiler::setSyncTimer(string timerName, bool sync){
	if (!_timers.count(timerName)) return ;
	_timers[timerName]._timer->set_sync(sync);
}

void TimerProfiler::print(string timerName, string outputPrefix){
	if (_checkTimer(timerName)) {
			global_log->info()<<outputPrefix<<getOutputString(timerName)<<_timers[timerName]._timer->get_etime()<<" sec\n";
	}
	else {
		_debugMessage(timerName);
	}
}

void TimerProfiler::printTimers(string startingTimerName, string outputPrefix){
	if (!_timers.count(startingTimerName)) return ;
	if (_checkTimer(startingTimerName)){
		print(startingTimerName, outputPrefix);
		outputPrefix += "\t";
	}
	for(size_t i=0; i<_timers[startingTimerName]._childTimerNames.size(); i++){
		printTimers(_timers[startingTimerName]._childTimerNames[i], outputPrefix);
	}
}

void TimerProfiler::start(string timerName){
	if (_checkTimer(timerName)){
		_timers[timerName]._timer->start();
	}
	else{
		_debugMessage(timerName);
	}
}

void TimerProfiler::stop(string timerName){
	if (_checkTimer(timerName)){
		_timers[timerName]._timer->stop();
	}
	else{
		_debugMessage(timerName);
	}
}

void TimerProfiler::reset(string timerName){
	if (_checkTimer(timerName)){
		_timers[timerName]._timer->reset();
	}
	else{
		_debugMessage(timerName);
	}
}

void TimerProfiler::resetTimers(string startingTimerName){
	if (startingTimerName.compare("_baseTimer")) {
		_numElapsedIterations = 0;
	}

	if (!_timers.count(startingTimerName)) return ;
	reset(startingTimerName);
	for(size_t i=0; i<_timers[startingTimerName]._childTimerNames.size(); i++){
		resetTimers(_timers[startingTimerName]._childTimerNames[i]);
	}
}

void TimerProfiler::readInitialTimersFromFile(string fileName){
	//temporary until read from .xml file is implemented

	//timer classes for grouping timers -> easier printing and resetting
	vector<string> timerClasses = {
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
	vector<tuple<string, vector<string>, bool>> timerAttrs = {
		make_tuple("VECTORIZATION_TUNER_TUNER", vector<string>{"TUNERS"}, true),
		make_tuple("AQUEOUS_NA_CL_GENERATOR_INPUT", vector<string>{"GENERATORS"}, true),
		make_tuple("CRYSTAL_LATTICE_GENERATOR_INPUT", vector<string>{"GENERATORS"}, true),
		make_tuple("CUBIC_GRID_GENERATOR_INPUT", vector<string>{"GENERATORS"}, true),
		make_tuple("DROPLET_GENERATOR_INPUT", vector<string>{"GENERATORS"}, true),
		make_tuple("MS2RST_GENERATOR_INPUT", vector<string>{"GENERATORS"}, true),
		make_tuple("REYLEIGH_TAYLOR_GENERATOR_INPUT", vector<string>{"GENERATORS"}, true),
		make_tuple("REPLICA_GENERATOR_VLE_INPUT", vector<string>{"GENERATORS"}, true),
		make_tuple("L2P_CELL_PROCESSOR_L2P", vector<string>{"CELL_PROCESSORS"}, true),
		make_tuple("P2M_CELL_PROCESSOR_P2M", vector<string>{"CELL_PROCESSORS"}, true),
		make_tuple("VECTORIZED_CHARGE_P2P_CELL_PROCESSOR_VCP2P", vector<string>{"CELL_PROCESSORS"}, true),
		make_tuple("VECTORIZED_LJP2P_CELL_PROCESSOR_VLJP2P", vector<string>{"CELL_PROCESSORS"}, true),
		make_tuple("BINARY_READER_INPUT", vector<string>{"IO"}, true),
		make_tuple("INPUT_OLDSTYLE_INPUT", vector<string>{"IO"}, true),
		make_tuple("MPI_IO_READER_INPUT", vector<string>{"IO"}, true),
		make_tuple("MPI_CHECKPOINT_WRITER_INPUT", vector<string>{"IO"}, true),
		make_tuple("SIMULATION_LOOP", vector<string>{"SIMULATION"}, true),
		make_tuple("SIMULATION_DECOMPOSITION", vector<string>{"SIMULATION_LOOP"}, true),
		make_tuple("SIMULATION_COMPUTATION", vector<string>{"SIMULATION_LOOP"}, true),
		make_tuple("SIMULATION_PER_STEP_IO", vector<string>{"SIMULATION_LOOP"}, true),
		make_tuple("SIMULATION_IO", vector<string>{"SIMULATION"}, true),
		make_tuple("SIMULATION_MPI_OMP_COMMUNICATION", vector<string>{"SIMULATION_DECOMPOSITION"}, true),
		make_tuple("SIMULATION_FORCE_CALCULATION", vector<string>{"SIMULATION_COMPUTATION"}, true),
		make_tuple("COMMUNICATION_PARTNER_INIT_SEND", vector<string>{"COMMUNICATION_PARTNER", "SIMULATION_MPI_OMP_COMMUNICATION"}, true),
		make_tuple("COMMUNICATION_PARTNER_TEST_RECV", vector<string>{"COMMUNICATION_PARTNER", "SIMULATION_MPI_OMP_COMMUNICATION"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_CELLS", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_GLOBAL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_LOKAL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_LOKAL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_FAR_FIELD", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMMUNICATION_HALOS", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_HALO_GATHER", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_CALCULATION", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_INIT", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_FINALIZE", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_M", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_LM", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE_ME", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true),
		make_tuple("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL", vector<string>{"UNIFORM_PSEUDO_PARTICLE_CONTAINER"}, true)
	};

	for (auto timerClass : timerClasses){
		registerTimer(timerClass, vector<string>());
	}
	for (auto timerAttr : timerAttrs){
		registerTimer(get<0>(timerAttr), get<1>(timerAttr), new Timer(), get<2>(timerAttr));
	}
}

void TimerProfiler::setOutputString(string timerName, string outputString){
	if (!_timers.count(timerName)) return ;
	if (outputString[outputString.length() - 1] != ' ') {
		outputString += " ";
	}
	_timers[timerName]._outputString = outputString;
}

string TimerProfiler::getOutputString(string timerName){
	if (!_timers.count(timerName)) throw invalid_argument("No timer with name "+timerName+".");
	string output = _timers[timerName]._outputString;
	if (output.compare("") == 0){
		output = "Timer "+timerName+" took: ";
	}
	return output;
}

double TimerProfiler::getTime(string timerName){
	if (!_timers.count(timerName)) throw invalid_argument("No timer with name "+timerName+".");
	if (!_checkTimer(timerName, false)) throw invalid_argument("Timer "+timerName+" is not a timer");
	if (!_checkTimer(timerName)) throw invalid_argument("Timer "+timerName+" is not active");
	return _timers[timerName]._timer->get_etime();
}

/* private Methods */

bool TimerProfiler::_checkTimer(string timerName, bool checkActive){
	return _timers.count(timerName) && _timers[timerName]._timer && (!checkActive || _timers[timerName]._timer->isActive());
}

void TimerProfiler::_debugMessage(string timerName){
	if(_timers.count(timerName)){
		global_log->debug()<<"Timer "<<timerName<<" is not "<<(!_checkTimer(timerName, false) ? "a timer" : "active")<<".\n";
	}
	else{
		global_log->debug()<<"Timer "<<timerName<<" is not registered."<< endl;
	}
}

void TimerProfiler::_clearTimers(string startingTimerName){
	//if the timer is not in the container it must have already been deleted -> return
	if (!_timers.count(startingTimerName)){
		return;
	}
	for(size_t i=0; i<_timers[startingTimerName]._childTimerNames.size(); i++){
		_clearTimers(_timers[startingTimerName]._childTimerNames[i]);
	}
	if (_checkTimer(startingTimerName, false)) {
		delete _timers[startingTimerName]._timer;
	}
	_timers[startingTimerName]._childTimerNames.clear();
	_timers[startingTimerName]._parentTimerNames.clear();
	_timers.erase(startingTimerName);
}
