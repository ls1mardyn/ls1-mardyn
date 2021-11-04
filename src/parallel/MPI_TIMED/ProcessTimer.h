/*
 * ProcessTimer.h
 *
 *  Created on: Apr 12, 2020
 *      Author: jeremyharisch
 */

#ifndef MARDYN_PROCESSTIMER_H
#define MARDYN_PROCESSTIMER_H

#include <array>
#include <fstream>
#include <map>
#include <vector>

#include <mpi.h>

class ProcessTimer {
public:
	static ProcessTimer& get() {
		static ProcessTimer processTimer;
		return processTimer;
	}

	/*! @brief Switches between the profiling Levels
	 *  This Fnction should be called right before and after a corresponding MPI-Function-Call
	 *  to exclude only this one
	 *
	 *  @param level ==1 Profiling is enabled(default); ==0 Profiling is disabled
	 */
	int switchProfiling(const int level = 1) {
		_profiling_switch = level;
		return _profiling_switch;
	}

	//! @brief Starts Timer for MPI-Measurement
	//!
	//! @param process Rank of process
	void startTimer() {
		if (_profiling_switch == 1) {
			double measurement_time = MPI_Wtime();
			_process_time -= measurement_time;
			_processes_debug.push_back(-measurement_time);
		}
	}

	//! @brief Stops Timer for MPI-Measurement
	void stopTimer() {
		if (_profiling_switch == 1) {
			double measurement_time = MPI_Wtime();
			_process_time += measurement_time;
			_processes_debug[_processes_debug.size() - 1] -= measurement_time;
		}
	}

	//! @brief Returns the Time saved with this Timer for given Process
	//! @param process Process for return time
	//! @param debug Prints time into csv
	//! @param reset Determines if Timer should be reseted
	double getTime(bool reset = false, bool debug = false) {
		double time = _process_time;
		if (debug) {
			int process;
			MPI_Comm_rank(MPI_COMM_WORLD, &process);
			writeProcessTimeLogSingle(process, _process_time, true);
		}
		if (reset) {
			resetTimer();
		}
		return time;
	}

	//! @brief Writes given Process and time spent in MPI-Calls into "processRuntime.txt"
	//! @param csv Time are written ether to a csv or txt file
	void writeProcessTimeLogSingle(int process, double time, bool csv) {
		if (!_processRuntimeSingle.is_open()) {
			if (csv) {
				_processRuntimeSingle.open("process" + std::to_string(process) + "MPI-Runtime.csv", std::ios::app);
			} else {
				_processRuntimeSingle.open("process" + std::to_string(process) + "MPI-Runtime.txt", std::ios::app);
			}
		}
		if (csv) {
			_processRuntimeSingle << time << ", ";
		} else {
			_processRuntimeSingle << "Rank: " << process << " Runtime: " << time << " seconds\n";
		}
	}

	//! @brief Writes whole debug-map into "processRuntime.txt" in current directory
	void writeProcessTimeLog() {
		if (!_processRuntime.is_open()) {
			_processRuntime.open("processDebugMPI-Runtime.txt");
		}
		for (auto const& iter : _processes_debug) {
			_processRuntime << " Runtime: " << iter << std::endl;
		}
	}

	//! @brief Resets the time of given process
	void resetTimer() { _process_time -= _process_time; }

private:
	double _process_time {0.};
	std::vector<double> _processes_debug;
	int _profiling_switch{1};
	std::ofstream _processRuntime;
	std::ofstream _processRuntimeSingle;
};

#endif  // MARDYN_PROCESSTIMER_H
