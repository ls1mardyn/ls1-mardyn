/*
 * ProcessTimer.h
 *
 *  Created on: Apr 12, 2020
 *      Author: jeremyharisch
 */

#ifndef MARDYN_PROCESSTIMER_H
#define MARDYN_PROCESSTIMER_H
#define _BSD_SOURCE

#include <fstream>
#include <vector>
#include <array>
#include <map>
#include <sys/time.h>

#include <mpi.h>


class ProcessTimer {

public:
	ProcessTimer() {}

	~ProcessTimer() {}

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
		if (_profiling_switch == 1){
			int process;
			MPI_Comm_rank(MPI_COMM_WORLD, &process);
			double measurement_time = MPI_Wtime();
			_process_time[process] -= measurement_time;
			_processes_debug[process].push_back(-measurement_time);
		}

	}

	//! @brief Stops Timer for MPI-Measurement
	void stopTimer() {
		if (_profiling_switch ==1){
			int process;
			MPI_Comm_rank(MPI_COMM_WORLD, &process);
			double measurement_time = MPI_Wtime();
			_process_time[process] += measurement_time;
			_processes_debug[process][_processes_debug[process].size()-1] -= measurement_time;
		}
	}

	//! @brief Returns the Time saved with this Timer for given Process
	//! @param process Process for return time
	//! @param debug Prints time into csv
	//! @param reset Determines if Timer should be reseted
	int getTime(int process, bool reset = false, bool debug = false){
		int time = _process_time[process];
		if (debug)
			writeProcessTimeLogSingle(process, _process_time[process], true);
		if (reset)
			_process_time[process] -= _process_time[process];
		return time;
	}

	//! @brief Writes given Process and time spent in MPI-Calls into "processRuntime.txt"
	//! @param csv Time are written ether to a csv or txt file
	void writeProcessTimeLogSingle(int process, double time, bool csv) {
		std::ofstream _processRuntime;
		if(csv){
			_processRuntime.open("process" + std::to_string(process) + "MPI-Runtime.csv", std::ios::app);
			_processRuntime << time << ", ";
		}
		else{
			_processRuntime.open("process" + std::to_string(process) + "MPI-Runtime.txt", std::ios::app);
			_processRuntime << "Rank: " << process << " Runtime: " << time << " seconds\n";
		}
        _processRuntime.close();
	}

	//! @brief Writes whole debug-map into "processRuntime.txt" in current directory
	void writeProcessTimeLog() {
		std::ofstream _processRuntime;
		_processRuntime.open("processDebugMPI-Runtime.txt");
		for (auto const &iter : _processes_debug) {
			for (auto innerIter = iter.second.begin(); innerIter != iter.second.end(); ++innerIter) {
				_processRuntime << "Rank: " << iter.first << " Runtime: " << *innerIter << std::endl;
			}
		}
		_processRuntime.close();
	}

	//! @brief Resets the time of given process
	//! @param process ID of process
	void resetTimer(int process) {
        _process_time[process] -= _process_time[process];
	}

protected:
private:
	std::map<int, double> _process_time;
	std::map<int, std::vector<double>> _processes_debug;
	int _profiling_switch = 1;
};

#endif //MARDYN_PROCESSTIMER_H
