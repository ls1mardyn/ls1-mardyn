/*
 * ProcessTimer.h
 *
 *  Created on: Apr 12, 2020
 *      Author: jeremyharisch
 */

#ifndef MARDYN_PROCESSTIMER_H
#define MARDYN_PROCESSTIMER_H


#include <fstream>
#include <vector>
#include <array>
#include <map>
#include <time.h>

#include <mpi.h>

#include "utils/Timer.h"

class ProcessTimer {

public:
	ProcessTimer() {}

	~ProcessTimer() {}


	//! @brief Starts Timer for MPI-Measurement
	//!
	//! @param process Rank of process
	void startTimer() {
		int process;
		MPI_Comm_rank(MPI_COMM_WORLD, &process);
		_process_time[process] -= _process_time[process];
#ifdef ENABLE_MPI
		double measurement_time = MPI_Wtime(); // TODO: Entfernen; da alte Struktur_process_time[process] -= measurement_time;
#else
		struct timeval tmp_time;
		gettimeofday(&tmp_time, NULL);
		double measurement_time = (1.0e6 * (double) tmp_time.tv_sec + (double) tmp_time.tv_usec) / 1.0e6; // TODO: Entfernen; da alte Struktur
#endif
		_process_time[process] -= measurement_time;
		_processes_debug[process].push_back(-measurement_time);
	}

	//! @brief Stops Timer for MPI-Measurement
	void stopTimer() {
		int process;
		MPI_Comm_rank(MPI_COMM_WORLD, &process);
#ifdef ENABLE_MPI
		double measurement_time = MPI_Wtime();
#else
		struct timeval tmp_time;
		gettimeofday(&tmp_time, NULL);
		double measurement_time = (1.0e6 * (double) tmp_time.tv_sec + (double) tmp_time.tv_usec) / 1.0e6;
#endif
		_process_time[process] += measurement_time;
		_processes_debug[process][_processes_debug[process].size()-1] -= measurement_time;
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
			_processRuntime.open("process" + std::to_string(process) + "Runtime.csv", std::ios::app);
			_processRuntime << time << ", ";
		}
		else{
			_processRuntime.open("process" + std::to_string(process) + "Runtime.txt", std::ios::app);
			_processRuntime << "Rank: " << process << " Runtime: " << time << " seconds\n";
		}
        _processRuntime.close();
	}

	//! @brief Writes whole debug-map into "processRuntime.txt" in current directory
	void writeProcessTimeLog() {
		std::ofstream _processRuntime;
		_processRuntime.open("processDebugRuntime.txt");
		for (auto const &iter : _processes_debug) {
			for (auto innerIter = iter.second.begin(); innerIter != iter.second.end(); ++innerIter) {
				_processRuntime << "Rank: " << iter.first << " Runtime: " << *innerIter << std::endl;
			}
		}
		_processRuntime.close();
	}

	double lastProcessTime = 0.0;

protected:
private:
	std::map<int, double> _process_time;
	std::map<int, std::vector<double>> _processes_debug;
};

#endif //MARDYN_PROCESSTIMER_H
