/*
 * ProcessTimer.h
 *
 *  Created on: Apr 12, 2020
 *      Author: andal
 */

#ifndef MARDYN_PROCESSTIMER_H
#define MARDYN_PROCESSTIMER_H

#include <fstream>
#include <vector>
#include <array>
#include <map>

#include "utils/Timer.h"

class ProcessTimer {

public:
	ProcessTimer() {}

	~ProcessTimer() {}

	//! @brief Inserts new time for given Process-Rank
	//!
	//! @param process Rank of process
	//! @param time Time of given process
	void insertTime(int process, double time) {
		_processes[process].push_back(time);
		writeProcessTimeLogSingle(process, time, true);
	}

	//! @brief Prints out whole map; Contains Ranks incl measured Times
	void printTimeOnConsole() {
		for (auto const &iter : _processes) {
			for (auto innerIter = iter.second.begin(); innerIter != iter.second.end(); ++innerIter) {
				std::cout << "Rank: " << iter.first << " Runtime: " << *innerIter << " seconds" << std::endl;
			}
		}
	}

	//! @brief Writes whole map into "processRuntime.txt" in current directory
	void writeProcessTimeLog() {
		std::ofstream _processRuntime;
		_processRuntime.open("processRuntime.txt");
		for (auto const &iter : _processes) {
			for (auto innerIter = iter.second.begin(); innerIter != iter.second.end(); ++innerIter) {
				_processRuntime << "Rank: " << iter.first << " Runtime: " << *innerIter << " seconds" << std::endl;
			}
		}
		_processRuntime.close();
	}

	//! @brief Writes given Process and time into "processRuntime.txt"
	//! @param csv Time are written ether to a csv or txt file
	void writeProcessTimeLogSingle(int process, double time, bool csv) {
		std::ofstream _processRuntime;
		if(csv){
			_processRuntime.open("process" + std::to_string(process) + "Runtime.csv", std::ios::app);
			_processRuntime << time << ", ";
		}
		else{
			// For each Process one File
			_processRuntime.open("process" + std::to_string(process) + "Runtime.txt", std::ios::app);
			// All Runtimes in one File, this solution works with MPI as well
			//_processRuntime.open("processRuntime.txt", std::ios::app);
			_processRuntime << "Rank: " << process << " Runtime: " << time << " seconds\n";
		}
        _processRuntime.close();
	}

protected:
private:
	std::map<int, std::vector<double>> _processes;


};

#endif //MARDYN_PROCESSTIMER_H
