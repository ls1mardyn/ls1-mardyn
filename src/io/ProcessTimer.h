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

    ProcessTimer(){}

    ~ProcessTimer(){}

    //! @brief Inserts new time for given Process-Rank
    //!
    //! @param process Rank of process
    //! @param time Time of given process
    void insertTime(int process, double time){
        if( _processes.count(process) > 0)
            _processes[process].push_back(time);
        else{
            std::vector<double> tempVec ={time};
            _processes.insert(std::pair<int, std::vector<double>>(process, tempVec));
        }
    }

    //! @brief Prints out whole map; Contains Ranks incl measured Times
    void printTimeOnConsole(){
        for (auto const& iter : _processes)
        {
            for (auto innerIter = iter.second.begin() ; innerIter != iter.second.end(); ++innerIter){
                std::cout << "Rank: " << iter.first << " Runtime: " << *innerIter <<" seconds" << std::endl;
            }
        }
    }

    //! @brief Writes whole map into "processRuntime.txt" in current directory
    void writeProcessTimeLog(){
        std::ofstream _processRuntime;
        _processRuntime.open("processRuntime.txt");
        for (auto const& iter : _processes)
        {
            for (auto innerIter = iter.second.begin() ; innerIter != iter.second.end(); ++innerIter){
                _processRuntime  << "Rank: " << iter.first << " Runtime: " << *innerIter <<" seconds" << std::endl;
            }
        }
        _processRuntime.close();
    }
protected:
private:
    std::map<int, std::vector<double>> _processes;



};
#endif //MARDYN_PROCESSTIMER_H
