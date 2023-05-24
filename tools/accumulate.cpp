
#include "utils/Accumulator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


/**
 * This tool computes accumulated and averaged values for pressure and potential 
 * energy, from Mardyn result files.
 */
int main(int argc, char** argsv) {

  if (argc != 3) {
    std::cout << "Usage: ./accumulate FILENAME NTIMES" << std::endl;
    return -1;
  }

  int numSteps = atoi(argsv[2]);
  std::cout << " Using file " << argsv[1] <<", averaging last " << numSteps << " steps." << std::endl;

  std::ifstream inputfilestream(argsv[1]);
  if (!inputfilestream.is_open()) {
    std::cerr << "Could not open file " << argsv[1] << std::endl;
    exit(1);
  }
  
  Accumulator<double> upot_acc(numSteps);
  Accumulator<double> p_acc(numSteps);

  std::string line;
  while(inputfilestream) {
    line.clear();
    getline(inputfilestream, line);
  
     if (line.empty() || line[0] == '#') {
       std::cout << " skippiing: " << line << std::endl;
       continue;
    }

    std::stringstream lineStream(line);
    double upot = 0; 
    double p = 0;
    double dummy = 0;

    lineStream >> dummy >> dummy; // ignore '#step' and 't'
    lineStream >> upot;
    lineStream >> dummy; // ignore 'p_avg'
    lineStream >> p;
    upot_acc.addEntry(upot);
    p_acc.addEntry(p);
  } 

  std::cout << "==================  UPot = " << upot_acc.getAverage() << ", p = " << p_acc.getAverage() << std::endl;
  return 0;
}

