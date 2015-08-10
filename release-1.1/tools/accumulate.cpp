
#include "utils/Accumulator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

/**
 * This tool computes accumulated and averaged values for pressure and potential 
 * energy, from Mardyn result files.
 */
int main(int argc, char** argsv) {

  if (argc != 3) {
    cout << "Usage: ./accumulate FILENAME NTIMES" << endl;
    return -1;
  }

  int numSteps = atoi(argsv[2]);
  cout << " Using file " << argsv[1] <<", averaging last " << numSteps << " steps." << endl;

  ifstream inputfilestream(argsv[1]);
  if (!inputfilestream.is_open()) {
    cerr << "Could not open file " << argsv[1] << endl;
    exit(1);
  }
  
  Accumulator<double> upot_acc(numSteps);
  Accumulator<double> p_acc(numSteps);

  string line;
  while(inputfilestream) {
    line.clear();
    getline(inputfilestream, line);
  
     if (line.empty() || line[0] == '#') {
       cout << " skippiing: " << line << endl;
       continue;
    }

    stringstream lineStream(line);
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

  cout << "==================  UPot = " << upot_acc.getAverage() << ", p = " << p_acc.getAverage() << endl;
  return 0;
}

