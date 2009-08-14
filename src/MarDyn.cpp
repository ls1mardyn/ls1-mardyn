#include <iostream>
#include <iomanip>
#include <ctime>
#include "Simulation.h"

using namespace std;

//! @page main
//! In this project, software for molecular dynamics simulation
//! with short-range forces is developed. The aim is to have a parallel code (MPI) 
//! for multi-centered molecules.
//!
//! The role of the main function is to run tests for all classes
//! and to instantiate an object of the Simulation class which
//! is actually responsible for the simulation
//!
int main(int argc, char** argv){

  cout.precision(6);

  Simulation simulation(&argc, &argv);
  
  double runtime = double(clock())/CLOCKS_PER_SEC;

  simulation.simulate();

  runtime=double(clock())/CLOCKS_PER_SEC-runtime;
  
  cout << "main: used " << fixed << setprecision(2) << runtime << " s" << endl;
}
