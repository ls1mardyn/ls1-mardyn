#include <iostream>
#include "Simulation.h"

using namespace std;

//! @page main
//! In this project, software for the simulation of Molecular Dynamics
//! with short-range forces is developed. The aim is to have a parallel code (MPI) 
//! for multi-centered molecules.
//!
//! The role of the main function is to run Tests for all classes
//! and to instantiate an Object of the Simulation class which
//! then is responsible for the simulation
//!
int main(int argc, char** argv){

  cout.precision(10);

  Simulation simulation(&argc, &argv);
  
  double runtime = double(clock())/CLOCKS_PER_SEC;

  simulation.simulate();

  runtime=double(clock())/CLOCKS_PER_SEC-runtime;
  
  cout << "main: used " << runtime << " s" << endl;

}
