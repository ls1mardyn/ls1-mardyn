#include "tests/MoleculeTest.h"
#include "molecules/SimpleMolecule.h"

#include <cmath>
#include <iostream>
using namespace std;

utils::Log MoleculeTest::_log("moleculeTest");

MoleculeTest::MoleculeTest():
  TestCase( "moleculeTest" ) {
}

MoleculeTest::~MoleculeTest(){
}


void MoleculeTest::get_and_set_Test(){

  // Test Constructor and get
  SimpleMolecule molecule1(4,5,0.01,0.02,0.03,0.11,0.12,0.13);
  
  validateEquals(molecule1.getId(),4, "get_and_set_Test()");
  validateEquals(molecule1.getType(),5, "get_and_set_Test()");
  validateEquals(molecule1.getPosition(0),0.01, "get_and_set_Test()");
  validateEquals(molecule1.getPosition(1),0.02, "get_and_set_Test()");
  validateEquals(molecule1.getPosition(2),0.03, "get_and_set_Test()");
  validateEquals(molecule1.getVelocity(0),0.11, "get_and_set_Test()");
  validateEquals(molecule1.getVelocity(1),0.12, "get_and_set_Test()");
  validateEquals(molecule1.getVelocity(2),0.13, "get_and_set_Test()");
  
  // Test set and get
  molecule1.setId(13);
  molecule1.setType(8);
  molecule1.setPosition(3.85,0);
  molecule1.setPosition(5.10,1);
  molecule1.setPosition(-3.94,2);
  molecule1.setVelocity(1.2,0);
  molecule1.setVelocity(-12.3,1);
  molecule1.setVelocity(8.0,2);
  molecule1.setForce(-0.01,0);
  molecule1.setForce(8.11,1);
  molecule1.setForce(-2.0,2);

  validateEquals(molecule1.getId(),13, "get_and_set_Test()");
  validateEquals(molecule1.getType(),8, "get_and_set_Test()");
  validateEquals(molecule1.getPosition(0),3.85, "get_and_set_Test()");
  validateEquals(molecule1.getPosition(1),5.10, "get_and_set_Test()");
  validateEquals(molecule1.getPosition(2),-3.94, "get_and_set_Test()");
  validateEquals(molecule1.getVelocity(0),1.2, "get_and_set_Test()");
  validateEquals(molecule1.getVelocity(1),-12.3, "get_and_set_Test()");
  validateEquals(molecule1.getVelocity(2),8.0, "get_and_set_Test()");
  validateEquals(molecule1.getForce(0),-0.01, "get_and_set_Test()");
  validateEquals(molecule1.getForce(1),8.11, "get_and_set_Test()");
  validateEquals(molecule1.getForce(2),-2.0, "get_and_set_Test()");
}


void MoleculeTest::distance_calc_Test(){
  // create two molecules and test whether the distance square is correctly calculatet
  SimpleMolecule molecule1(1,3,3.01,2.42,-1.6,0.11,0.12,0.13);
  SimpleMolecule molecule2(2,5,2.00,2.51,0.1,1.11,2.12,3.13);
  validateEquals(SimpleMolecule::calcDistanceSquare(molecule1,molecule2),3.9182, "distance_calc_Test()");
}

void MoleculeTest::force_calc_Test(){
  // create two molecules and test whether the force is correctly calculated
  SimpleMolecule molecule1(1,3,3.01,2.42,-1.6,0.11,0.12,0.13);
  SimpleMolecule molecule2(2,5,2.00,2.51,0.1,1.11,2.12,3.13);
  
  // Set the forces on both molecules to zero
  molecule1.setForce(0,0);
  molecule1.setForce(0,1);
  molecule1.setForce(0,2);
  molecule2.setForce(0,0);
  molecule2.setForce(0,1);
  molecule2.setForce(0,2);
  
  // squared distance
  double distance_square = 3.9182;
  
  // correct calculation of (Fi/rij)
  double Fi_by_rij = 24 * 5 * 1/distance_square * pow((1/distance_square),3) 
                     * (1-2*pow((1/distance_square),3));

  // calculate forces
  SimpleMolecule::calcForce(molecule1,molecule2);
  
  // Force on molecule1
  double F_m1_x = molecule1.getForce(0);
  double F_m1_y = molecule1.getForce(1);
  double F_m1_z = molecule1.getForce(2);
  double F_m2_x = molecule2.getForce(0);
  double F_m2_y = molecule2.getForce(1);
  double F_m2_z = molecule2.getForce(2);
  

  // compare calculated forces with correct forces
  // (Fi/rij) has to be multiplied with the distance in the corresponding dimension
  validateNumericalEquals(F_m1_x,Fi_by_rij*-1.01, "force_calc_Test()");
  validateNumericalEquals(F_m1_y,Fi_by_rij*0.09, "force_calc_Test()");
  validateNumericalEquals(F_m1_z,Fi_by_rij*1.7, "force_calc_Test()");
  validateNumericalEquals(F_m2_x,Fi_by_rij*1.01, "force_calc_Test()");
  validateNumericalEquals(F_m2_y,Fi_by_rij*-0.09, "force_calc_Test()");
  validateNumericalEquals(F_m2_z,Fi_by_rij*-1.7, "force_calc_Test()");
}


void MoleculeTest::run(){
  get_and_set_Test();
  distance_calc_Test();
  force_calc_Test();
}
