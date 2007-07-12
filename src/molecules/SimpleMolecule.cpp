#include "molecules/SimpleMolecule.h"
#include <iostream>
#include <cmath>


SimpleMolecule::SimpleMolecule(int id, int type,
                   double xPos, double yPos, double zPos,
                   double xVel, double yVel, double zVel){
  this->id = id;
  this->type = type;
  this->position[0] = xPos;
  this->position[1] = yPos;
  this->position[2] = zPos;
  this->velocity[0] = xVel;
  this->velocity[1] = yVel;
  this->velocity[2] = zVel;
  //this->nextPotentialNeighbour = NULL;
}

int SimpleMolecule::getId(){
  return this->id;
}

int SimpleMolecule::getType(){
  return this->type;
}

double SimpleMolecule::getPosition(int dimension){
  return this->position[dimension];
}

double SimpleMolecule::getVelocity(int dimension){
  return this->velocity[dimension];
}

double SimpleMolecule::getForce(int dimension){
  return this->force[dimension];
}

void SimpleMolecule::setId(int id){
  this->id = id;
}

void SimpleMolecule::setType(int type){
  this->type = type;
}

void SimpleMolecule::setPosition(double position, int dimension){
  this->position[dimension] = position;	
}

void SimpleMolecule::setVelocity(double velocity, int dimension){
  this->velocity[dimension] = velocity;
}

void SimpleMolecule::setForce(double force, int dimension){
  this->force[dimension] = force;
}

double SimpleMolecule::calcDistanceSquare(SimpleMolecule& molecule1, SimpleMolecule& molecule2){
  double distance_square = 0;
  for(int i=0; i<3; i++){
    distance_square += pow(molecule1.position[i]-molecule2.position[i],2);
  }
  return distance_square;
}

void SimpleMolecule::calcForce(SimpleMolecule& molecule1, SimpleMolecule& molecule2){
  double sigma = 1.0;
  double epsilon = 5.0;
  double r_2 = calcDistanceSquare(molecule1, molecule2);
  double s = sigma*sigma / r_2;
  s = pow(s,3);
  double f = 24 * epsilon * s / r_2 * (1-2*s);
  for(int i=0; i<3; i++){
  	double temp = f*(molecule2.position[i]-molecule1.position[i]);
    molecule1.force[i] += temp;
    molecule2.force[i] -= temp;
  }
}
