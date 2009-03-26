#include "datastructures/ParticleContainer.h"

#include <iostream>
using namespace std;

ParticleContainer::ParticleContainer(ParticlePairsHandler& partPairsHandler, double bBoxMin[3], double bBoxMax[3] ):
  _particlePairsHandler(partPairsHandler){
  for(int i=0; i<3; i++){
    _boundingBoxMin[i] = bBoxMin[i];
    _boundingBoxMax[i] = bBoxMax[i];
  }
}

ParticleContainer::~ParticleContainer(){
}

void ParticleContainer::rebuild(double bBoxMin[3], double bBoxMax[3]){
  cout << "REBUILD OF PARTICLE CONTAINER" << endl;
  for(int i=0; i<3; i++){
    _boundingBoxMin[i] = bBoxMin[i];
    _boundingBoxMax[i] = bBoxMax[i];
  }
}

double ParticleContainer::getBoundingBoxMin(int dimension){
  return this->_boundingBoxMin[dimension];
}

double ParticleContainer::getBoundingBoxMax(int dimension){
  return this->_boundingBoxMax[dimension];
}

double ParticleContainer::get_halo_L(int index){   
  cerr << "ERROR: ParticleContainer::get_halo_L(...) has to be implemented in derived class" << endl;
  return 0;
}


