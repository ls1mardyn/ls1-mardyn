#include "integrators/Integrator.h"
#include <iostream>

utils::Log integrators::Integrator::_log("Integrator");

integrators::Integrator::Integrator(){
}

integrators::Integrator::~Integrator(){
  _log.error("removeAllMolecules()","Has to be implemented in derived class");
}


//void integrators::Integrator::registerInformation(std::string information, PhaseSpace* phaseSpace){
//  _log.error("removeAllMolecules()","Has to be implemented in derived class");
//}


//void integrators::Integrator::euler_vel_acc(std::list<SimpleMolecule>& molecules, double timesteplength){
//
//  std::list<SimpleMolecule>::iterator mol_iter;
//  for(mol_iter=molecules.begin(); mol_iter!=molecules.end(); mol_iter++){
//    for(int i=0; i<3; i++){
//      // a = F/M
//      // vnew = vold + dt*a
//      // rnew = roel + dt*v
//      mol_iter->setVelocity(mol_iter->getVelocity(i) + mol_iter->getForce(i)*timesteplength,i);
//      mol_iter->setPosition(mol_iter->getPosition(i) + mol_iter->getVelocity(i)*timesteplength,i);
//    }
//  }
//}
