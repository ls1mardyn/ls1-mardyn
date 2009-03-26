#include "integrators/Leapfrog.h"
#include "datastructures/ParticleContainer.h"
#include "Domain.h"
#include "molecules/Molecule.h"

utils::Log integrators::Leapfrog::_log("Leapfrog");

integrators::Leapfrog::Leapfrog(double timestepLength){
  // set starting state
  this->_state = 3; 
  
  this->_timestepLength = timestepLength;
}

integrators::Leapfrog::~Leapfrog(){
}


void integrators::Leapfrog::eventForcesCalculated(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state == 2){
    transition2to3(molCont, domain);
  }
}

void integrators::Leapfrog::eventNewTimestep(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state == 3){
    transition3to1(molCont, domain);
    transition1to2(molCont, domain);
  }
}


void integrators::Leapfrog::transition1to2(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state==1){
    Molecule* tempMolecule;
    double vcorr=2.-1./domain->getGlobalBetaTrans();
    double Dcorr=2.-1./domain->getGlobalBetaRot();
    for(tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()){
      tempMolecule->upd_preF(_timestepLength, vcorr, Dcorr);
    }

    this->_state = 2;
  }
  else { _log.error("transition1to2(...)","Wrong state for state transition");}
}


void integrators::Leapfrog::transition2to3(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state==2){
    Molecule* tempMolecule;
    double summv2 = 0.0;
    double sumIw2 = 0.0;
    double dt_halve=.5*_timestepLength;
    
    for(tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()){
      tempMolecule->upd_postF(dt_halve,summv2,sumIw2); 
    }
    
    domain->setLocalSummv2(summv2);
    domain->setLocalSumIw2(sumIw2);
    
    this->_state = 3;
  }
  else { _log.error("transition2to3(...)","Wrong state for state transition");}
}


void integrators::Leapfrog::transition3to1(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state==3){ this->_state = 1;}
  else{ _log.error("transition3to1(...)","Wrong state for state transition");}
}
