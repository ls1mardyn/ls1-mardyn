#include "integrators/Leapfrog.h"
#include "datastructures/ParticleContainer.h"
#include "Domain.h"
#include "molecules/Molecule.h"

Leapfrog::Leapfrog(double timestepLength){
  // set starting state
  this->_state = 3; 
  
  this->_timestepLength = timestepLength;
}

Leapfrog::~Leapfrog(){
}


void Leapfrog::eventForcesCalculated(ParticleContainer* molCont, Domain* domain){
  if(this->_state == 2){
    transition2to3(molCont, domain);
  }
}

void Leapfrog::eventNewTimestep(ParticleContainer* molCont, Domain* domain){
  if(this->_state == 3){
    transition3to1(molCont, domain);
    transition1to2(molCont, domain);
  }
}


void Leapfrog::transition1to2(ParticleContainer* molCont, Domain* domain){
  if(this->_state==1){
    Molecule* tempMolecule;
    double vcorr=2.-1./domain->getGlobalBetaTrans();
    double Dcorr=2.-1./domain->getGlobalBetaRot();
    for(tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()){
      tempMolecule->upd_preF(_timestepLength, vcorr, Dcorr);
    }
		
    this->_state = 2;
  }
  else { cerr << "Leapfrog::transition1to2(...): Wrong state for state transition" << endl;}
}


void Leapfrog::transition2to3(ParticleContainer* molCont, Domain* domain){
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
  else { cerr << "Leapfrog::transition2to3(...): Wrong state for state transition" << endl;}
}


void Leapfrog::transition3to1(ParticleContainer* molCont, Domain* domain){
  if(this->_state==3){ this->_state = 1;}
  else{ cerr << "Leapfrog::transition3to1(...): Wrong state for state transition" << endl;}
}
