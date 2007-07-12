#include "parallel/DomainDecompDummy.h"
#include "molecules/Molecule.h"
#include "datastructures/ParticleContainer.h"
#include <iostream>
#include <ctime>
//#include <list>

using namespace std;

utils::Log parallel::DomainDecompDummy::_log("DomainDecompDummy");

parallel::DomainDecompDummy::DomainDecompDummy(){
}

parallel::DomainDecompDummy::~DomainDecompDummy(){
}

const char* parallel::DomainDecompDummy::getProcessorName() const{
  _log.error("getProcessorName()", "This method is not implemented yet");
  return "main";
}

void parallel::DomainDecompDummy::exchangeMolecules(datastructures::ParticleContainer<Molecule>* moleculeContainer, 
                                 const vector<Component>& components, Domain* domain){
  
  double rmin[3]; // lower corner of the process-specific domain //PARALLEL
  double rmax[3];
  double halo_L[3]; // width of the halo strip //PARALLEL
  for(int i=0; i<3; i++){
    rmin[i] =  moleculeContainer->getBoundingBoxMin(i);
    rmax[i] =  moleculeContainer->getBoundingBoxMax(i);
    halo_L[i] =  moleculeContainer->get_halo_L(i);
  }

  Molecule* currentMolecule;
  // molecules that have to be copied (because of halo), get a new position
  double new_position[3];

  double phaseSpaceSize[3];

  double low_limit; // particles below this limit have to be copied or moved to the lower process 
  double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process 

  for(unsigned short d=0;d<3;++d)
  {
    phaseSpaceSize[d] = rmax[d]-rmin[d];
    // set limits (outside "inner" region) 
    low_limit = rmin[d]+halo_L[d];
    high_limit = rmax[d]-halo_L[d];
    currentMolecule = moleculeContainer->begin();
    
    //cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << endl;
    //cout << "halo_L: " << halo_L[0] << " / " << halo_L[1] << " / " << halo_L[2] << endl;
    //cout << "proc_domain_L: " << proc_domain_L[0] << " / " << proc_domain_L[1] << " / " << proc_domain_L[2] << endl;
    while(currentMolecule!=moleculeContainer->end()){
      
      const double& rd=currentMolecule->r(d);
      if(rd<low_limit){ 
        // determine the position for the copy of the molecule
        for(unsigned short d2=0; d2<3; d2++){
          // when moving parallel to the coordinate d2 to another process, the
          // local coordinates in d2 change
          if(d2==d) new_position[d2] = rd+phaseSpaceSize[d2];
          else new_position[d2] = currentMolecule->r(d2);
        }

        Molecule m1 = Molecule(currentMolecule->id(),currentMolecule->componentid(),
                                       new_position[0], new_position[1], new_position[2],
                                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2), 
                                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
                                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2), &components);
        moleculeContainer->addParticle(m1);
        currentMolecule = moleculeContainer->next();
      }
      else if(rd>=high_limit){
        // determine the position for the copy of the molecule
        for(unsigned short d2=0; d2<3; d2++){
          // when moving parallel to the coordinate d2 to another process, the
          // local coordinates in d2 change
          if(d2==d) new_position[d2] = rd-phaseSpaceSize[d2];
          else new_position[d2] = currentMolecule->r(d2);
        }    
        Molecule m1 = Molecule(currentMolecule->id(),currentMolecule->componentid(),
                                       new_position[0], new_position[1], new_position[2],
                                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2), 
                                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
                                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2), &components);      
        moleculeContainer->addParticle(m1);                             
        currentMolecule = moleculeContainer->next();
      }
      else currentMolecule = moleculeContainer->next();
    }
  }
}

//void parallel::DomainDecompDummy::gather_molecules(list<Molecule>& m_molecules, list<Molecule>& allMolecules, 
//      const double *proc_domain_L, const double *m_rmin, const vector<Component>& components){
//}

void parallel::DomainDecompDummy::writeMoleculesToFile(string filename, datastructures::ParticleContainer<Molecule>* moleculeContainer){
 
  ofstream checkpointfilestream(filename.c_str(), ios::app);
  
  Molecule* tempMolecule;
  for(tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()){
    tempMolecule->write(checkpointfilestream);
  }    
  checkpointfilestream.close();
}


void parallel::DomainDecompDummy::reducevalues(double *Upot, double *Virial, double *summv2, double *sumIw2, unsigned long *num_mol) {
}

int parallel::DomainDecompDummy::getRank(int x, int y, int z) {
  return 0;
}

int parallel::DomainDecompDummy::getRank() {
  return 0;
}

int parallel::DomainDecompDummy::getGridSize(int dimension) {
  return 1;
}
 
int parallel::DomainDecompDummy::getCoords(int dimension) {
  return 0;
}

void parallel::DomainDecompDummy::barrier() {
}

void parallel::DomainDecompDummy::setGridSize(int num_procs){
}

double parallel::DomainDecompDummy::getTime(){
  return double(clock())/CLOCKS_PER_SEC;
}
