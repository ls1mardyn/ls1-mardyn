
#ifndef CAVITYENSEMBLE_H_
#define CAVITYENSEMBLE_H_

#include <map>

#include "utils/Random.h"
#include "molecules/Molecule.h"
#include "molecules/Component.h"

using namespace std;

class DomainDecompBase;
class ParticleContainer;

class CavityEnsemble {
 public:
   CavityEnsemble();

   unsigned getInterval() { return this->interval; }
   void setInterval(unsigned delta) { this->interval = delta; }
   void setSystem(double x, double y, double z);
   double systemSize(int d) { return this->system[d]; }
   
   void setSubdomain(int rank, double x0, double x1, double y0, double y1, double z0, double z1);
   void submitTemperature(double T_in) { this->T = T_in; }
   void setControlVolume(
      double x0, double y0, double z0, double x1, double y1, double z1
   );

   void setMaxNeighbours(unsigned n, double rr) { this->maxNeighbours = n; this->r2n = rr; }
   double getRR() const { return this->r2n; }

   void init(Component* component, unsigned Nx, unsigned Ny, unsigned Nz);
   void preprocessStep();
   bool decideActivity(unsigned neighbours, unsigned long tmid);
   bool decideActivity(double uPotTilde, unsigned long tmid);
   
   void setIdOffset(unsigned long offset) { this->idoffset = offset; }
   unsigned long communicateNumCavities(DomainDecompBase* comm);
   unsigned long numCavities() { return this->globalActive; }

   map<unsigned long, Molecule*>* particleContainer() { return &(this->reservoir); }
   map<unsigned long, Molecule*> activeParticleContainer();
   
   void determineBoundary();
   void processBoundary();

   void cavityStep(ParticleContainer * globalMoleculeContainer);

 private:
   /**
    * iterate a box of center this molecule and side-lengths two times the search-radius
    * and count the neighbours contained in it.
    */
   unsigned countNeighbours(ParticleContainer * container, Molecule * m1) const;

   int ownrank;  // for debugging purposes (indicate rank in console output)
   bool initialized;
   bool rotated;
   unsigned interval;  // how often?

   unsigned componentid;
   double T;

   double globalV;
   double system[3];  // extent of the system
   float minredco[3];  // minimal coordinates of the subdomain reduced w. r. t. the system size
   float maxredco[3];   // maximal coordinates of the subdomain reduced w. r. t. the system size

   bool restrictedControlVolume;
   double control_bottom[3];
   double control_top[3];

   unsigned long idoffset;
   set<unsigned long> active;
   map<unsigned long, Molecule*> reservoir;
   unsigned long globalActive;
   
   bool boundarySpecified;
   double init_bottom[3];
   double init_top[3];

   unsigned maxNeighbours;
   double r2n;

   Random async;
};

#endif
