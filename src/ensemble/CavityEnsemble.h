
#ifndef CAVITYENSEMBLE_H_
#define CAVITYENSEMBLE_H_

#define NUM_THRESHOLD 14
#define MIN_THRESHOLD 3
#define MAX_THRESHOLD 3000

#include <map>

#include "utils/Random.h"
#include "molecules/Molecule.h"
#include "molecules/Component.h"

using namespace std;

class DomainDecompBase;

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

   void setPotential(double u) { this->insertionPotential = u; }
   void setMaxNeighbours(unsigned n, double rr) { this->maxNeighbours = n; this->r2n = rr; if(0.0 >= r2vicinity) r2vicinity = rr; }
   double getRR() { return this->r2n; }

   void init(Component* component, unsigned Nx, unsigned Ny, unsigned Nz);
   void preprocessStep();
   bool decideActivity(unsigned neighbours, unsigned long tmid);
   bool decideActivity(double uPotTilde, unsigned long tmid);
   
   void setIdOffset(unsigned long offset) { this->idoffset = offset; }
   unsigned long communicateNumCavities(DomainDecompBase* comm);
   unsigned long numCavities() { return this->globalActive; }

   map<unsigned long, Molecule*>* particleContainer() { return &(this->reservoir); }
   map<unsigned long, Molecule*> activeParticleContainer();
   
   map<unsigned long, Molecule*> exportBottom(int d);
   map<unsigned long, Molecule*> exportTop(int d);
   
   map<unsigned long, unsigned long> exportClusterBottom(int d);
   map<unsigned long, unsigned long> exportClusterTop(int d);

   void haloClear();  // deletes all halo cavities and clears the halo map
   void haloInsert(Molecule* m, bool active);
   void haloCluster(unsigned long molid, unsigned long clusterid);
   
   void determineBoundary();
   void processBoundary();
   void detectClusters();
   void connect(unsigned long i, unsigned long j);
   void attach(unsigned long vertex, unsigned long cluster, bool halo, bool detach_previous);
   
   unsigned getThreshold(unsigned i) { return ((i < NUM_THRESHOLD)? globalThreshold[i]: 0); }
   unsigned getThresholdPopulation(unsigned i) { return ((i < NUM_THRESHOLD)? globalThresholdPopulation[i]: 0); }
   unsigned getLargestCavity() { return this->globalMaxClusterSize; }
   
   bool isUltra() { return this->ultra; }
   void enableUltra() { this->ultra = true; }

 private:
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
   
   set<unsigned long> activeHalo;
   map<unsigned long, Molecule*> reservoirHalo;
   
   bool boundarySpecified;
   double init_bottom[3];
   double init_top[3];
   set<unsigned long> export_bottom[3];  // which molecule IDs to export down in dimension d
   set<unsigned long> export_top[3];   // which molecule IDs to export up in dimension d

   double insertionPotential;
   unsigned maxNeighbours;
   double r2n;

   Random async;

   double r2vicinity;
   map< unsigned long, set<unsigned long> > vicinity;  // for all active cavities, contains the active neighbour cavities (INCLUDING cavity halo!)

   map<unsigned long, unsigned long> clusterID;  // contains cluster ID of all active cavities (i.e. smallest cavity ID from the cluster, INCLUDING halo)
   map< unsigned long, set<unsigned long> > clusterVertices;  // all (local) cavities within the cluster (INCLUDING cavity halo!)
   map<unsigned long, unsigned> localClusterSize;  // number of cavities within the cluster (EXCLUDING cavity halo!)

   unsigned globalMaxClusterSize;
   map<unsigned, unsigned> globalExactPopulation;  // number of clusters with exactly N cavities
   unsigned globalThreshold[NUM_THRESHOLD];
   unsigned globalThresholdPopulation[NUM_THRESHOLD];  // number of clusters with at least globalThreshold[N] cavities
   
   bool ultra;  // deactivate all unnecessary cluster analysis -> for ultra-large simulations
};

#endif
