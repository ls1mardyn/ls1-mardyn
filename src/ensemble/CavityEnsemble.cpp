
#include "CavityEnsemble.h"

#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"
#include "molecules/Quaternion.h"

#define COMMUNICATION_THRESHOLD 3

using Log::global_log;

CavityEnsemble::CavityEnsemble()
{
   this->ownrank = -1;
   this->initialized = false;
   this->rotated = false;
   this->interval = 1000;
   this->componentid = (unsigned)-1;
   this->T = -1.0;

   this->globalV = 0.0;
   this->system[0] = 0.0; system[1] = 0.0; system[2] = 0.0;
   this->minredco[0] = 0.0; minredco[1] = 0.0; minredco[2] = 0.0;
   this->maxredco[0] = 1.0; maxredco[1] = 1.0; maxredco[2] = 1.0;

   this->restrictedControlVolume = false;
   this->control_bottom[0] = 0.0; control_bottom[1] = 0.0; control_bottom[2] = 0.0;
   this->control_top[0] = 1.0; control_top[1] = 1.0; control_top[2] = 1.0;

   this->active = set<unsigned long>();
   this->reservoir = map<unsigned long, Molecule*>();
   this->globalActive = 0;

   this->boundarySpecified = false;
   for(int d = 0; d < 3; d++)
   {
      this->export_bottom[d] = set<unsigned long>();
      this->export_top[d] = set<unsigned long>();
   }
   
   this->insertionPotential = -1.0e+10;
   this->maxNeighbours = 1024;
   this->r2n = 0.0;
   this->idoffset = 0;

   this->r2vicinity = -1.0;
   this->vicinity = map< unsigned long, set<unsigned long> >();

   this->clusterID = map<unsigned long, unsigned long>();
   this->clusterVertices = map< unsigned long, set<unsigned long> >();
   this->localClusterSize = map<unsigned long, unsigned>();

   this->globalMaxClusterSize = 0;
   this->globalExactPopulation = map<unsigned, unsigned>();
   double thrfactor = pow((double)MAX_THRESHOLD/MIN_THRESHOLD, 1.0/(NUM_THRESHOLD - 1.0));
   double tthreshold = MIN_THRESHOLD;
   for(unsigned i=0; i < NUM_THRESHOLD; i++)
   {
      this->globalThreshold[i] = (unsigned)round(tthreshold);
      // cout << "globalThreshold[" << i << "] = " << this->globalThreshold[i] << "\n";
      this->globalThresholdPopulation[i] = 0;
      tthreshold *= thrfactor;
   }
   this->ultra = false;
}

void CavityEnsemble::setSystem(double x, double y, double z)
{
   this->system[0] = x; this->system[1] = y; this->system[2] = z;
   if(!this->restrictedControlVolume)
   {
      this->globalV = x*y*z;
	    
      for(int d=0; d < 3; d++)
      {
         this->control_bottom[d] = 0.0;
         this->control_top[d] = this->system[d];
      }
   }
}

void CavityEnsemble::setSubdomain(int rank, double x0, double x1, double y0, double y1, double z0, double z1)
{
   this->ownrank = rank;
   this->async.init(8623);

   for(int d=0; d < 3; d++)
      if(control_bottom[d] >= control_top[d])
         this->restrictedControlVolume = false;

   if(!this->restrictedControlVolume)
   {
      this->globalV = this->system[0] * system[1] * system[2];
    
      for(int d=0; d < 3; d++)
      {
         this->control_bottom[d] = 0.0;
         this->control_top[d] = this->system[d];
      }
   }
   this->minredco[0] = (x0 - control_bottom[0]) /
                          (control_top[0] - control_bottom[0]);
   this->minredco[1] = (y0 - control_bottom[1]) /
                          (control_top[1] - control_bottom[1]);
   this->minredco[2] = (z0 - control_bottom[2]) /
                          (control_top[2] - control_bottom[2]);
   this->maxredco[0] = (x1 - control_bottom[0]) /
                          (control_top[0] - control_bottom[0]);
   this->maxredco[1] = (y1 - control_bottom[1]) /
                          (control_top[1] - control_bottom[1]);
   this->maxredco[2] = (z1 - control_bottom[2]) /
                          (control_top[2] - control_bottom[2]);
}

void CavityEnsemble::setControlVolume(double x0, double y0, double z0, double x1, double y1, double z1)
{
   if((x0 >= x1) || (y0 >= y1) || (z0 >= z1))
   {
      global_log->error() << "\nInvalid control volume (" << x0 << " / " << y0 
                          << " / " << z0 << ") to (" << x1 << " / " << y1 << " / "
                          << z1 << ")." << std::endl;
			exit(711);
   }
	 
   this->restrictedControlVolume = true;
   this->globalV = (x1-x0) * (y1-y0) * (z1-z0);
   this->control_bottom[0] = x0;
   this->control_top[0]    = x1;
   this->control_bottom[1] = y0;
   this->control_top[1]    = y1;
   this->control_bottom[2] = z0;
   this->control_top[2]    = z1;
}

void CavityEnsemble::init(Component* component, unsigned Nx, unsigned Ny, unsigned Nz)
{
   if(this->ownrank < 0)
   {
      global_log->error() << "\nInvalid rank " << ownrank << ".\n";
      exit(712);
   }
   if(this->initialized)
   {
      global_log->error() << "\nCavity ensemble initialized twice.\n";
      exit(713);
   }
   if(0.0 >= this->T)
   {
      global_log->error() << "\nInvalid temperature T = " << T << ".\n";
      exit(714);
   }
   if(0.0 >= this->globalV)
   {
      global_log->error() << "\nInvalid control volume V_ctrl = " << globalV << ".\n";
      exit(715);
   }

   this->componentid = component->ID();
   unsigned rotdof = component->getRotationalDegreesOfFreedom();
   
   double nun[3];
   nun[0] = Nx; nun[1] = Ny; nun[2] = Nz;

   // reduced lattice in each dimension given by epsilon + 0.5/N, epsilon + 1.5/N, ..., epsilon + (N-0.5)/N
   unsigned minlu[3]; unsigned maxlu[3];
   for(int d = 0; d < 3; d++)
   {
      minlu[d] = (unsigned)round(nun[d] * this->minredco[d] - 0.0009756);
      maxlu[d] = (unsigned)round(nun[d] * this->maxredco[d] - 0.0009756) - 1;
   }
   /*
   cout << "nun: " << nun[0] << " / " << nun[1] << " / " << nun[2] << "\n";
   cout << "minredco: " << minredco[0] << " / " << minredco[1] << " / " << minredco[2] << "\n";
   cout << "minlu: " << minlu[0] << " / " << minlu[1] << " / " << minlu[2] << "\n";
   cout << "maxredco: " << maxredco[0] << " / " << maxredco[1] << " / " << maxredco[2] << "\n";
   cout << "maxlu: " << maxlu[0] << " / " << maxlu[1] << " / " << maxlu[2] << "\n";
   */

   double max_spacing = 0.0;
   double grid_spacing[3];
   for(int d = 0; d < 3; d++)
   {
      grid_spacing[d] = (control_top[d] - control_bottom[d])/nun[d];
      if(grid_spacing[d] > max_spacing) max_spacing = grid_spacing[d];
   }
   this->r2vicinity = 2.001 * max_spacing * max_spacing;
   if(!this->ownrank)
   {
      cout << "Cavities closer than " << sqrt(r2vicinity) << " are regarded as neighbours.\n";
   }

   unsigned tid;
   unsigned tlu[3];
   double tq[3];
   Molecule* tm;
   for(tlu[0] = minlu[0]; maxlu[0] >= tlu[0]; tlu[0]++)
   {
      for(tlu[1] = minlu[1]; maxlu[1] >= tlu[1]; tlu[1]++)
      {
         for(tlu[2] = minlu[2]; maxlu[2] >= tlu[2]; tlu[2]++)
         {
            tid = tlu[0]*Ny*Nz + tlu[1]*Nz + tlu[2] + 1; // + this->idoffset;
            for(int d=0; d < 3; d++)
               tq[d] = control_bottom[d] + (0.5009756 + tlu[d])*grid_spacing[d];

            double v[3];
            double vv = 0.0;
            for(int d=0; d < 3; d++)
            {
               v[d] = -0.5 + this->async.rnd();
               vv += v[d] * v[d];
            }
            double vnorm = sqrt(3.0*T / (vv * component->m()));

            double qtr[4];
            double qqtr = 0.0;
            for(int d=0; d < 4; d++)
            {
               qtr[d] = -0.5 + this->async.rnd();
               qqtr += qtr[d] * qtr[d];
            }
            double qtrnorm = sqrt(1.0 / qqtr);

            double D[3];
            double Dnorm = 0.0;
            if(rotdof > 0)
            {
               for(int d=0; d < 3; d++) D[d] = -0.5 + this->async.rnd();

               double w[3];
               Quaternion tqtr = Quaternion(qtr[0]*qtrnorm, qtr[1]*qtrnorm, qtr[2]*qtrnorm, qtr[3]*qtrnorm);
               tqtr.rotate(D, w);
               double Iw2 = w[0]*w[0]*component->I11()
                          + w[1]*w[1]*component->I22()
                          + w[2]*w[2]*component->I33();

               Dnorm = sqrt(T*rotdof / Iw2);
            }
            else
            {
               D[0] = 0.0; D[1] = 0.0; D[2] = 0.0;
            }

            tm = new Molecule(tid, component, tq[0], tq[1], tq[2], v[0]*vnorm, v[1]*vnorm, v[2]*vnorm, qtr[0]*qtrnorm, qtr[1]*qtrnorm, qtr[2]*qtrnorm, qtr[3]*qtrnorm, D[0]*Dnorm, D[1]*Dnorm, D[2]*Dnorm);

            this->reservoir[tid] = tm;
         }
      }
   }

   this->initialized = true;
}

unsigned long CavityEnsemble::communicateNumCavities(DomainDecompBase* comm)
{
   comm->collCommInit(1);
   comm->collCommAppendUnsLong(this->active.size());
   comm->collCommAllreduceSum();
   this->globalActive = comm->collCommGetUnsLong();
   comm->collCommFinalize();
   
   if(this->ultra) return this->globalActive;
   
   map<unsigned long, unsigned> significantClusterSize;
   // cout << "\t\t\tClusters at rank " << ownrank << ":\n";
   for(map<unsigned long, unsigned>::iterator lcsit = localClusterSize.begin(); lcsit != localClusterSize.end(); lcsit++)
   {
      // cout << "\t\t\t\t" << lcsit->first << " (size " << lcsit->second << ")";
      if(lcsit->second >= COMMUNICATION_THRESHOLD)
      {
         significantClusterSize[lcsit->first] = lcsit->second;
         // cout << ": significant";
      }
      // cout << "\n";
   }
   this->globalMaxClusterSize = comm->gatherClusters(&significantClusterSize, &this->globalExactPopulation);
   
   if(this->ownrank == 0)
   {
      for(unsigned i = 0; i < NUM_THRESHOLD; i++) this->globalThresholdPopulation[i] = 0;
         
      map<unsigned, unsigned>::iterator gepit;
      for(gepit = globalExactPopulation.begin(); gepit != globalExactPopulation.end(); gepit++)
      {
         // cout << "\t\tclstat\t" << gepit->first << "\t" << gepit->second << "\n";
         
         for(unsigned i = 0; i < NUM_THRESHOLD; i++)
            if(gepit->first >= globalThreshold[i]) globalThresholdPopulation[i] += gepit->second;
      }
   }
        
   return this->globalActive;
}

void CavityEnsemble::preprocessStep()
{
   if(this->ultra && this->rotated) return;
   if(this->reservoir.size() == 0) return;
   
   map<unsigned long, Molecule*>::iterator resit = this->reservoir.begin();
   
   double qtr[4];
   Component* tc = resit->second->component();
   unsigned rotdof = tc->getRotationalDegreesOfFreedom();
   if(rotdof > 0)
   {
      for(resit = reservoir.begin(); resit != reservoir.end(); resit++)
      {
            double qqtr = 0.0;
            for(int d=0; d < 4; d++)
            {
               qtr[d] = -0.5 + this->async.rnd();
               qqtr += qtr[d] * qtr[d];
            }
            double qtrnorm = sqrt(1.0 / qqtr);
            resit->second->setq(Quaternion(qtrnorm*qtr[0], qtrnorm*qtr[1], qtrnorm*qtr[2], qtrnorm*qtr[3]));
      }
   }
   this->rotated = true;
}

bool CavityEnsemble::decideActivity(unsigned neighbours, unsigned long tmid)
{
   bool isActive = (this->maxNeighbours >= neighbours);
   if(isActive) this->active.insert(tmid);
   else this->active.erase(tmid);
   return isActive;
}

bool CavityEnsemble::decideActivity(double uPotTilde, unsigned long tmid)
{
   /*** ADD CAVITY ENSEMBLE ***/
   
   this->active.erase(tmid);
   
   /*** ADD CAVITY ENSEMBLE ***/
   
   return false;
}

map<unsigned long, Molecule*> CavityEnsemble::activeParticleContainer()
{
   map<unsigned long, Molecule*> retv;
   set<unsigned long>::iterator resit;
   for(resit = this->active.begin(); resit != active.end(); resit++)
   {
      retv[*resit] = this->reservoir[*resit];
   }
   return retv;
}

map<unsigned long, Molecule*> CavityEnsemble::exportBottom(int d)
{
   Molecule* tm;
   map<unsigned long, Molecule*> retv;
   if(ultra) return retv;
   
   set<unsigned long>::iterator resit;
   for(resit = this->export_bottom[d].begin(); resit != this->export_bottom[d].end(); resit++)
   {
      if(this->active.count(*resit) > 0)
      {
         assert(this->reservoir.count(*resit) > 0);
         tm = this->reservoir[*resit];
         if(tm != NULL) retv[*resit] = tm;
         // cout << "\t\t\t\tExport " << *resit << " (i.e. " << tm->id() << ") down in dimension no. " << d << ".\n";
      }
   }
   return retv;
}

map<unsigned long, Molecule*> CavityEnsemble::exportTop(int d)
{
   Molecule* tm;
   map<unsigned long, Molecule*> retv;
   if(ultra) return retv;
   
   set<unsigned long>::iterator resit;
   for(resit = this->export_top[d].begin(); resit != this->export_top[d].end(); resit++)
   {
      if(this->active.count(*resit) > 0)
      {
         assert(this->reservoir.count(*resit) > 0);
         tm = this->reservoir[*resit];
         if(tm != NULL) retv[*resit] = tm;
         // cout << "\t\t\t\tExport " << *resit << " (i.e. " << tm->id() << ") up in dimension no. " << d << ".\n";
      }
   }
   return retv;
}

map<unsigned long, unsigned long> CavityEnsemble::exportClusterBottom(int d)
{
   map<unsigned long, unsigned long> retv;
   if(ultra) return retv;
   
   set<unsigned long>::iterator resit;
   for(resit = this->export_bottom[d].begin(); resit != this->export_bottom[d].end(); resit++)
   {
      if(this->active.count(*resit) > 0)
      {
         retv[*resit] = (clusterID.count(*resit) > 0)? clusterID[*resit]: *resit;
      }
   }
   return retv;
}

map<unsigned long, unsigned long> CavityEnsemble::exportClusterTop(int d)
{
   map<unsigned long, unsigned long> retv;
   if(ultra) return retv;
   
   set<unsigned long>::iterator resit;
   for(resit = this->export_top[d].begin(); resit != this->export_top[d].end(); resit++)
   {
      if(this->active.count(*resit) > 0)
      {
         retv[*resit] = (clusterID.count(*resit) > 0)? clusterID[*resit]: *resit;
      }
   }
   return retv;
}

void CavityEnsemble::haloClear()
{
   this->activeHalo.clear();
   for(map<unsigned long, Molecule*>::iterator rhit = reservoirHalo.begin(); rhit != reservoirHalo.end(); rhit++)
   {
      delete rhit->second;
   }
   this->reservoirHalo.clear();
}

void CavityEnsemble::haloInsert(Molecule* m, bool active)
{
   this->reservoirHalo[m->id()] = m;
   if(active) this->activeHalo.insert(m->id());
}

void CavityEnsemble::haloCluster(unsigned long molid, unsigned long clusterid)
{
   // cout << "\t\t\t\tInserting halo vertex " << molid << " in cluster " << clusterid << ".\n";
   assert(molid >= clusterid);
   if(activeHalo.count(molid) > 0) this->clusterID[molid] = clusterid;
}

void CavityEnsemble::determineBoundary()
{
   double rvic = sqrt(r2vicinity);
   
   for(int d = 0; d < 3; d++)
   {
      this->init_bottom[d] = minredco[d]*system[d] + rvic;
      this->init_top[d] = maxredco[d]*system[d] - rvic;
      cout << "rank " << ownrank << " dim " << d << ": comm bottom at " << init_bottom[d] << " (i.e. " << minredco[d]*system[d] << " + " << rvic << "), comm top at " << init_top[d] << " (i.e. " << maxredco[d]*system[d] << " - " << rvic << ").\n";
   }
   this->boundarySpecified = true;
}

void CavityEnsemble::processBoundary()
{
   for(int d = 0; d < 3; d++)
   {
      this->export_bottom[d].clear();
      this->export_top[d].clear();
   }
   if(this->ultra) return;
   
   if(!this->boundarySpecified) this->determineBoundary();
   
   Molecule* tm;
   for(set<unsigned long>::iterator tacit = active.begin(); tacit != active.end(); tacit++)
   {
      assert(this->reservoir.count(*tacit) == 1);
      tm = this->reservoir[*tacit];
      assert(tm->id() == *tacit);
      
      if(tm != NULL)
      {
         for(int d = 0; d < 3; d++)
         {
            if(tm->r(d) < this->init_bottom[d]) this->export_bottom[d].insert(*tacit);
            if(tm->r(d) > this->init_top[d]) this->export_top[d].insert(*tacit);
         }
      }
   }
}

void CavityEnsemble::detectClusters()
{
   if(this->ultra) return;
   
   double distanceVector[3];
   set<unsigned long>::iterator acti, actj;
   unsigned long cli;
   Molecule *cavi, *cavj;
   
   unsigned long lowlink, tnode;
   set<unsigned long> processed_nodes;
   set<unsigned long> present_nodes;
   set<unsigned long> unprocessed_nodes;
   
   stack<unsigned long> dfs_stack;
   map<unsigned long, set<unsigned long>::iterator> edgeit;
   
   this->vicinity = map< unsigned long, set<unsigned long> >();
   this->clusterVertices = map< unsigned long, set<unsigned long> >();
   this->localClusterSize = map<unsigned long, unsigned>();
   for(acti = this->active.begin(); acti != this->active.end(); acti++)
   {
      this->vicinity[*acti] = set<unsigned long>();
      this->attach(*acti, *acti, false, false);
      unprocessed_nodes.insert(*acti);
   }
   for(acti = this->activeHalo.begin(); acti != this->activeHalo.end(); acti++)
   {
      this->vicinity[*acti] = set<unsigned long>();
      cli = (this->clusterID.count(*acti) > 0)? clusterID[*acti]: *acti;
      this->attach(*acti, cli, true, false);
      unprocessed_nodes.insert(*acti);
   }
   
   for(acti = this->active.begin(); acti != this->active.end(); acti++)
   {
      assert(this->reservoir.count(*acti) > 0);
      cavi = this->reservoir[*acti];
      assert(cavi != NULL);
      
      actj = acti;
      for(actj++; actj != this->active.end(); actj++)
      {
         assert(this->reservoir.count(*actj) > 0);
         cavj = this->reservoir[*actj];
         assert(cavj != NULL);
         
         double dd = cavj->dist2(*cavi, distanceVector);
         if (dd < this->r2vicinity) this->connect(*acti, *actj);
      }
      for(actj = this->activeHalo.begin(); actj != this->activeHalo.end(); actj++)
      {
         assert(this->reservoirHalo.count(*actj) > 0);
         cavj = this->reservoirHalo[*actj];
         assert(cavj != NULL);
         
         double dd = cavj->dist2(*cavi, distanceVector);
         if (dd < this->r2vicinity) this->connect(*acti, *actj);
      }
   }
   
   while(!unprocessed_nodes.empty())
   {
      present_nodes = set<unsigned long>();
      dfs_stack.push( *(unprocessed_nodes.begin()) );
      lowlink = dfs_stack.top();
      
      while(!dfs_stack.empty())
      {
         tnode = dfs_stack.top();
         if(unprocessed_nodes.count(tnode) > 0)
         {
            if(this->clusterID[tnode] < lowlink)
            {
               lowlink = this->clusterID[tnode];
            }
            unprocessed_nodes.erase(tnode);
            present_nodes.insert(tnode);
            edgeit[tnode] = this->vicinity[tnode].begin();
         }
         
         while((edgeit[tnode] != this->vicinity[tnode].end()) && (unprocessed_nodes.count(*edgeit[tnode]) == 0))
         {
            edgeit[tnode] ++;
         }
         if(edgeit[tnode] != this->vicinity[tnode].end())
         {
            tnode = *edgeit[tnode];
            dfs_stack.push(tnode);
         }
         else
         {
            dfs_stack.pop();
         }
      }
      for(acti = present_nodes.begin(); acti != present_nodes.end(); acti++)
      {
         processed_nodes.insert(*acti);
         attach(*acti, lowlink, (active.count(*acti) == 0), true);
      }
   }
}

void CavityEnsemble::connect(unsigned long i, unsigned long j)
{
   this->vicinity[i].insert(j);
   this->vicinity[j].insert(i);
   // cout << "\t\t\t\t" << "Connect " << i << " to " << j << ".\n";
}

void CavityEnsemble::attach(unsigned long vertex, unsigned long cluster, bool halo, bool detach_previous)
{
   // cout << "\t\t\t\tRank " << ownrank << " attaches vertex " << vertex << " (halo? " << halo << ") to cluster " << cluster;
   // if(detach_previous) cout << " (detaching from " << clusterID[vertex] << ")";
   assert(vertex >= cluster);

   if(detach_previous && (clusterID.count(vertex) > 0))
   {
      if(clusterVertices.count(clusterID[vertex]) > 0)
      {
         this->clusterVertices[clusterID[vertex]].erase(vertex);
      }
      if(!halo && (localClusterSize.count(clusterID[vertex]) > 0))
      {
         this->localClusterSize[this->clusterID[vertex]] --;
         if(localClusterSize[clusterID[vertex]] == 0) localClusterSize.erase(clusterID[vertex]);
      }
   }
   
   this->clusterID[vertex] = cluster;
   
   if(this->clusterVertices.count(cluster) == 0)
   {
      this->clusterVertices[cluster] = set<unsigned long>();
   }
   if(clusterVertices[cluster].count(vertex) == 0)
   {
      this->clusterVertices[cluster].insert(vertex);
      if(!halo)
      {
         if(this->localClusterSize.count(cluster) == 0)
         {
            this->localClusterSize[cluster] = 0;
         }
         this->localClusterSize[cluster] ++;
      }
   }
   // cout << ", reaching size " << localClusterSize[cluster] << " locally.\n";
}
