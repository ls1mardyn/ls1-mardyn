#include "DomainDecompDummy.h"

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "ensemble/CavityEnsemble.h"

using namespace std;
using Log::global_log;

DomainDecompDummy::DomainDecompDummy() {
}

DomainDecompDummy::~DomainDecompDummy() {
}

void DomainDecompDummy::readXML(XMLfileUnits& xmlconfig) {
	/* no parameters */
}

void DomainDecompDummy::exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain) {

	double rmin[3]; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax[3];
	double halo_L[3]; // width of the halo strip //ENABLE_MPI
	for (int i = 0; i < 3; i++) {
		rmin[i] = moleculeContainer->getBoundingBoxMin(i);
		rmax[i] = moleculeContainer->getBoundingBoxMax(i);
		halo_L[i] = moleculeContainer->get_halo_L(i);
	}

	Molecule* currentMolecule;
	// molecules that have to be copied (because of halo), get a new position
	double new_position[3];

	double phaseSpaceSize[3];

	double low_limit; // particles below this limit have to be copied or moved to the lower process
	double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process

	for (unsigned short d = 0; d < 3; ++d) {
		phaseSpaceSize[d] = rmax[d] - rmin[d];
		// set limits (outside "inner" region)
		low_limit = rmin[d] + halo_L[d];
		high_limit = rmax[d] - halo_L[d];
		currentMolecule = moleculeContainer->begin();

		//cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << endl;
		//cout << "halo_L: " << halo_L[0] << " / " << halo_L[1] << " / " << halo_L[2] << endl;
		//cout << "proc_domain_L: " << proc_domain_L[0] << " / " << proc_domain_L[1] << " / " << proc_domain_L[2] << endl;
		while (currentMolecule != moleculeContainer->end()) {
			const double& rd = currentMolecule->r(d);
			if (rd < low_limit) {
				// determine the position for the copy of the molecule
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					// when moving parallel to the coordinate d2 to another process, the
					// local coordinates in d2 change
					if (d2 == d)
						new_position[d2] = rd + phaseSpaceSize[d2];
					else
						new_position[d2] = currentMolecule->r(d2);
				}
				Component* component = _simulation.getEnsemble()->component(currentMolecule->componentid());
				Molecule m1 = Molecule(currentMolecule->id(),component,
				                       new_position[0], new_position[1], new_position[2],
				                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
				                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
				                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2));
				moleculeContainer->addParticle(m1);
				currentMolecule = moleculeContainer->next();
			}
			else if (rd >= high_limit) {
				// determine the position for the copy of the molecule
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					// when moving parallel to the coordinate d2 to another process, the
					// local coordinates in d2 change
					if (d2 == d)
						new_position[d2] = rd - phaseSpaceSize[d2];
					else
						new_position[d2] = currentMolecule->r(d2);
				}
				Component* component = _simulation.getEnsemble()->component(currentMolecule->componentid());
				Molecule m1 = Molecule(currentMolecule->id(),component,
				                       new_position[0], new_position[1], new_position[2],
				                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
				                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
				                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2));
				moleculeContainer->addParticle(m1);
				currentMolecule = moleculeContainer->next();
			}
			else
				currentMolecule = moleculeContainer->next();
		}
	}
}

void DomainDecompDummy::exchangeCavities(CavityEnsemble* cavityEnsemble, Domain* domain)
{
   cavityEnsemble->haloClear();
   
   Molecule* currentMolecule;
   Molecule* newMolecule;
   double new_position[3];
   
   map<unsigned long, Molecule*> texp;
   map<unsigned long, unsigned long> texp_cluster;
   map<unsigned long, Molecule*>::iterator texit;
   
   cavityEnsemble->processBoundary();
   
   for(int d = 0; d < 3; d++)
   {
      texp = cavityEnsemble->exportBottom(d);
      texp_cluster = cavityEnsemble->exportClusterBottom(d);
      
      for(texit = texp.begin(); texit != texp.end(); texit++)
      {
         currentMolecule = texit->second;
         assert(currentMolecule->id() == texit->first);
         for (int d2 = 0; d2 < 3; d2++)
         {
            // For DomainDecompDummy (sequential linked cells), the position ALWAYS needs to be shifted!
            //
            if(d2 == d) new_position[d2] = currentMolecule->r(d2) + cavityEnsemble->systemSize(d2);
            else new_position[d2] = currentMolecule->r(d2);
         }
         Component* component = _simulation.getEnsemble()->component(currentMolecule->componentid());
         newMolecule = new Molecule(
            texit->first, component,
            new_position[0], new_position[1], new_position[2],
            currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
            currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
            currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2)
         );
         cavityEnsemble->haloInsert(newMolecule, true);
         if(texp_cluster.count(texit->first) > 0)
         {
            cavityEnsemble->haloCluster(texit->first, texp_cluster[texit->first]);
         }
      }
      
      texp = cavityEnsemble->exportTop(d);
      texp_cluster = cavityEnsemble->exportClusterTop(d);
      
      for(texit = texp.begin(); texit != texp.end(); texit++)
      {
         currentMolecule = texit->second;
         assert(currentMolecule->id() == texit->first);
         for (int d2 = 0; d2 < 3; d2++)
         {
            // For DomainDecompDummy (sequential linked cells), the position ALWAYS needs to be shifted!
            //
            if(d2 == d) new_position[d2] = currentMolecule->r(d2) - cavityEnsemble->systemSize(d2);
            else new_position[d2] = currentMolecule->r(d2);
         }
         Component* component = _simulation.getEnsemble()->component(currentMolecule->componentid());
         Molecule* newMolecule = new Molecule(
            texit->first, component,
            new_position[0], new_position[1], new_position[2],
            currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
            currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
            currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2)
         );
         cavityEnsemble->haloInsert(newMolecule, true);
         if(texp_cluster.count(texit->first) > 0)
         {
            cavityEnsemble->haloCluster(texit->first, texp_cluster[texit->first]);
         }
      }
   }
}

unsigned DomainDecompDummy::gatherClusters(map<unsigned long, unsigned>* localClusterSize, map<unsigned, unsigned>* globalSizePopulation)
{
   unsigned max_grid_points = 0;
   
   (*globalSizePopulation) = map<unsigned, unsigned>();
   map<unsigned long, unsigned>::iterator lcsit;
   for(lcsit = localClusterSize->begin(); lcsit != localClusterSize->end(); lcsit++)
   {
      if(lcsit->second > max_grid_points) max_grid_points = lcsit->second;
      
      if(globalSizePopulation->count(lcsit->second) == 0)
      {
         (*globalSizePopulation)[lcsit->second] = 1;
      }
      else (*globalSizePopulation)[lcsit->second] ++;
   }
   
   return max_grid_points;
}

void DomainDecompDummy::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, Domain* domain) {
	exchangeMolecules(moleculeContainer, domain);
}


double DomainDecompDummy::getBoundingBoxMin(int dimension, Domain* domain) {
	return 0.0;
}
double DomainDecompDummy::getBoundingBoxMax(int dimension, Domain* domain) {
	return domain->getGlobalLength(dimension);
}


double DomainDecompDummy::getTime() {
	return double(clock()) / CLOCKS_PER_SEC;
}

unsigned DomainDecompDummy::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	*minrnd = 0.0;
	*maxrnd = 1.0;
	return localN;
}

void DomainDecompDummy::assertIntIdentity(int IX) {
}

void DomainDecompDummy::assertDisjunctivity(TMoleculeContainer* mm) {
}

void DomainDecompDummy::printDecomp(std::string filename, Domain* domain) {
	global_log->warning() << "printDecomp useless in serial mode" << endl;
}
