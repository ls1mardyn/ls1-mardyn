/*
 * Martin Horsch, LS1/Mardyn project moderated by Martin Bernreuther
 * (C)2009 GNU General Public License
 */
#include "GrandCanonical.h"

#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"

using namespace std;
using Log::global_log;

Random::Random() { this->init(8624); }

void Random::init(int seed)
{
	 this->ix_muVT = (seed ^ (int)888889999) | (int)1;
	 this->iy_muVT = seed ^ (int)777755555;
	 // Calculate normalization factor
	 this->am_muVT = 2.0 / (1.0 + (unsigned)((int)-1));
}

float Random::rnd_muVT()
{
	 float rnd;
	 const int IA = 16807;
	 const int IM = 2147483647;
	 const int IQ = 127773; 
	 const int IR = 2836;

	 ix_muVT ^= (ix_muVT >> 13);
	 ix_muVT ^= (ix_muVT << 17);
	 ix_muVT ^= (ix_muVT >> 5);

	 int k = iy_muVT / IQ;
	 iy_muVT = IA * (iy_muVT - k*IQ) - IR*k;
	 if(iy_muVT < 0) iy_muVT += IM;
	 rnd = am_muVT * ((IM & (ix_muVT ^ iy_muVT)) | (int)1);
	 return rnd;
}

ChemicalPotential::ChemicalPotential()
{
	 this->ownrank = -1;
	 this->muTilde = 0.0;
	 this->interval = (unsigned)((int)-1);
	 this->instances = 0;
	 for(int d=0; d<3; d++)
	 {
	    this->system[d] = 1.0;
	    this->minredco[d] = 0.0;
	    this->maxredco[d] = 1.0;
	    this->control_bottom[d] = 0.0;
	    this->control_top[d] = 1.0;
	 }
	 this->nextid = 10000000;
	 this->globalN = 1;
	 this->globalV = 1.0;
	 this->restrictedControlVolume = false;

	 this->remainingDeletions = list<unsigned>();
	 for(int d=0; d<3; d++) this->remainingInsertions[d] = list<double>();
	 this->remainingInsertionIDs = list<unsigned long>();
	 this->remainingDecisions = list<float>(); 
	 this->reservoir = list<Molecule>();
	 this->id_increment = 1;
}

void ChemicalPotential::setSubdomain(int rank, double x0, double x1, double y0, double y1, double z0, double z1)
{
	 this->ownrank = rank;
	 this->rnd.init(8624);
	 if(!this->restrictedControlVolume)
	 {
	    this->globalV = this->system[0] * this->system[1] * this->system[2];
	    
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

void ChemicalPotential::setSystem(double x, double y, double z, double m)
{   
	 this->system[0] = x; this->system[1] = y; this->system[2] = z;
	 this->molecularMass = m;
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

// note that *C must not contain the halo
// but when the decisions are evaluated, the halo must be taken into account!
//
void ChemicalPotential::prepareTimestep(TMoleculeContainer* cell, DomainDecompBase* comm)
{
	 this->remainingDeletions.clear();
#ifndef NDEBUG
	 for(int d=0; d<3; d++) 
	   assert(this->remainingInsertions[d].empty());
#endif
	 this->remainingDecisions.clear();

	 // get information on the system decomposition
	 //
	 unsigned localN;
	 if( (maxredco[0] < 0.0) || (maxredco[1] < 0.0) ||
	     (maxredco[2] < 0.0) || (minredco[0] > 1.0) ||
	     (minredco[1] > 1.0) || (minredco[2] > 1.0)    ) localN = 0;
	 else if( (minredco[0] < 0.0) || (minredco[1] < 0.0) ||
	          (minredco[2] < 0.0) || (maxredco[0] > 1.0) ||
	    (maxredco[1] > 1.0) || (maxredco[2] > 1.0)    )
			localN = cell->countParticles(
				 this->componentid, this->control_bottom, this->control_top
			);
	 else localN = cell->countParticles(this->componentid);
	 float minrnd = 0.0;
	 float maxrnd = 1.0;
	 this->globalN = comm->Ndistribution(localN, &minrnd, &maxrnd);
#ifndef NDEBUG
	 global_log->debug() << " believes N(" << componentid << ")=" << globalN << ", rho=" << globalN/globalV
	         << ", the decisive density quotient equals " << (float)globalN/globalReducedVolume << "\n";
#endif

	 // construct deletions
	 //
	 float sel, dec;
	 unsigned localIndex;
	 for(unsigned i=0; i < this->instances; i++)
	 {
	    sel = this->rnd.rnd_muVT();
	    dec = this->rnd.rnd_muVT();
#ifndef NDEBUG
			// if(!ownrank) cout << "global index " << sel << " chosen for deletion.\n";
#endif
			if((sel >= minrnd) && (sel < maxrnd))
			{
				 localIndex = (unsigned)floor(localN*(sel-minrnd)/(maxrnd-minrnd));
#ifndef NDEBUG
				 // cout << "rank " << ownrank << " will try to delete index " << localIndex << ".\n";  // \\ //
#endif
				 this->remainingDeletions.push_back(localIndex);
				 this->remainingDecisions.push_back(dec);
			}
	 }

	 int insertions = this->instances;
#ifndef NDEBUG
	 global_log->debug() << "Number of insertions: " << insertions << ".\n";
#endif

	 // construct insertions
	 //
	 float redc[3];
	 double tc[3];
	 for(int i=0; i < insertions; i++)
	 {
	    for(int d=0; d < 3; d++) redc[d] = this->rnd.rnd_muVT();
	    dec = this->rnd.rnd_muVT();
	    if(    (redc[0] >= minredco[0]) && (redc[1] >= minredco[1]) && (redc[2] >= minredco[2]) 
	        && (redc[0] <  maxredco[0]) && (redc[1] <  maxredco[1]) && (redc[2] <  maxredco[2]) )
	    {
	 for(int d=0; d < 3; d++)
	 {
	    tc[d] = control_bottom[d]
	       + redc[d]*(control_top[d] - control_bottom[d]);
	 }
#ifndef NDEBUG
				 cout << "rank " << ownrank << " will try to insert ID "
				      << nextid << " (" << tc[0] << "/" << tc[1]
				      << "/" << tc[2] << ").\n";  // \\ //
#endif
				 for(int d=0; d < 3; d++)
	 {
	    this->remainingInsertions[d].push_back(tc[d]);
	 }
				 this->remainingDecisions.push_back(dec);
				 this->remainingInsertionIDs.push_back(this->nextid);
			}
			this->nextid += id_increment;
	 }
}

bool ChemicalPotential::getDeletion(TMoleculeContainer* cell, double* minco, double* maxco)
{
	 if(this->remainingDeletions.empty()) return false;
	 if(cell->getNumberOfParticles() == 0) return false;
	 
	 unsigned idx = *this->remainingDeletions.begin();
	 this->remainingDeletions.erase(this->remainingDeletions.begin());
	 double tminco[3];
	 double tmaxco[3];
	 if(restrictedControlVolume) for(int d=0; d < 3; d++)
	 {
	    tminco[d] = (minco[d] > control_bottom[d])? minco[d]
	                                              : control_bottom[d];
	    tmaxco[d] = (maxco[d] < control_top[d])? maxco[d]
	                                           : control_top[d];
	 }
	 else for(int d=0; d < 3; d++)
	 {
	    tminco[d] = minco[d];
	    tmaxco[d] = maxco[d];
	 }
	 
	 Molecule* m = cell->begin();
	 int j=0;
	 for(unsigned i=0; (i < idx); i++)
	 {
	    while(( (m->r(0) > tmaxco[0]) || (m->r(1) > tmaxco[1]) ||
	      (m->r(2) > tmaxco[2]) || (m->r(0) < tminco[0]) ||
	      (m->r(1) < tminco[1]) || (m->r(2) < tminco[2]) ||
							(m->componentid() != this->componentid) )
					 && (m != cell->end()))
			{
				 m = cell->next();
				 if(m == cell->end())
	 {
	    if(j == 0) return false;
	    m = cell->begin();
	    j = 0;
	 }
			}
			m = cell->next();
			j++;
			if(m == cell->end())
			{
	 if(j == 0) return false;
	 m = cell->begin();
	 j = 0;
			}
	 }
	 while( (m->r(0) > tmaxco[0]) || (m->r(1) > tmaxco[1]) ||
	        (m->r(2) > tmaxco[2]) || (m->r(0) < tminco[0]) ||
	  (m->r(1) < tminco[1]) || (m->r(2) < tminco[2]) ||
					(m->componentid() != this->componentid) )
	 {
	    m = cell->next();
	    if(m == cell->end())
	    {
	 if(j == 0) return false;
	 m = cell->begin();
			}
	 }
#ifndef NDEBUG
	 global_log->debug() << "ID " << m->id() << " selected for deletion (index " << idx << ")." << std::endl;
#endif
	 assert(m->id() < nextid);
	 return true;
}

// returns 0 if no insertion remains for this subdomain
unsigned long ChemicalPotential::getInsertion(double* ins)
{
	 if(this->remainingInsertionIDs.empty()) return 0;

	 for(int d=0; d<3; d++)
	 {
	    ins[d] = *this->remainingInsertions[d].begin();
	    this->remainingInsertions[d].erase(this->remainingInsertions[d].begin());
	 }
	 unsigned long nextid = *this->remainingInsertionIDs.begin();
	 this->remainingInsertionIDs.erase(this->remainingInsertionIDs.begin());
	 return nextid;
}

bool ChemicalPotential::decideDeletion(double deltaUTilde)
{
	 if(this->remainingDecisions.empty())
	 {
	    global_log->error() << "No decision is possible." << std::endl;
	    exit(1);
	 }
	 float dec = *this->remainingDecisions.begin();
	 this->remainingDecisions.erase(this->remainingDecisions.begin());
	 float acc = ((float)(this->globalN)) * exp(-muTilde-deltaUTilde) / this->globalReducedVolume;
	 bool ans;
	 if(dec < 0.000001) ans = true;
	 else if(dec > 0.999999) ans = false;
	 else ans = (acc > dec);
	 // ans = (acc > 1.0e-05)? (acc > dec): false;
#ifndef NDEBUG
	 // cout << "rank " << ownrank << (ans? " accepted ": " rejected ")
	 //      << "deletion with deltaUtilde = " << deltaUTilde << " (P = "
	 //      << ((acc > 1.0)? 1.0: acc) << ").\n"; // \\ //
#endif
	 if(ans) this->globalN -= (2*ownrank + 1);  // estimate, the precise value is communicated later
	 return ans;
}

bool ChemicalPotential::decideInsertion(double deltaUTilde)
{
	 if(this->remainingDecisions.empty())
	 {
	    global_log->error() << "No decision is possible." << std::endl;
	    exit(1);
	 }
	 float dec = *this->remainingDecisions.begin();
	 this->remainingDecisions.erase(this->remainingDecisions.begin());
	 double acc = this->globalReducedVolume * exp(muTilde - deltaUTilde) / (1.0 + (double)(this->globalN));
	 bool ans;
	 if(dec < 0.000001) ans = true;
	 else if(dec > 0.999999) ans = false;
	 else ans = (acc > dec);
	 // ans = (acc > 1.0e-05)? (acc > dec): false;
#ifndef NDEBUG
	 // cout << "rank " << ownrank << (ans? " accepted ": " rejected ")
	 //      << "insertion with deltaUtilde = " << deltaUTilde
	 //      << " (P = " << ((acc > 1.0)? 1.0: acc) << ").\n"; // \\ //
#endif
	 if(ans) this->globalN += (2*ownrank + 1);  // estimate, the precise value is communicated later
	 return ans;
}

void ChemicalPotential::submitTemperature(double T)
{
	 this->muTilde = this->mu / T;
	 double lambda = 0.39894228 * h / sqrt(molecularMass*T);
	 globalReducedVolume = globalV / (lambda*lambda*lambda);
	 double doOutput = this->rnd.rnd_muVT();
#ifdef NDEBUG
	 if(ownrank) return;
#endif
	 if(doOutput >= 0.01) return;
	 cout << "rank " << ownrank << " sets mu~ <- " << muTilde;
	 cout << ", lambda <- " << lambda;
	 cout << ", and Vred <- " << globalReducedVolume << "\n";
}

void ChemicalPotential::assertSynchronization(DomainDecompBase* comm)
{
	 comm->assertIntIdentity(this->rnd.getIX());
}

void ChemicalPotential::setControlVolume(
	 double x0, double y0, double z0, double x1, double y1, double z1
) {
	 if((x0 >= x1) || (y0 >= y1) || (z0 >= z1))
	 {
	       global_log->error() << "\nInvalid control volume (" << x0 << " / " << y0 
	      << " / " << z0 << ") to (" << x1 << " / " << y1 << " / "
	      << z1 << ")." << std::endl;
			exit(611);
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

Molecule ChemicalPotential::loadMolecule()
{
			assert(!this->reservoir.empty());
			Molecule tmp = this->reservoir.front();
			this->reservoir.pop_front();
			if(this->reservoir.empty())
			{
	 tmp.scale_v(-1.0);
	 this->reservoir.push_back( tmp );
			}
			assert(tmp.componentid() == componentid);
#ifndef NDEBUG
			tmp.check(tmp.id());
#endif
			return tmp;
}
