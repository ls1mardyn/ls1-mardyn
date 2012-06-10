/*
 * GNU GPL version 2
 */

#include "TcTS.h"
#include "utils/Random.h"
#include "molecules/Molecule.h"
#include "particleContainer/LinkedCells.h"
#include "io/ResultWriter.h"
#include "io/XyzWriter.h"
#include "io/CheckpointWriter.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecomposition.h"
#endif

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <math.h>
#include <string.h>

using namespace std;

#define BINS 2048
#define DT 0.002
#define PRECISION 5
#define VARFRACTION 0.125

TcTS::TcTS(int argc, char** argv, Domain* domain, DomainDecompBase** domainDecomposition, Integrator** integrator, ParticleContainer** moleculeContainer, std::list<OutputBase*>* outputPlugins, RDF* irdf, Simulation* isimulation) 
{
   this->_domain = domain;
   this->_domainDecomposition = domainDecomposition;
   this->_integrator = integrator;
   this->_moleculeContainer = moleculeContainer;
   this->_outputPlugins = outputPlugins;
   this->_mrdf = irdf;
   this->_msimulation = isimulation;

   const char* usage = "usage: mkTcTS <prefix> -c <density> [-d <second density>] [-h <height>] [-m <chemical potential>] -N <particles> [-p <pair correlation cutoff>] [-R <cutoff>] [-S] -T <temperature> [-U]\n\n-S\tshift (active by default)\n-U\tunshift\n";
   if((argc < 9) || (argc > 24))
   {
      cout << "There are " << argc
           << " arguments where 9 to 24 should be given.\n\n";
      cout << usage;
      exit(1);
   }

   bool do_shift = true;
   bool in_h = false;
   bool use_mu = false;
   bool gradient = false;

   double cutoff = 2.5;
   double dRDF = 12.0;
   double h = 0.0;
   double mu = 0.0;
   unsigned N = 131072;
   double rho = 0.319;
   double rho2 = 0.319;
   double T = 1.0779;

   char* prefix = argv[2];

   for(int i=3; i < argc; i++)
   {
      if(*argv[i] != '-')
      {
         cout << "Flag expected where '" << argv[i]
              << "' was given.\n\n";
         cout << usage;
         exit(2);
      }
      for(int j=1; argv[i][j]; j++)
      {
         if(argv[i][j] == 'c')
         {
            i++;
            rho = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'd')
         {
            gradient = true;
            i++;
            rho2 = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'h')
         {
            in_h = true;
            i++;
            h = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'm')
         {
            use_mu = true;
            i++;
            mu = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'N')
         {
            i++;
            N = atoi(argv[i]);
            break;
         }
         else if(argv[i][j] == 'p')
         {
            i++;
            dRDF = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'R')
         {
            i++;
            cutoff = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'S') do_shift = true;
         else if(argv[i][j] == 'T')
         {
            i++;
            T = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'U') do_shift = false;
         else
         {
            cout << "Invalid flag '-" << argv[i][j]
                 << "' was detected.\n\n" << usage;
            exit(2); 
         }
      }
   }

   if(in_h && !gradient)
   {
      cout << "The box dimension can only be specified for "
           << "systems with a density gradient.\n\n" << usage;
      exit(5);
   }

   if(!in_h) h = pow((double)N/rho, 1.0/3.0);

   if(gradient) this->domain(h, N, rho, rho2);
   else this->domain(N, rho, dRDF);
   this->write(prefix, cutoff, mu, T, do_shift, use_mu);
}

void TcTS::domain(double t_h, unsigned t_N, double t_rho, double t_rho2)
{
   this->gradient = true;
   this->rho = t_rho;
   this->rho2 = t_rho2;
   this->dRDF = 0.0;

   double V = 2.0 * t_N / (rho + rho2);
   this->box[0] = sqrt(V / t_h);
   this->box[1] = t_h;
   this->box[2] = sqrt(V / t_h);
}

void TcTS::domain(unsigned t_N, double t_rho, double t_RDF)
{
   this->gradient = false;
   this->rho = t_rho;
   this->rho2 = t_rho;
   this->dRDF = t_RDF;

   this->box[0] = pow((double)t_N/rho, 1.0/3.0);
   this->box[1] = this->box[0];
   this->box[2] = this->box[0];
}

void TcTS::write(char* prefix, double cutoff, double mu, double T, bool do_shift, bool use_mu)
{
   unsigned fl_units[3][2];
   double fl_unit[3][2];
   double N_id[2];
   for(int i=0; i < 2; i++)
   {
      N_id[i] = box[0]*(0.5*box[1])*box[2] * ((i == 0)? rho: rho2);
      double N_boxes = N_id[i] / 3.0;
      fl_units[1][i] = (unsigned int) round(
                          pow(
                             (N_boxes * (0.5*box[1]) * (0.5*box[1]))
                                      / (this->box[0] * this->box[2]), 1.0/3.0
                          )
                       );
      if(fl_units[1][i] == 0) fl_units[1][i] = 1;
      double bxbz_id = N_boxes / fl_units[1][i];
      fl_units[0][i] = (unsigned int) round(sqrt(this->box[0] * bxbz_id / this->box[2]));
      if(fl_units[0][i] == 0) fl_units[0][i] = 1;
      fl_units[2][i] = (unsigned int) ceil(bxbz_id / fl_units[0][i]);
      for(int d=0; d < 3; d++) fl_unit[d][i] = ((d == 1)? 0.5: 1.0) * box[d] / (double)fl_units[d][i];
      cout << "Elementary cell " << i << ": " << fl_unit[0][i] << " x " << fl_unit[1][i] << " x " << fl_unit[2][i] << ".\n";
   }

   Random* r = new Random();
   r->init(
      (int)(10000.0*box[0]) - (int)(3162.3*cutoff)
        + (int)(1000.0*T) - (int)(316.23*mu)
        + (int)(100.0*box[1])
   );

   bool fill0[fl_units[0][0]][fl_units[1][0]][fl_units[2][0]][3];
   bool fill1[fl_units[0][1]][fl_units[1][1]][fl_units[2][1]][3];
   unsigned N[2];
   unsigned slots[2];
   for(unsigned l=0; l < 2; l++)
   {
      for(unsigned i=0; i < fl_units[0][l]; i++)
         for(unsigned j=0; j < fl_units[1][l]; j++)
            for(unsigned k=0; k < fl_units[2][l]; k++)
               for(unsigned d=0; d < 3; d++)
               {
                  if(l == 0) fill0[i][j][k][d] = true;
                  else fill1[i][j][k][d] = true;
               }
      slots[l] = 3 * fl_units[0][l] * fl_units[1][l] * fl_units[2][l];
      N[l] = slots[l];
   }
   bool tswap;
   double pswap;
   for(unsigned l=0; l < 2; l++)
   {
      for(unsigned m=0; m < PRECISION; m++)
      {
         tswap = (N[l] < N_id[l]);
         pswap = (N_id[l] - (double)N[l]) / ((tswap? slots[l]: 0) - (double)N[l]);
         // cout << "N = " << N[l] << ", N_id = " << N_id[l] << " => tswap = " << tswap << ", pswap = " << pswap << "\n";
         for(unsigned i=0; i < fl_units[0][l]; i++)
            for(unsigned j=0; j < fl_units[1][l]; j++)
               for(unsigned k=0; k < fl_units[2][l]; k++)
                  for(unsigned d=0; d < 3; d++)
                     if(pswap >= r->rnd())
                     {
                        if(((l == 0) && fill0[i][j][k][d]) || ((l == 1) && fill1[i][j][k][d])) N[l] --;
                        if(l == 0) fill0[i][j][k][d] = tswap;
                        else fill1[i][j][k][d] = tswap;
                        if(tswap) N[l] ++;
                     }
      }
      cout << "Filling " << N[l] << " of 3*"
           << fl_units[0][l] << "*" << fl_units[1][l] << "*" << fl_units[2][l]
           << " = " << slots[l] << " slots (ideally " << N_id[l] << ").\n";
   }

   *this->_integrator = new Leapfrog(DT);
   double ecutoff = (dRDF > cutoff)? 1.031623*dRDF: cutoff;
   this->_msimulation->setcutoffRadius(ecutoff);
   this->_msimulation->setLJCutoff(cutoff);

   this->_domain->setCurrentTime(0.0);
   for( int d = 0; d < 3; d++ ) this->_domain->setGlobalLength(d, box[d]);

   vector<Component>& dcomponents = this->_domain->getComponents();
   dcomponents.resize(1);
   dcomponents[0].setID(0);
   dcomponents[0].addLJcenter(0, 0, 0, 1, 1, 1, cutoff, do_shift);
   dcomponents[0].setI11(0);
   dcomponents[0].setI22(0);
   dcomponents[0].setI33(0);
   this->_domain->setepsilonRF(1.0e+10);

   this->_domain->setGlobalTemperature(T);
   this->_domain->setglobalNumMolecules(N[0]+N[1]);

#ifdef ENABLE_MPI
   *(this->_domainDecomposition) = (DomainDecompBase*) new DomainDecomposition();
#endif
   double bBoxMin[3];
   double bBoxMax[3];
   for (int i = 0; i < 3; i++) {
      bBoxMin[i] = (*_domainDecomposition)->getBoundingBoxMin(i, _domain);
      bBoxMax[i] = (*_domainDecomposition)->getBoundingBoxMax(i, _domain);
   }
   *(this->_moleculeContainer) = new LinkedCells(bBoxMin, bBoxMax, ecutoff, cutoff, 0.5, 1);
   stringstream opstream;
   opstream << prefix << "_1R";
   this->_msimulation->setOutputPrefix(opstream.str().c_str());
   this->_outputPlugins->push_back(new ResultWriter(1500, opstream.str().c_str()));
   this->_outputPlugins->push_back(new XyzWriter(60000, opstream.str().c_str(), 100000000, true));
   this->_outputPlugins->push_back(new CheckpointWriter(60000, opstream.str().c_str(), 100000000, true));
   this->_msimulation->initCanonical(10);
   this->_msimulation->initStatistics(3003003);
   if(gradient)
   {
      this->_domain->setupProfile(1, BINS, 1);
      this->_msimulation->profileSettings(1, 3000000, opstream.str());
      this->_domain->considerComponentInProfile(1);
   }
   else
   {
      this->_mrdf = new RDF(dRDF/(double)BINS, BINS, dcomponents);
      this->_mrdf->setOutputTimestep(3000000);
      this->_mrdf->setOutputPrefix(opstream.str());
   }

   double v = sqrt(3.0 * T);
   double loffset[3][2];
   loffset[0][0] = 0.1; loffset[1][0] = 0.3; loffset[2][0] = 0.1;
   loffset[0][1] = 0.1; loffset[1][1] = 0.8; loffset[2][1] = 0.1;
   double goffset[3][3];
   goffset[0][0] = 0.0; goffset[1][0] = 0.5; goffset[2][0] = 0.5;
   goffset[0][1] = 0.5; goffset[1][1] = 0.0; goffset[2][1] = 0.5;
   goffset[0][2] = 0.5; goffset[1][2] = 0.5; goffset[2][2] = 0.0;

   this->_domain->setglobalRotDOF(0);
   // unsigned maxid = 1;
   unsigned ID = 1;
   for(unsigned l=0; l < 2; l++)
      for(unsigned i=0; i < fl_units[0][l]; i++)
         for(unsigned j=0; j < fl_units[1][l]; j++)
            for(unsigned k=0; k < fl_units[2][l]; k++)
               for(unsigned d=0; d < 3; d++)
               {
                  if(((l == 0) && fill0[i][j][k][d]) || ((l == 1) && fill1[i][j][k][d]))
                  {
                     double q[3];
                     q[0] = i * fl_unit[0][l];
                     q[1] = j * fl_unit[1][l];
                     q[2] = k * fl_unit[2][l];
                     for(int m=0; m < 3; m++)
                     {
                        q[m] += box[m]*loffset[m][l] + fl_unit[m][l]*goffset[m][d];
                        q[m] += VARFRACTION * fl_unit[m][l] * (r->rnd() - 0.5);
                        if(q[m] > box[m]) q[m] -= box[m];
                     }
                     double phi = 6.283185 * r->rnd();
                     double omega = 6.283185 * r->rnd();

                     Molecule m1 = Molecule(ID, 0, q[0], q[1], q[2], v*cos(phi)*cos(omega), v*cos(phi)*sin(omega), v*sin(phi), 1, 0, 0, 0, 0, 0, 0, &dcomponents);
                     (*this->_moleculeContainer)->addParticle(m1);
                     dcomponents[0].incNumMolecules();
                     // domain->setglobalRotDOF(dcomponents[0].getRotationalDegreesOfFreedom() + domain->getglobalRotDOF());
		
                     // if(id > maxid) maxid = id;

                     /*
                     std::list<ChemicalPotential>::iterator cpit;
                     for(cpit = lmu->begin(); cpit != lmu->end(); cpit++) {
                        if( !cpit->hasSample() && (componentid == cpit->getComponentID()) ) {
                           cpit->storeMolecule(m1);
                        }
                     }
                     */

                     ID++;
                  }
               }
}

