/*
 * GNU GPL version 2
 */

#ifndef TCTS_H
#define TCTS_H

#include "Simulation.h"
#include "Domain.h"
#include "RDF.h"
#include "integrators/Leapfrog.h"
#include "parallel/DomainDecompBase.h"
#include "io/OutputBase.h"
#include "utils/OptionParser.h"
using optparse::Values;

using namespace std;

class TcTS
{
 public:
   TcTS(Values &options, Domain* domain, DomainDecompBase** domainDecomposition, Integrator** integrator, ParticleContainer** moleculeContainer, std::list<OutputBase*>* outputPlugins, RDF* irdf, Simulation* isimulation);
   void write(
      char* prefix, double cutoff, double mu, double T, bool do_shift, bool use_mu
   );

 private:
   void domain(
      double t_h, unsigned t_N, double t_rho, double t_rho2
   );
   void domain(
      unsigned t_N, double t_rho, double t_RDF
   );

   double box[3];

   bool gradient;
   double rho;
   double rho2;
   double dRDF;

   Domain* _domain;
   DomainDecompBase** _domainDecomposition;
   Integrator** _integrator;
   ParticleContainer** _moleculeContainer;
   std::list<OutputBase*>* _outputPlugins;
   RDF* _mrdf;
   Simulation* _msimulation;
};

#endif

