/*
 * GNU GPL version 2
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;

#define FORMAT_BUCHHOLZ 0
#define FORMAT_BRANCH 1
#define FORMAT_BERNREUTHER 2

class Domain
{
 public:
   Domain(
      double t_h, unsigned t_N, double t_rho, double t_rho2
   );
   Domain(
      unsigned t_N, double t_rho, double t_RDF
   );
   void write(
      char* prefix, double cutoff, double mu, double T, bool do_shift, bool use_mu, int format
   );

 private:
   double box[3];

   bool gradient;
   double rho;
   double rho2;
   double RDF;
};

