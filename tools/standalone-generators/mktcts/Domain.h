#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>


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
      char* prefix, double cutoff, double mu, double T, bool do_shift, bool use_mu, bool compute_autocorr, int format
   );
   void hato(double in1, double in2, double inc1, double inc2) { use_hato = true; p1 = in1; p2 = in2; mu_low = inc1; mu_high = inc2; }

 private:
   double box[3];

   bool gradient;
   double rho;
   double rho2;
   double RDF;
   
   bool use_hato;
   double p1, p2;
   double mu_low, mu_high;
};
