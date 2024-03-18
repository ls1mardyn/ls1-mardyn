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
   Domain(double tRin, double tRout, double trhoin, double trhoout) {
      R_i = tRin; R_o = tRout; rho_i = trhoin; rho_o = trhoout;
   };

   void write(char* prefix, double cutoff, double T, bool do_shift, int format);

 private:
   double box[3];

   double R_i, R_o, rho_i, rho_o;
};

