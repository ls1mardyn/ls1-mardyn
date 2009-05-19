/*
 * GNU GPL version 2
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;

#define TIME 20090511

#define FLUID_AR 0
#define FLUID_CH4 1
#define FLUID_C2H6 2
#define FLUID_N2 3
#define FLUID_CO2 4

#define EPS_AR 4.36704e-04
#define SIG_AR 6.40920
#define ARMASS 0.039948
#define EPS_CH4 5.54383e-04
#define SIG_CH4 7.03753
#define CH4MASS 0.016042
#define EPS_C2H6 4.33822e-04
#define SIG_C2H6 6.59439
#define C2H6MASS 0.03007
#define C2H6LONG 4.49037
#define QDR_C2H6 -0.61537
#define EPS_N2 1.10512e-04
#define SIG_N2 6.27597
#define N2MASS 0.0280134
#define N2LONG 1.97741
#define QDR_N2 -1.07038
#define EPS_CO2 4.21883e-04
#define SIG_CO2 5.64027
#define CO2MASS 0.0440095
#define CO2LONG 4.56860
#define QDR_CO2 -2.82059

#define FORMAT_BUCHHOLZ 0
#define FORMAT_BRANCH 1
#define FORMAT_BERNREUTHER 2

#define FLOW_NONE 0
#define FLOW_COUETTE 1
#define FLOW_POISEUILLE 2

class Domain
{
 public:
   Domain(
      int sp_flow, double sp_bondlength, double sp_rho, int sp_d,
      int sp_fluid, double sp_h, double sp_ETA, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      bool sp_nanotube, double sp_m_per_n, unsigned sp_N, double sp_T,
      double sp_XI
   );
   void write(
      char* prefix, double a, bool empty, int format, double mu,
      double TAU, double U, bool original
   );

 private:
   void specifyGraphite(double rho, unsigned N);
   void specifyNanotube(double rho, double m_per_n, unsigned N);
   void writeGraphite(
      char* prefix, double a, bool empty, int format, double mu, 
      double TAU, double U, bool original
   );
   void writeNanotube(
      char* prefix, double a, bool empty, int format, double mu,
      double TAU, double U, bool original
   );

   bool muVT, nanotube;
   int flow, fluid;
   unsigned d, m, n;
   double bondlength, h, ETA, SIG_REF, EPS_REF, REFMASS, T, XI;

   double box[3], r_tub,  // box dimensions (Couette == two boxes)
          eff[3], r_eff,  // effective space for the fluid
          off[3], r_off;  // offset coordinates for the fluid
   double shielding;  // shielding of the wall due to ETA*SIGMA

   int fl_unit[3], fl_units[3], fl_unit_r, fl_units_r,
       fl_unit_phi, fl_units_phi;
   double pfill;  // probability of placing a molecule in a fluid unit
};
