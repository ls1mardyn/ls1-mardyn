/*
 * GNU GPL version 2
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;

#define TIME 20110209

#define FLUID_CH4 0
#define FLUID_AR 1
#define FLUID_C2H6 2
#define FLUID_N2 3
#define FLUID_CO2 4
#define FLUID_EOX 5

#define EPS_AR 3.69853e-04
#define SIG_AR 6.41600
#define ARMASS 0.039948
#define CUT_AR 38.496
#define EPS_CH4 4.70431e-04
#define SIG_CH4 7.04599
#define CH4MASS 0.016042
#define CUT_CH4 42.276

#define EPS_C2H6 4.33822e-04
#define SIG_C2H6 6.59439
#define C2H6MASS 0.03007
#define C2H6LONG 4.49037
#define QDR_C2H6 -0.61537
#define CUT_C2H6 40.689
#define EPS_N2 1.10512e-04
#define SIG_N2 6.27597
#define N2MASS 0.0280134
#define N2LONG 1.97741
#define QDR_N2 -1.07038
#define CUT_N2 38.150
#define EPS_CO2 4.21883e-04
#define SIG_CO2 5.64027
#define CO2MASS 0.0440095
#define CO2LONG 4.56860
#define QDR_CO2 -2.82059
#define CUT_CO2 34.984

#define EPS_CEOX 2.68353e-04
#define SIG_CEOX 6.66431
#define CEOXMASS 0.014027
#define EPS_OEOX 1.96742e-04
#define SIG_OEOX 5.84473
#define OEOXMASS 0.015999
#define DIPOLEOX -0.96744
#define R0_C1EOX +1.47399
#define R1_C1EOX +0.0
#define R2_C1EOX -0.83729
#define R0_C2EOX -1.47399
#define R1_C2EOX +0.0
#define R2_C2EOX -0.83729
#define R0_O_EOX +0.0
#define R1_O_EOX +0.0
#define R2_O_EOX +1.46818
#define R0DIPEOX +0.0
#define R1DIPEOX +0.0
#define R2DIPEOX +0.07792
#define CUTLJEOX 40.833
#define I_XX_EOX  54.154
#define I_YY_EOX 115.11
#define I_ZZ_EOX  60.951

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
      int sp_fluid, double* sp_box, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      unsigned sp_N, double sp_T
   );
   void write(char* prefix, int format, double mu);

 private:
   bool muVT;
   int fluid;
   unsigned N;
   double SIG_REF, EPS_REF, REFMASS, T;

   double box[3];  // offset coordinates for the fluid
};

