#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>


#define TIME 20111116

#define FLUID_NIL -1
#define FLUID_AR 0
#define FLUID_CH4 1
#define FLUID_C2H6 2
#define FLUID_N2 3
#define FLUID_CO2 4
#define FLUID_AVE 5
#define FLUID_H2O 6
#define FLUID_CH3OH 7
#define FLUID_C6H14 8

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
#define SIG_AVE 5.89406
#define EPS_AVE 0.00189803
#define AVEMASS 0.018015

#define EPS_OH2O 6.58951e-04
#define SIG_OH2O 5.89277
#define OH2OMASS 0.0159994
#define HH2OMASS 0.0010079
#define CHG_HH2O +0.419548
#define CHG_EH2O -0.839096
#define R0_O_H2O +0.0
#define R1_O_H2O -0.14948
#define R2_O_H2O +0.0
#define R0_H1H2O -1.72579
#define R1_H1H2O +1.18644
#define R2_H1H2O +0.0
#define R0_H2H2O +1.72579
#define R1_H2H2O +1.18644
#define R2_H2H2O +0.0
#define R0_E_H2O +0.0
#define R1_E_H2O +0.23765
#define R2_E_H2O +0.0
#define I_XX_H2O 0.0031950
#define I_YY_H2O 0.0060038
#define I_ZZ_H2O 0.0091988
#define H2O_LONG 0.29896

#define EPS_CCH3OH 3.81893e-04
#define SIG_CCH3OH 7.09460
#define CCH3OHMASS 0.0150348
#define CHG_CCH3OH +0.24746
#define EPS_OCH3OH 2.78297e-04
#define SIG_OCH3OH 5.72587
#define OCH3OHMASS 0.0159994
#define CHG_OCH3OH -0.67874
#define HCH3OHMASS 0.0010079
#define CHG_HCH3OH +0.43128
#define R0_C_CH3OH -1.44677
#define R1_C_CH3OH -0.05327
#define R2_C_CH3OH +0.0
#define R0_O_CH3OH +1.24534
#define R1_O_CH3OH -0.05327
#define R2_O_CH3OH +0.0
#define R0_H_CH3OH +1.81292
#define R1_H_CH3OH +1.64012
#define R2_H_CH3OH +0.0
#define I_XX_CH3OH 0.0027993
#define I_YY_CH3OH 0.0595954
#define I_ZZ_CH3OH 0.0623947
#define CH3OH_LONG 2.6942

#define EPS_FC6H14 3.10e-04
#define SIG_FC6H14 7.086
#define FC6H14MASS 0.01503
#define EPS_MC6H14 1.46e-04
#define SIG_MC6H14 7.464
#define MC6H14MASS 0.01403
#define R0_F1C6H14 -6.139
#define R1_F1C6H14 +0.418
#define R2_F1C6H14 +0.0
#define R0_M1C6H14 -3.606
#define R1_M1C6H14 -1.015
#define R2_M1C6H14 +0.0
#define R0_M2C6H14 -1.266
#define R1_M2C6H14 +0.716
#define R2_M2C6H14 +0.0
#define R0_M3C6H14 +1.266
#define R1_M3C6H14 -0.716
#define R2_M3C6H14 +0.0
#define R0_M4C6H14 +3.606
#define R1_M4C6H14 +1.015
#define R2_M4C6H14 +0.0
#define R0_F2C6H14 +6.139
#define R1_F2C6H14 -0.418
#define R2_F2C6H14 +0.0
#define I_XX_C6H14 0.04855
#define I_YY_C6H14 1.54288
#define I_ZZ_C6H14 1.59144
#define C6H14_LONG 12.31

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
      int sp_fluid, int sp_fluid2, double sp_h, double sp_ETA, double sp_ETA2, double sp_ETAF, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      bool sp_nanotube, double sp_m_per_n, unsigned sp_N, double sp_T,
      double sp_XI, double sp_XI2, double sp_XIF, double sp_wo_wall
   );
   void write(
      char* prefix, double a, bool empty, int format, double mu,
      double TAU, double U, bool original, double wo_acceleration,
      double polarity, bool WLJ, bool symmetric, bool widom, double x
   );

 private:
   void specifyGraphite(double rho, unsigned N);
   void specifyNanotube(double rho, double m_per_n, unsigned N);
   void writeGraphite(
      char* prefix, double a, bool empty, int format, double mu, 
      double TAU, double U, bool original, double wo_acceleration,
      double polarity, bool WLJ, bool symmetric, bool widom, double x
   );
   void writeNanotube(
      char* prefix, double a, bool empty, int format, double mu,
      double TAU, double U, bool original, double wo_acceleration,
      double polarity, bool WLJ, bool symmetric, bool widom, double x
   );

   bool muVT, nanotube;
   int flow, fluid, fluid2;
   unsigned d, m, n;
   double bondlength, h, ETA, ETA2, ETAF, SIG_REF, EPS_REF, REFMASS, T, XI, XI2, XIF, wo_wall;

   double box[3],  // box dimensions (Couette == two boxes)
          eff[3],  // effective space for the fluid
          off[3];  // offset coordinates for the fluid
   double shielding;  // shielding of the wall due to ETA*SIGMA

   double fl_unit[3];
   unsigned fl_units[3];
   double pfill;  // probability of placing a molecule in a fluid unit

   bool do_fill_ext;
   double ext[3], off_ext[3];
   double fl_unit_ext[3];
   unsigned fl_units_ext[3];
   double pfill_ext;  // probability of placing a molecule in a fluid unit
};

