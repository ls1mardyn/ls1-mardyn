#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>


#define TIME 20130802

#define FLUID_NIL -1
#define FLUID_CH4 0
#define FLUID_AR 1
#define FLUID_C2H6 2
#define FLUID_N2 3
#define FLUID_CO2 4
#define FLUID_EOX 5
#define FLUID_JES 6
#define FLUID_MER 7
#define FLUID_TOL 8
#define FLUID_VEG 9

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
#define I_XX_EOX 0.054154
#define I_YY_EOX 0.11511
#define I_ZZ_EOX 0.60951

#define EPS_OJES 6.58951e-04
#define SIG_OJES 5.89277
#define OJESMASS 0.0159994
#define HJESMASS 0.0010079
#define CHG_HJES +0.419548
#define CHG_EJES -0.839096
#define R0_O_JES +0.0
#define R1_O_JES -0.14948
#define R2_O_JES +0.0
#define R0_H1JES -1.72579
#define R1_H1JES +1.18644
#define R2_H1JES +0.0
#define R0_H2JES +1.72579
#define R1_H2JES +1.18644
#define R2_H2JES +0.0
#define R0_E_JES +0.0
#define R1_E_JES +0.23765
#define R2_E_JES +0.0
#define I_XX_JES 0.0031950
#define I_YY_JES 0.0060038
#define I_ZZ_JES 0.0091988
#define JES_LONG 1.0

#define EPS_OMER 3.18243e-04
#define SIG_OMER 5.62288
#define OMERMASS 0.0159994
#define EPS_CMER 3.9181e-05
#define SIG_CMER 5.31712
#define QDR_CMER -3.02883
#define CMERMASS 0.012011
#define MER_LONG 4.8638

#define R1_CH3_TOL +0.0
#define R2_CH3_TOL -5.1999
#define R1_CTR_TOL +0.0
#define R2_CTR_TOL -1.8130
#define R1_CHA_TOL +2.9707
#define R2_CHA_TOL -0.8715
#define R1_CHB_TOL -2.9707
#define R2_CHB_TOL -0.8715
#define R1_CHC_TOL +2.9767
#define R2_CHC_TOL +2.5625
#define R1_CHD_TOL -2.9767
#define R2_CHD_TOL +2.5625
#define R1_CHE_TOL +0.0
#define R2_CHE_TOL +4.2958
#define CH3TOLMASS 0.015035
#define CTRTOLMASS 0.012011
#define CH_TOLMASS 0.013019
#define EPS_CH3TOL 3.9107e-04
#define EPS_CTRTOL 3.4645e-05
#define EPS_CH_TOL 3.18328e-04
#define SIG_CH3TOL 6.77656
#define SIG_CTRTOL 5.2799
#define SIG_CH_TOL 6.1907
#define DIP_CTRTOL 0.173144
#define QDR_CH_TOL -1.2548
#define I_XX_TOL 1.33752
#define I_YY_TOL 0.87701
#define I_ZZ_TOL 0.46050
#define TOL_LONG 9.4957

#define EPS_OVEG 2.9515e-04
#define SIG_OVEG 5.9695
#define OVEGMASS 0.0159994
#define HVEGMASS 0.0010079
#define CHG_HVEG +0.5564
#define CHG_EVEG -1.1128
#define R0_O_VEG +0.0
#define R1_O_VEG -0.12389
#define R2_O_VEG +0.0
#define R0_H1VEG -1.43052
#define R1_H1VEG +0.98330
#define R2_H1VEG +0.0
#define R0_H2VEG +1.43052
#define R1_H2VEG +0.98330
#define R2_H2VEG +0.0
#define R0_E_VEG +0.0
#define R1_E_VEG +0.16826
#define R2_E_VEG +0.0
#define I_XX_VEG 0.0021946
#define I_YY_VEG 0.0041251
#define I_ZZ_VEG 0.0063197
#define VEG_LONG 1.0

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
      int sp_fluid, int sp_fluid2, double* sp_box, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      unsigned sp_N, double sp_T, double sp_ETAF, double sp_XIF
   );
   void write(char* prefix, int format, double mu, double x);

 private:
   bool muVT;
   int fluid, fluid2;
   unsigned N;
   double SIG_REF, EPS_REF, REFMASS, T, ETAF, XIF;

   double box[3];  // offset coordinates for the fluid
};

