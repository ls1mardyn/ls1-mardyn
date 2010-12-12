/*
 * GNU GPL version 2
 */

#include "Domain.h"
#include "Random.h"
#include "Graphit.h"
#include <cmath>

#define DT 0.030620
#define TAU_ZERO 0.30620
#define ZETA 0.0036143
#define TAUPRIME 61.240

Domain::Domain(
      int sp_flow, double sp_bondlength, double sp_rho, int sp_d,
      int sp_fluid, double sp_h, double sp_ETA, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      bool sp_nanotube, double sp_m_per_n, unsigned sp_N, double sp_T,
      double sp_XI, double sp_wo_wall
) {
   this->flow = sp_flow;
   this->bondlength = sp_bondlength;
   this->d = sp_d;
   this->fluid = sp_fluid;
   this->h = sp_h;
   this->ETA = sp_ETA;
   this->SIG_REF = sp_SIG_REF;
   this->EPS_REF = sp_EPS_REF;
   this->REFMASS = sp_REFMASS;
   this->muVT = sp_muVT;
   this->nanotube = sp_nanotube;
   this->pfill = 0.0;
   this->T = sp_T;
   this->XI = sp_XI;
   this->wo_wall = sp_wo_wall;

   if(fluid == FLUID_AR) this->shielding = this->ETA * SIG_AR;
   else if(fluid == FLUID_CH4) this->shielding = ETA * SIG_CH4;
   /*
    * 2CLJQ fluids: shielding = elongation/4 + eta*sigma
    */
   else if(fluid == FLUID_C2H6)
   {
      this->shielding = 0.25*C2H6LONG + this->ETA*SIG_C2H6;
   }
   else if(fluid == FLUID_N2)
   {
      this->shielding = 0.25*N2LONG + this->ETA*SIG_N2;
   }
   else if(fluid == FLUID_CO2)
   {
      this->shielding = 0.25*CO2LONG + this->ETA*SIG_CO2;
   }
   else if(fluid == FLUID_AVE) this->shielding = ETA * SIG_AVE;
   else if(fluid == FLUID_H2O)
   {
      this->shielding = 0.25*H2O_LONG + this->ETA*SIG_OH2O;
   }
   else if(fluid == FLUID_CH3OH)
   {
      this->shielding = 0.25*CH3OH_LONG + this->ETA*SIG_CCH3OH;
   }
   else if(fluid == FLUID_C6H14)
   {
      this->shielding = 0.25*C6H14_LONG + this->ETA*SIG_MC6H14;
   }

   if(nanotube) this->specifyNanotube(sp_rho, sp_m_per_n, sp_N);
   else this->specifyGraphite(sp_rho, sp_N);
}

void Domain::write(
   char* prefix, double a, bool empty, int format, double mu,
   double TAU, double U, bool original, double wo_acceleration,
   double polarity, bool WLJ
) {

   if(this->nanotube) this->writeNanotube(
      prefix, a, empty, format, mu, TAU, U, original,
      wo_acceleration, polarity, WLJ
   );
   else this->writeGraphite(
      prefix, a, empty, format, mu, TAU, U, original,
      wo_acceleration, polarity, WLJ
   );
}

void Domain::writeGraphite(
   char* prefix, double a, bool empty, int format, double mu,
   double TAU, double U, bool original, double wo_acceleration,
   double polarity, bool WLJ
) {
   ofstream xdr, txt;
   stringstream strstrm, txtstrstrm;
   strstrm << prefix << ".xdr";
   xdr.open(strstrm.str().c_str(), ios::trunc);
   txtstrstrm << prefix << "_1R.txt";
   txt.open(txtstrstrm.str().c_str(), ios::trunc);

   Random* r = new Random();
   r->init(
      d + (int)(1000000.0*eff[2]) + (int)(10000.0*U)
        + (int)(100.0*off[1]) - (int)(10000.0*TAU) + (int)(31623.0*a)
        - (int)(316.23*box[1])
   );
   double REFTIME = SIG_REF * sqrt(REFMASS / EPS_REF);
   double VEL_REF = SIG_REF / REFTIME;
   cout << "Velocity unit 1 = " << VEL_REF << " * 1620.34 m/s = "
        << 1620.34 * VEL_REF << " m/s.\n";
   double ACC_REF = VEL_REF / REFTIME;
   double REFCARG = sqrt(EPS_REF * SIG_REF);
   cout << "Charge unit 1 = " << REFCARG << " e.\n";
   double QDR_REF = SIG_REF*SIG_REF * REFCARG;
   double REFOMGA = 1.0 / REFTIME;

   int repl = (this->flow == FLOW_COUETTE)? 2: 1;
   bool fill[fl_units[0]][fl_units[1]][fl_units[2]][repl][3];
   unsigned N1 = 0;
   if(!empty)
   {
      for(int i=0; i < this->fl_units[0]; i++)
         for(int j=0; j < this->fl_units[1]; j++)
            for(int k=0; k < this->fl_units[2]; k++)
               for(int l=0; l < repl; l++)
                  for(int d=0; d < 3; d++)
                  {
                     bool tfill = (pfill >= r->rnd());
                     fill[i][j][k][l][d] = tfill;
                     if(tfill) N1++;
                  }
      cout << "Filling " << N1 << " of " << repl << " x 3*"
           << fl_units[0] << "*" << fl_units[1] << "*" << fl_units[2]
           << " = " << 3*repl*fl_units[0]*fl_units[1]*fl_units[2]
           << " slots (ideally "
           << pfill*3*repl*fl_units[0]*fl_units[1]*fl_units[2]
           << ").\n";
   }
   Graphit gra;
   gra.calculateCoordinatesOfAtoms(
      1, this->box[0], this->box[2], this->bondlength, this->wo_wall
   );
   unsigned Ngraphene = gra.getNumberOfAtoms();
   cout << "Inserting " << repl*d << " x " << Ngraphene
        << " carbon atoms.\n";
   unsigned Ntotal = N1 + repl*d*Ngraphene;

   double LJ_CUTOFF;
   double EL_CUTOFF;
   double FLUIDMASS, EPS_FLUID, SIG_FLUID, FLUIDLONG, QDR_FLUID;
   if(fluid == FLUID_AR)
   {
      FLUIDMASS = ARMASS;
      EPS_FLUID = EPS_AR;
      SIG_FLUID = SIG_AR;
      LJ_CUTOFF = 2.5*SIG_FLUID;
   }
   else if(fluid == FLUID_CH4)
   {
      FLUIDMASS = CH4MASS;
      EPS_FLUID = EPS_CH4;
      SIG_FLUID = SIG_CH4;
      LJ_CUTOFF = 2.5*SIG_FLUID;
   }
   else if(fluid == FLUID_C2H6)
   {
      FLUIDMASS = C2H6MASS;
      EPS_FLUID = EPS_C2H6;
      SIG_FLUID = SIG_C2H6;
      FLUIDLONG = C2H6LONG;
      QDR_FLUID = QDR_C2H6;
      LJ_CUTOFF = 4.0*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_N2)
   {
      FLUIDMASS = N2MASS;
      EPS_FLUID = EPS_N2;
      SIG_FLUID = SIG_N2;
      FLUIDLONG = N2LONG;
      QDR_FLUID = QDR_N2;
      LJ_CUTOFF = 4.0*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_CO2)
   {
      FLUIDMASS = CO2MASS;
      EPS_FLUID = EPS_CO2;
      SIG_FLUID = SIG_CO2;
      FLUIDLONG = CO2LONG;
      QDR_FLUID = QDR_CO2;
      LJ_CUTOFF = 4.0*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_AVE)
   {
      FLUIDMASS = AVEMASS;
      EPS_FLUID = EPS_AVE;
      SIG_FLUID = SIG_AVE;
      LJ_CUTOFF = 2.5*SIG_FLUID;
   }
   else if(fluid == FLUID_H2O)
   {
      FLUIDMASS = OH2OMASS + 2.0*HH2OMASS;
      EPS_FLUID = EPS_OH2O;
      SIG_FLUID = SIG_OH2O;
      FLUIDLONG = H2O_LONG;
      LJ_CUTOFF = 4.0*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_CH3OH)
   {
      FLUIDMASS = CCH3OHMASS + OCH3OHMASS + HCH3OHMASS;
      FLUIDLONG = CH3OH_LONG;
      LJ_CUTOFF = 4.0*SIG_CCH3OH + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_C6H14)
   {
      FLUIDMASS = 2.0*FC6H14MASS + 4.0*MC6H14MASS;
      FLUIDLONG = C6H14_LONG;
      LJ_CUTOFF = 4.0*SIG_MC6H14 + 0.5*FLUIDLONG;
   }
   else
   {
      cout << "Unavailable fluid ID " << fluid << ".\n";
      exit(20);
   }

   if((fluid == FLUID_AR) || (fluid == FLUID_CH4) || (fluid == FLUID_AVE))
   {
      EL_CUTOFF = LJ_CUTOFF;
   }
   else if((fluid == FLUID_C2H6) || (fluid == FLUID_N2) || (fluid == FLUID_CO2))
   {
      EL_CUTOFF = 1.08*LJ_CUTOFF;
   }
   else EL_CUTOFF = 1.1664*LJ_CUTOFF;

   if(format == FORMAT_BRANCH)
   {
      unsigned wallcomp_sng = (polarity == 0.0)? d: NCOMP_POLAR;
      unsigned wallcomp = repl*wallcomp_sng;

      xdr.precision(7);
      xdr << "mardyn " << TIME << " tersoff\n"
          << "# mardyn input file, ls1 project\n"
          << "# generated by the mkcp tool\n";

      xdr << "t\t0.0\ndt\t" << DT/REFTIME << "\n" << "# rho" << "\t"
          << Ntotal * SIG_REF*SIG_REF*SIG_REF
                    / (repl * box[0]*box[1]*box[2])
          << "\n" << "L" << "\t" << box[0]/SIG_REF << "\t"
          << repl*box[1]/SIG_REF << "\t" << box[2]/SIG_REF
          << "\n";
      xdr << "C" << "\t" << wallcomp + 1 << "\n";

      if((fluid == FLUID_AR) || (fluid == FLUID_CH4) || (fluid == FLUID_AVE))
      {
         xdr << "1 0 0 0 0\n"  // LJ, C, Q, D, Tersoff
             << "0.0 0.0 0.0\t"
             << FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF << " "
             << SIG_FLUID/SIG_REF << "\t0.0 0.0 0.0\n";
      }
      else if((fluid == FLUID_C2H6) || (fluid == FLUID_N2) || (fluid == FLUID_CO2))
      {
         xdr << "2 0 1 0 0\n"  // LJ, C, Q, D, Tersoff
             << "0.0 0.0 " << -0.5*FLUIDLONG/SIG_REF << "\t"
             << 0.5*FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF
             << " " << SIG_FLUID/SIG_REF << "\n"
             << "0.0 0.0 " << +0.5*FLUIDLONG/SIG_REF << "\t"
             << 0.5*FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF
             << " " << SIG_FLUID/SIG_REF << "\n"
             << "0.0 0.0 0.0\t0.0 0.0 1.0\t" << QDR_FLUID/QDR_REF
             << "\n0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_H2O)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, Q, D, Tersoff

         xdr << R0_O_H2O/SIG_REF << " " << R1_O_H2O/SIG_REF << " " << R2_O_H2O/SIG_REF << "\t"
             << OH2OMASS/REFMASS << " " << EPS_OH2O/EPS_REF << " " << SIG_OH2O/SIG_REF << "\n";

         xdr << R0_H1H2O/SIG_REF << " " << R1_H1H2O/SIG_REF << " " << R2_H1H2O/SIG_REF << "\t"
             << HH2OMASS/REFMASS << " " << CHG_HH2O/REFCARG << "\n";
         xdr << R0_H2H2O/SIG_REF << " " << R1_H2H2O/SIG_REF << " " << R2_H2H2O/SIG_REF << "\t"
             << HH2OMASS/REFMASS << " " << CHG_HH2O/REFCARG << "\n";
         xdr << R0_E_H2O/SIG_REF << " " << R1_E_H2O/SIG_REF << " " << R2_E_H2O/SIG_REF << "\t"
             << "0.0 " << CHG_EH2O/REFCARG << "\n";

         xdr << "0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_CH3OH)
      {
         xdr << "2 3 0 0 0\n";  // LJ, C, Q, D, Tersoff

         xdr << R0_C_CH3OH/SIG_REF << " " << R1_C_CH3OH/SIG_REF << " " << R2_C_CH3OH/SIG_REF << "\t"
             << CCH3OHMASS/REFMASS << " " << EPS_CCH3OH/EPS_REF << " " << SIG_CCH3OH/SIG_REF << "\n";
         xdr << R0_O_CH3OH/SIG_REF << " " << R1_O_CH3OH/SIG_REF << " " << R2_O_CH3OH/SIG_REF << "\t"
             << OCH3OHMASS/REFMASS << " " << EPS_OCH3OH/EPS_REF << " " << SIG_OCH3OH/SIG_REF << "\n";

         xdr << R0_C_CH3OH/SIG_REF << " " << R1_C_CH3OH/SIG_REF << " " << R2_C_CH3OH/SIG_REF << "\t"
             << "0.0 " << CHG_CCH3OH/REFCARG << "\n";
         xdr << R0_O_CH3OH/SIG_REF << " " << R1_O_CH3OH/SIG_REF << " " << R2_O_CH3OH/SIG_REF << "\t"
             << "0.0 " << CHG_OCH3OH/REFCARG << "\n";
         xdr << R0_H_CH3OH/SIG_REF << " " << R1_H_CH3OH/SIG_REF << " " << R2_H_CH3OH/SIG_REF << "\t"
             << HCH3OHMASS/REFMASS << " " << CHG_HCH3OH/REFCARG << "\n";

         xdr << "0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_C6H14)
      {
         xdr << "6 0 0 0 0\n";  // LJ, C, Q, D, Tersoff

         xdr << R0_F1C6H14/SIG_REF << " " << R1_F1C6H14/SIG_REF << " " << R2_F1C6H14/SIG_REF << "\t"
             << FC6H14MASS/REFMASS << " " << EPS_FC6H14/EPS_REF << " " << SIG_FC6H14/SIG_REF << "\n";
         xdr << R0_M1C6H14/SIG_REF << " " << R1_M1C6H14/SIG_REF << " " << R2_M1C6H14/SIG_REF << "\t"
             << MC6H14MASS/REFMASS << " " << EPS_MC6H14/EPS_REF << " " << SIG_MC6H14/SIG_REF << "\n";
         xdr << R0_M2C6H14/SIG_REF << " " << R1_M2C6H14/SIG_REF << " " << R2_M2C6H14/SIG_REF << "\t"
             << MC6H14MASS/REFMASS << " " << EPS_MC6H14/EPS_REF << " " << SIG_MC6H14/SIG_REF << "\n";
         xdr << R0_M3C6H14/SIG_REF << " " << R1_M3C6H14/SIG_REF << " " << R2_M3C6H14/SIG_REF << "\t"
             << MC6H14MASS/REFMASS << " " << EPS_MC6H14/EPS_REF << " " << SIG_MC6H14/SIG_REF << "\n";
         xdr << R0_M4C6H14/SIG_REF << " " << R1_M4C6H14/SIG_REF << " " << R2_M4C6H14/SIG_REF << "\t"
             << MC6H14MASS/REFMASS << " " << EPS_MC6H14/EPS_REF << " " << SIG_MC6H14/SIG_REF << "\n";
         xdr << R0_F2C6H14/SIG_REF << " " << R1_F2C6H14/SIG_REF << " " << R2_F2C6H14/SIG_REF << "\t"
             << FC6H14MASS/REFMASS << " " << EPS_FC6H14/EPS_REF << " " << SIG_FC6H14/SIG_REF << "\n";

         xdr << "0.0 0.0 0.0\n";
      }
      else
      {
         cout << "Fluid code " << fluid << ": Not yet implemented.\n";
         exit(1000+fluid);
      }

      double crga[repl*NCOMP_POLAR];
      for(unsigned i=0; i < wallcomp; i++) crga[i] = 0.0;
      for(int r=0; r < repl; r++)
      {
         crga[1 + r*NCOMP_POLAR] = polarity/3.0;
         crga[2 + r*NCOMP_POLAR] = -2.0*polarity/3.0;      
         crga[3 + r*NCOMP_POLAR] = -1.0*polarity/3.0;
      }

      for(unsigned i=0; i < wallcomp; i++)
      {
         if(polarity == 0.0) xdr << "1 0 0 0 1\n";  // LJ, C, Q, D, Tersoff
         else xdr << "0 1 0 0 1\n";
         
         if(polarity == 0.0)
            xdr << "0.0 0.0 0.0\t" << 0.5*ATOMIC_MASS_C/REFMASS << " "
                << EPS_FLUID/EPS_REF << " " << SIG_FLUID/SIG_REF << "\n";
         else
            xdr << "0.0 0.0 0.0\t" << 0.5*ATOMIC_MASS_C/REFMASS
                << " " << crga[i]/REFCARG << "\n";

         xdr << "0.0 0.0 0.0\t" << 0.5*ATOMIC_MASS_C/REFMASS << " "
             << TERSOFF_A/EPS_REF << " " << TERSOFF_B/EPS_REF << "\t"
             << (original? TERSOFF_LAMBDA_ORIG
                         : TERSOFF_LAMBDA) * SIG_REF << " "
             << (original? TERSOFF_MU_ORIG: TERSOFF_MU) * SIG_REF
             << " " << (original? TERSOFF_R_ORIG: TERSOFF_R)/SIG_REF
             << " " << (original? TERSOFF_S_ORIG: TERSOFF_S)/SIG_REF
             << "\t" << TERSOFF_C << " " << TERSOFF_D << " "
             << TERSOFF_H << " " << TERSOFF_N << " " << TERSOFF_BETA
             << "\n0.0 0.0 0.0\n";
      }
      for(unsigned j=2; wallcomp + 1 >= j; j++)
         xdr << this->XI << " " << this->ETA << "   ";
      xdr << "\n";
      for(unsigned i=2; wallcomp >= i; i++)
      {
         for(unsigned j = i+1; wallcomp + 1 >= j; j++)
            if(WLJ) xdr << "1.0 1.0   ";
            else xdr << "0.0 1.0   ";
         xdr << "\n";
      }
      xdr << "1.0e+10\n";

      double wallthermostats = (polarity == 0.0)? wallcomp: repl;
      for(unsigned i=2; wallcomp + 1 >= i; i++)
      {
         xdr << "CT\t" << i << " " << ((polarity == 0.0)? i-1: ((wallcomp_sng+1 >= i)? 1: 2)) << "\n";
      }
      for(unsigned i=1; wallthermostats >= i; i++)
      {
         xdr << "ThT\t" << i << " " << T/EPS_REF << "\n";
         if(flow == FLOW_COUETTE) xdr << "U\t" << i << "\n";
      }
      xdr << "CT\t1 " << wallthermostats + 1 << "\n"
          << "ThT\t" << wallthermostats + 1 << " " << T/EPS_REF
          << "\nU\t" << wallthermostats + 1 << "\n";
      if(flow == FLOW_COUETTE)
      {
         if(polarity == 0.0)
         {
            xdr << "S\t2 1\nS\t" << this->d + 1 << " 2\n"
                << "A\t1\t0.0 0.0 " << U/VEL_REF << "\t" << TAU/REFTIME
                << "\t0.0 0.0 " << a/ACC_REF << "\n"
                << "A\t2\t0.0 0.0 " << U/VEL_REF << "\t" << TAU/REFTIME
                << "\t0.0 0.0 " << a/ACC_REF << "\n"
                << "S\t" << this->d + 2 << " 3\n"
                << "S\t" << 2*this->d + 1 << " 4\n"
                << "A\t3\t0.0 0.0 0.0\t" << TAU_ZERO/REFTIME
                << "\t0.0 0.0 " << -1.0*a/ACC_REF << "\n"
                << "A\t4\t0.0 0.0 0.0\t" << TAU_ZERO/REFTIME
                << "\t0.0 0.0 " << -1.0*a/ACC_REF << "\n";
            for(unsigned i=3; this->d >= i; i++)
               xdr << "S\t" << i << " 5\n";
            for(unsigned i=3; this->d >= i; i++)
               xdr << "S\t" << d+i << " 6\n";
            xdr << "A\t" << 5 << "\t0.0 0.0 " << U/VEL_REF << "\t"
                << TAU/REFTIME << "\t0.0 0.0 0.0\n";
            xdr << "A\t" << 6 << "\t0.0 0.0 0.0\t" << TAU_ZERO/REFTIME
                << "\t0.0 0.0 0.0\n";
         }
         else
         {
            for(unsigned i=1; i < wallcomp_sng; i++)
            {
               xdr << "S\t" << i+1 << " 1\nS\t" << wallcomp_sng+i+1 << "\t 2\n";
            }
            xdr << "A\t1\t0.0 0.0 " << U/VEL_REF << "\t" << TAU/REFTIME
                << "\t0.0 0.0 " << a/ACC_REF << "\n"
                << "A\t2\t0.0 0.0 0.0\t" << TAU/REFTIME
                << "\t0.0 0.0 " << -1.0*a/ACC_REF << "\n";
         }
      }
      else if(flow == FLOW_POISEUILLE)
      {
         if(polarity == 0.0)
         {
            xdr << "S\t1 1\nA\t1\t0.0 0.0 " << U/VEL_REF << "\t"
                << TAU/REFTIME << "\t0.0 0.0 " << a/ACC_REF << "\n"
                << "S\t2 2\nS\t" << d+1 << " 3\n"
                << "A\t2\t0.0 0.0 0.0\t" << TAU_ZERO/REFTIME
                << "\t0.0 0.0 0.0\n"
                << "A\t3\t0.0 0.0 0.0\t" << TAU_ZERO/REFTIME
                << "\t0.0 0.0 0.0\n";
            for(unsigned i=3; this->d >= i; i++)
               xdr << "S\t" << i << " " << 4 << "\n";
            xdr << "A\t4\t0.0 0.0 0.0\t" << TAU_ZERO/REFTIME
                << "\t0.0 0.0 0.0\n";
         }
         else
         {
            xdr << "S\t1 1\nA\t1\t0.0 0.0 " << U/VEL_REF << "\t"
                << TAU/REFTIME << "\t0.0 0.0 " << a/ACC_REF << "\n";
            for(unsigned i=1; i < wallcomp; i++)
            {
               xdr << "S\t" << i+1 << " 2\n";
            }
            xdr << "A\t2\t0.0 0.0 0.0\t" << TAU/REFTIME
                << "\t0.0 0.0 0.0\n";
         }
      }
      xdr << "N" << "\t" << Ntotal << "\nM" << "\t" << "ICRVQD\n\n";

      txt.precision(6);
      txt << "mardynconfig\n# \ntimestepLength\t" << DT/REFTIME
          << "\ncutoffRadius\t" << EL_CUTOFF/SIG_REF
          << "\nLJCutoffRadius\t" << LJ_CUTOFF/SIG_REF
          << "\ntersoffCutoffRadius\t"
          << 1.0001*(original? TERSOFF_S_ORIG: TERSOFF_S) / SIG_REF
          << "\nconstantAccelerationTimesteps\t25\n"
          << "zetaFlow\t" << ZETA*sqrt(REFTIME) << "\n"
          << "tauPrimeFlow\t" << TAUPRIME/REFTIME << "\n"
          << "initCanonical\t25001\n";
      if(muVT)
      {
         txt.precision(9);
         txt << "chemicalPotential " << mu/EPS_REF
             << " component 1 control 0.0 0.0 0.0 to "
             << this->box[0]/SIG_REF << " " << 0.5*this->h/SIG_REF
             << " " << this->box[2]/SIG_REF << " conduct "
             << 1 + (int)round(((flow == FLOW_COUETTE)? 0.0001
                                                      : 0.0002) * N1)
             << " tests every 2 steps\n";
         if(flow == FLOW_COUETTE)
         {
            txt << "chemicalPotential " << mu/EPS_REF
                << " component 1 control 0.0 " << this->box[1]/SIG_REF
                << " 0.0 to " << this->box[0]/SIG_REF << " "
                << (this->box[1] + 0.5*this->h)/SIG_REF << " "
                << box[2]/SIG_REF << " conduct "
                << 1 + (int)round(((flow == FLOW_COUETTE)? 0.0001
                                                         : 0.0002) * N1)
                << " tests every 2 steps\n";
         }
	 txt << "planckConstant\t" 
	     << sqrt(6.28319 * T/EPS_REF) << "\n";  // sqrt(2 pi kT)
         txt << "initGrandCanonical\t50001\n";
      }
      txt.precision(5);
      txt << "initStatistics\t100001\nphaseSpaceFile\t" << prefix
          << ".xdr\n# for LinkedCells, the cellsInCutoffRadius has to"
          << " be provided\ndatastructure\tLinkedCells\t1\noutput\t"
          << "ResultWriter\t40\t" << prefix
          << "_1R\noutput\tXyzWriter\t10000\t" << prefix
          << "_1R.buxyz\noutput\tVisittWriter\t10000000\t" << prefix
          << "_1R\nprofile\t6 1024 8\nprofileRecordingTimesteps\t3\n"
          << "profileOutputTimesteps\t50000"
          << "\nprofiledComponent\t1\nprofileOutputPrefix\t" << prefix
          << "_1R\n";
      if(WLJ)
      {
         txt << "wallLJ\ton\n";
         if(flow == FLOW_COUETTE)
         {
             txt << "zOscillator 1024\n";
         }
         else
         {
             txt << "oscillator 1024\n";
         }
      }
      else
      {
         txt << "wallLJ\toff\n";
         if(flow == FLOW_COUETTE)
         {
             txt << "zOscillator 128\n";
         }
         else
         {
             txt << "oscillator 128\n";
         }
      }
      if(flow == FLOW_NONE)
         txt << "nomomentum 16384\n";
      if(wo_acceleration != 0.0)
      {
         txt << "flowControl\t0 0 " << 0.5*wo_acceleration*box[2]/SIG_REF
             << " to " << box[0]/SIG_REF << " " << repl*box[1]/SIG_REF << " "
             << (1.0 - 0.5*wo_acceleration)*box[2]/SIG_REF << "\n";
      }
   }

   double I[3];
   for(int k=0; k < 3; k++) I[k] = 0.0;
   if((fluid == FLUID_C2H6) || (fluid == FLUID_N2) || (fluid == FLUID_CO2))
   {
      I[0] = 0.25 * FLUIDMASS * FLUIDLONG * FLUIDLONG;
      I[1] = I[0];
   }
   else if(fluid == FLUID_H2O)
   {
      I[0] = I_XX_H2O;
      I[1] = I_YY_H2O;
      I[2] = I_ZZ_H2O;
   }
   else if(fluid == FLUID_CH3OH)
   {
      I[0] = I_XX_CH3OH;
      I[1] = I_YY_CH3OH;
      I[2] = I_ZZ_CH3OH;
   }
   else if(fluid == FLUID_C6H14)
   {
      I[0] = I_XX_C6H14;
      I[1] = I_YY_C6H14;
      I[2] = I_ZZ_C6H14;
   }

   unsigned id = 1;
   double tr[3];
   int ii[3];
   if(!empty) for(ii[0]=0; ii[0] < this->fl_units[0]; (ii[0]) ++)
      for(ii[1]=0; ii[1] < this->fl_units[1]; (ii[1]) ++)
         for(ii[2]=0; ii[2] < this->fl_units[2]; (ii[2]) ++)
            for(int j=0; j < repl; j++)
	    {
               for(int d=0; d < 3; d++)
               {
                  if(fill[ ii[0] ][ ii[1] ][ ii[2] ][ j ][ d ])
                  {
                     for(int k=0; k < 3; k++)
                     {
                        tr[k] = off[k] + fl_unit[k] * (
                                   ii[k] + 0.02*r->rnd() + ((k == d)? 0.24: 0.74)
                                );
                     }
                     tr[1] += j*box[1];
		     box[1] *= repl;
                     for(int k=0; k < 3; k++)
                     {
                        if(tr[k] > box[k]) tr[k] -= box[k];
                        else if(tr[k] < 0.0) tr[k] += box[k];
                     }
                     box[1] /= repl;
                     double tv = sqrt(3.0*T / FLUIDMASS);
                     double phi = 6.283185 * r->rnd();
                     double omega = 6.283185 * r->rnd();
                     double w[3];
                     for(int k=0; k < 3; k++)
                        w[k] = (I[k] == 0)? 0.0: ((r->rnd() > 0.5)? 1: -1) * sqrt(2.0*r->rnd()*T / I[k]);
                     xdr << id << " " << 1 << "\t" << tr[0]/SIG_REF
                         << " " << tr[1]/SIG_REF << " " << tr[2]/SIG_REF
                         << "\t" << tv*cos(phi)*cos(omega)/VEL_REF << " "
                         << tv*cos(phi)*sin(omega)/VEL_REF << " "
                         << tv*sin(phi)/VEL_REF << "\t1.0 0.0 0.0 0.0\t"
                         << w[0]/REFOMGA << " " << w[1]/REFOMGA << " "
                         << w[2]/REFOMGA << "\n";
                     id++;
                  }
	          else xdr << "\n";
               }
	    }
   xdr << "\n\n";

   for(int j=0; j < repl; j++)
   {
      double layerU = ((flow == FLOW_COUETTE) && (j == 0))
                    ? U
                    : 0.0;
      gra.calculateVelocities(T, layerU);
      for(unsigned k=2; this->d+1 >= k; k++)
      {
         unsigned layerid = j*d + k;
         double yoffset = 0.5*h + (k-2.0)*Z;
         for(unsigned l=0; l < Ngraphene; l++)
         {
            tr[0] = gra.getX(l);
            tr[1] = gra.getY(l) + yoffset;
            tr[2] = gra.getZ(l);
            unsigned tcid = gra.getComponent(l) + j*NCOMP_POLAR;
            for(int m=0; m < 3; m++)
            {
               tr[m] += (0.004*r->rnd() - 0.002) * bondlength;
            }
            if(j == 1) tr[1] += box[1];
            box[1] *= repl;
            for(int m=0; m < 3; m++)
            {
               if(tr[m] > box[m]) tr[m] -= box[m];
               else if(tr[m] < 0.0) tr[m] += box[m];
            }
            box[1] /= repl;

            xdr << id << " " << ((polarity == 0.0)? layerid: tcid)
                << "\t" << tr[0]/SIG_REF
                << " " << tr[1]/SIG_REF << " " << tr[2]/SIG_REF
                << "\t" << gra.getVelocityX(l)/VEL_REF << " "
                << gra.getVelocityY(l)/VEL_REF << " "
                << gra.getVelocityZ(l)/VEL_REF
                << "\t1.0 0.0 0.0 0.0\t0.0 0.0 0.0\n";

            id++;
         }
         xdr << "\n";
      }
      xdr << "\n";
   }

   xdr.close();
   txt.close();
}

void Domain::writeNanotube(
   char* prefix, double a, bool empty, int format, double mu,
   double TAU, double U, bool original, double wo_acceleration,
   double polarity, bool WLJ
) {
   cout << "Cannot create the nanotube - implementation missing.\n";
   exit(19);
}

void Domain::specifyGraphite(double rho, unsigned N)
{
   /*
    * y direction
    */
   this->box[1] = this->h + ((double)(this->d)-1.0)*Z;
   if(1.1 * this->shielding > 0.5 * this->h)
   {
      cout << "Warning: shielding = " << shielding
           << " versus h = " << h << ", ";
      this->shielding = 0.5*this->h - 0.1 * this->shielding;
      if(this->shielding < 0.0) this->shielding = 0.0;
      cout << "corrected to shielding = " << shielding << ".\n";
   }
   this->eff[1] = this->h - 2.0*this->shielding;
   this->off[1] = 0.5*this->h + ((double)(this->d)-1.0)*Z
                              + this->shielding;

   double V_id = (double)(N/rho) / (wo_wall + (1.0 - wo_wall)*eff[1]/box[1]);
   if(this->flow == FLOW_COUETTE) V_id *= 0.5;
   cout << "Carbon-carbon bond length: " << bondlength
        << " * 0.05291772 nm.\n";
   cout << "Total volume should approach " << V_id
        << " * 1.4818e-04 nm^3.\n";

   /*
    * x and z direction
    */
   double A = V_id / this->box[1];
   double tX = this->bondlength * 3.0;
   double tZ = this->bondlength * 12.124356 / (1.0 - wo_wall);  // effektiv 14 * sin(pi/3)
   int zeta = round(sqrt(2.0*A) / tZ);  // bewirkt, dass etwa <z> = 2<x> wird
   if(zeta == 0) zeta = 1;
   int xi = round(A / (tX*tZ*zeta));
   if(xi == 0)
   {
      xi = 1;
      zeta = round(A / (tX*tZ*xi));
      if(zeta == 0)
      {
         cout << "Warning: The generated box will be larger than specified, due to technical reasons.\n\n";
         zeta = 1;
      }
   }
   this->box[0] = xi*tX;
   this->eff[0] = this->box[0];
   this->box[2] = zeta*tZ;
   this->eff[2] = this->box[2];

   /*
    * fluid unit box dimensions
    */
   double V_eff = this->eff[0] * this->eff[1] * this->eff[2];
   cout << "Symmetry volume " << box[0]*box[1]*box[2]
        << " * 1.4818e-04 nm^3, effectively " << V_eff
        << " * 1.4818e-04 nm^3.\n";
   double V_ext = 0.0;
   if(wo_wall > 2.0*shielding)
   {
      V_ext = (wo_wall - 2.0*shielding) * (this->box[1] - this->eff[1]);
      cout << "Additionally available: " << V_ext << " * 1.4818e-04 nm^3.\n";
   }
   double V_total = V_eff + V_ext;
   double N_id = V_total*rho;
   double N_boxes = N_id / 3.0;
   fl_units[1] = round(
                    pow(
                       (N_boxes * this->eff[1] * this->eff[1])
                                / (this->eff[0] * this->eff[2]), 1.0/3.0
                    )
                 );
   if(fl_units[1] == 0) fl_units[1] = 1;
   double bxbz_id = N_boxes / fl_units[1];
   fl_units[0] = round(sqrt(this->eff[0] * bxbz_id / this->eff[2]));
   if(fl_units[0] == 0) this->fl_units[0] = 1;
   fl_units[2] = ceil(bxbz_id / fl_units[0]);
   cout << "Elementary cell: " << this->eff[0]/fl_units[0] << " a0 x " << this->eff[1]/fl_units[1] << " a0 x " << this->eff[2]/fl_units[2] << " a0.\n";

   for(int i=0; i < 3; i++)
      this->fl_unit[i] = this->eff[i] / (double)fl_units[i];
   this->pfill = N_boxes / ((double)fl_units[0]*fl_units[1]*fl_units[2]);
}

void Domain::specifyNanotube(double rho, double m_per_n, unsigned N)
{
   cout << "Nanotubes are not yet implemented.\n";
   exit(16);
}
