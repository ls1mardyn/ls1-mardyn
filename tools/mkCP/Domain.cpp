/*
 * GNU GPL version 2
 */

#include "Domain.h"
#include "Random.h"
#include "Graphit.h"
#include <cmath>

#define DT 0.030620
#define TAU_ZERO 15.3100

Domain::Domain(
      int sp_flow, double sp_bondlength, double sp_rho, int sp_d,
      int sp_fluid, double sp_h, double sp_ETA, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      bool sp_nanotube, double sp_m_per_n, unsigned sp_N, double sp_T,
      double sp_XI
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

   if(fluid == FLUID_AR) this->shielding = this->ETA * SIG_AR;
   else if(fluid == FLUID_CH4) this->shielding = ETA * SIG_CH4;
   /*
    * 2CLJQ fluids: shielding = elongation/4 + eta*sigma
    */
   else if(fluid == FLUID_C2H6)
   {
      this->shielding = 0.25*CO2LONG + this->ETA*SIG_CO2;
   }
   else if(fluid == FLUID_N2)
   {
      this->shielding = 0.25*N2LONG + this->ETA*SIG_N2;
   }
   else if(fluid == FLUID_CO2)
   {
      this->shielding = 0.25*CO2LONG + this->ETA*SIG_CO2;
   }

   if(nanotube) this->specifyNanotube(sp_rho, sp_m_per_n, sp_N);
   else this->specifyGraphite(sp_rho, sp_N);
}

void Domain::write(
   char* prefix, double a, bool empty, int format, double mu,
   double TAU, double U, bool original
) {

   if(this->nanotube) this->writeNanotube(
      prefix, a, empty, format, mu, TAU, U, original
   );
   else this->writeGraphite(
      prefix, a, empty, format, mu, TAU, U, original
   );
}

void Domain::writeGraphite(
   char* prefix, double a, bool empty, int format, double mu,
   double TAU, double U, bool original
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
   bool fill[fl_units[0]][fl_units[1]][fl_units[2]][repl];
   unsigned N1 = 0;
   if(!empty)
   {
      for(int i=0; i < this->fl_units[0]; i++)
         for(int j=0; j < this->fl_units[1]; j++)
            for(int k=0; k < this->fl_units[2]; k++)
               for(int l=0; l < repl; l++)
               {
                  bool tfill = (pfill >= r->rnd());
                  fill[i][j][k][l] = tfill;
                  if(tfill) N1++;
               }
      cout << "Filling " << N1 << " of " << repl << " x "
           << fl_units[0] << "*" << fl_units[1] << "*" << fl_units[2]
           << " = " << repl*fl_units[0]*fl_units[1]*fl_units[2]
           << " slots (ideally "
           << pfill*repl*fl_units[0]*fl_units[1]*fl_units[2]
           << ").\n";
   }
   Graphit gra;
   gra.calculateCoordinatesOfAtoms(
      1, this->box[0], this->box[2], this->bondlength
   );
   unsigned Ngraphene = gra.getNumberOfAtoms();
   cout << "Inserting " << repl*d << " x " << Ngraphene
        << " carbon atoms.\n";
   unsigned Ntotal = N1 + repl*d*Ngraphene;

   double LJ_CUTOFF;
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
      LJ_CUTOFF = 5.5*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_N2)
   {
      FLUIDMASS = N2MASS;
      EPS_FLUID = EPS_N2;
      SIG_FLUID = SIG_N2;
      FLUIDLONG = N2LONG;
      QDR_FLUID = QDR_N2;
      LJ_CUTOFF = 5.5*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_CO2)
   {
      FLUIDMASS = CO2MASS;
      EPS_FLUID = EPS_CO2;
      SIG_FLUID = SIG_CO2;
      FLUIDLONG = CO2LONG;
      QDR_FLUID = QDR_CO2;
      LJ_CUTOFF = 5.5*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else
   {
      cout << "Unavailable fluid ID " << fluid << ".\n";
      exit(20);
   }

   if(format == FORMAT_BRANCH)
   {
      xdr.precision(6);
      xdr << "mardyn " << TIME << " tersoff\n"
          << "# mardyn input file, ls1 project\n"
          << "# generated by the mkcp tool\n";

      xdr << "t\t0.0\ndt\t" << DT/REFTIME << "\n" << "# rho" << "\t"
          << Ntotal * SIG_REF*SIG_REF*SIG_REF
                    / (repl * box[0]*box[1]*box[2])
          << "\n" << "L" << "\t" << box[0]/SIG_REF << "\t"
          << repl*box[1]/SIG_REF << "\t" << box[2]/SIG_REF
          << "\n" << "C" << "\t" << repl*d + 1 << "\n";

      if((fluid == FLUID_AR) || (fluid == FLUID_CH4))
      {
         xdr << "1 0 0 0 0\n"  // LJ, C, Q, D, Tersoff
             << "0.0 0.0 0.0\t"
             << FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF << " "
             << SIG_FLUID/SIG_REF << "\t0.0 0.0 0.0\n";
      }
      else
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
      for(unsigned i=0; i < repl * this->d; i++)
      {
         xdr << "1 0 0 0 1\n"  // LJ, C, Q, D, Tersoff
             << "0.0 0.0 0.0\t" << 0.5*ATOMIC_MASS_C/REFMASS << " "
             << EPS_FLUID/EPS_REF << " " << SIG_FLUID/SIG_REF << "\n"
             << "0.0 0.0 0.0\t" << 0.5*ATOMIC_MASS_C/REFMASS << " "
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
      for(unsigned j=2; repl*this->d + 1 >= j; j++)
         xdr << this->XI << " " << this->ETA << "   ";
      xdr << "\n";
      for(unsigned i=2; repl*this->d >= i; i++)
      {
         for(unsigned j = i+1; repl*d + 1 >= j; j++)
            xdr << "1.0 1.0   ";
         xdr << "\n";
      }
      xdr << "1.0e+10\n";

      for(unsigned i=2; repl*d + 1 >= i; i++)
         xdr << "CT\t" << i << " " << i-1
             << "\nThT\t" << i-1 << " " << T/EPS_REF << "\n";
      if(flow == FLOW_COUETTE)
         for(unsigned i=1; d >= i; i++) xdr << "U\t" << i << "\n";
      xdr << "CT\t1 " << repl*d + 1 << "\n"
          << "ThT\t" << repl*d + 1 << " " << T/EPS_REF
          << "\nU\t" << repl*d + 1 << "\n";
      if(flow == FLOW_COUETTE)
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
      else if(flow == FLOW_POISEUILLE)
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
      xdr << "N" << "\t" << Ntotal << "\nM" << "\t" << "ICRVQD\n";

      txt.precision(6);
      txt << "mardynconfig\n# \ntimestepLength\t" << DT/REFTIME
          << "\ncutoffRadius\t" << LJ_CUTOFF/SIG_REF
          << "\ntersoffCutoffRadius\t"
          << 1.0001*(original? TERSOFF_S_ORIG: TERSOFF_S) / SIG_REF
          << "\nconstantAccelerationTimesteps\t25\n"
          << "initCanonical\t10000\n";
      if(muVT)
      {
         txt.precision(9);
         txt << "chemicalPotential " << mu/EPS_REF
             << " component 1 control 0.0 0.0 0.0 to "
             << this->box[0]/SIG_REF << " " << 0.5*this->h/SIG_REF
             << " " << this->box[2]/SIG_REF << " conduct "
             << 1 + (int)round(((flow == FLOW_COUETTE)? 0.0005
                                                      : 0.0010) * N1)
             << " tests every 8 steps\n";
         if(flow == FLOW_COUETTE)
         {
            txt << "chemicalPotential " << mu/EPS_REF
                << " component 1 control 0.0 " << this->box[1]/SIG_REF
                << " 0.0 to " << this->box[0]/SIG_REF << " "
                << (this->box[1] + 0.5*this->h)/SIG_REF << " "
                << box[2]/SIG_REF << " conduct "
                << 1 + (int)round(((flow == FLOW_COUETTE)? 0.0005
                                                       : 0.0010) * N1)
                << " tests every 8 steps\n";
         }
	 txt << "planckConstant\t" 
	     << sqrt(6.28319 * T/EPS_REF) << "\n";  // sqrt(2 pi kT)
         txt << "initGrandCanonical\t30000\n";
      }
      txt.precision(5);
      txt << "initStatistics\t60000\nphaseSpaceFile\t" << prefix
          << ".xdr\n# for LinkedCells, the cellsInCutoffRadius has to"
          << " be provided\ndatastructure\tLinkedCells\t1\noutput\t"
          << "ResultWriter\t40\t" << prefix
          << "_1R\noutput\tXyzWriter\t10000\t" << prefix
          << "_1R.buxyz\noutput\tVisittWriter\t10000000\t" << prefix
          << "_1R\nprofile\t23 128 32\nprofileRecordingTimesteps\t2\n"
          << "profileOutputTimesteps\t75000"
          << "\nprofiledComponent\t1\nprofileOutputPrefix\t" << prefix
          << "_1R\nzOscillator 512\n";
   }

   double I_xx_yy;
   if((fluid == FLUID_AR) || (fluid == FLUID_CH4)) I_xx_yy = 0.0;
   else I_xx_yy = 0.25 * FLUIDMASS * FLUIDLONG * FLUIDLONG;

   unsigned id = 1;
   double tr[3];
   int ii[3];
   if(!empty) for(ii[0]=0; ii[0] < this->fl_units[0]; ii[0] ++)
      for(ii[1]=0; ii[1] < this->fl_units[1]; ii[1] ++)
         for(ii[2]=0; ii[2] < this->fl_units[2]; ii[2] ++)
            for(int j=0; j < repl; j++)
               if(fill[ii[0]][ii[1]][ii[2]][j])
               {
                  for(int k=0; k < 3; k++)
                  {
                     tr[k] = off[k]
                           + (ii[k] + (1.0 + r->rnd())/3.0)*fl_unit[k];
                     if(tr[k] > box[k]) tr[k] -= box[k];
                     else if(tr[k] < 0.0) tr[k] += box[k];
                  }
                  tr[1] += j*box[1];
                  double tv = sqrt(3.0*T / FLUIDMASS);
                  double phi = 6.283185 * r->rnd();
                  double omega = 6.283185 * r->rnd();
                  double xxomga, yyomga;
                  if(I_xx_yy > 0.0)
                  {
                     xxomga = sqrt(r->rnd()*T/I_xx_yy);
                     if(r->rnd() > 0.5) xxomga *= -1.0;
                     yyomga = sqrt(r->rnd()*T/I_xx_yy);
                     if(r->rnd() > 0.5) yyomga *= -1.0;
                  }
                  else
                  {
                     xxomga = 0.0;
                     yyomga = 0.0;
                  }

                  xdr << id << " " << 1 << "\t" << tr[0]/SIG_REF
                      << " " << tr[1]/SIG_REF << " " << tr[2]/SIG_REF
                      << "\t" << tv*cos(phi)*cos(omega)/VEL_REF << " "
                      << tv*cos(phi)*sin(omega)/VEL_REF << " "
                      << tv*sin(phi)/VEL_REF << "\t1.0 0.0 0.0 0.0\t"
                      << xxomga/REFOMGA << " " << yyomga/REFOMGA
                      << " 0.0\n";
                  id++;
               }

   for(int j=0; j < repl; j++)
   {
      double layerU = ((flow == FLOW_COUETTE) && (j == 0))
                    ? U
                    : 0.0;
      gra.calculateVelocities(T, layerU);
      for(unsigned k=2; this->d+1 >= k; k++)
      {
         double yoffset = 0.5*h + (k-2.0)*Z;
         for(unsigned l=0; l < Ngraphene; l++)
         {
            tr[0] = gra.getX(l);
            tr[1] = gra.getY(l) + yoffset;
            tr[2] = gra.getZ(l);
            for(int m=0; m < 3; m++)
            {
               tr[m] += (0.004*r->rnd() - 0.002) * bondlength;
               if(tr[m] > box[m]) tr[m] -= box[m];
               else if(tr[m] < 0.0) tr[m] += box[m];
            }
            if(j == 1) tr[1] += box[1];

            xdr << id << " " << j*d + k << "\t" << tr[0]/SIG_REF
                << " " << tr[1]/SIG_REF << " " << tr[2]/SIG_REF
                << "\t" << gra.getVelocityX(l)/VEL_REF << " "
                << gra.getVelocityY(l)/VEL_REF << " "
                << gra.getVelocityZ(l)/VEL_REF
                << "\t1.0 0.0 0.0 0.0\t0.0 0.0 0.0\n";

            id++;
         }
      }
   }

   xdr.close();
   txt.close();
}

void Domain::writeNanotube(
   char* prefix, double a, bool empty, int format, double mu,
   double TAU, double U, bool original
) {
   cout << "Cannot create the nanotube - implementation missing.\n";
   exit(19);
}

void Domain::specifyGraphite(double rho, unsigned N)
{
   double V_id = (double)N/rho;
   if(this->flow == FLOW_COUETTE) V_id *= 0.5;
   cout << "Carbon-carbon bond length: " << bondlength
        << " * 0.05291772 nm.\n";
   cout << "Effective symmetry volume should approach " << V_id
        << " * 1.4818e-04 nm^3.\n";

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

   /*
    * x and z direction
    */
   double LxLz_id = V_id / this->eff[1];
   double C = 0.86602540378 * this->bondlength;  // sin(pi/3)
   double lxlz_id = LxLz_id / (6.0 * this->bondlength * C);
   unsigned lx = round(
                    sqrt(0.28867513 * lxlz_id)  // sin(pi/3) / 3
                 );
   if(lx == 0) lx = 1;
   unsigned lz = round(lxlz_id / lx);
   if(lz == 0) lz = 1;
   this->box[0] = 3.0*this->bondlength * (double)lx;
   this->eff[0] = this->box[0];
   this->off[0] = 0.0;
   this->box[2] = 2.0*C * (double)lz;
   this->eff[2] = this->box[2];
   this->off[2] = 0.0;

   /*
    * fluid unit box dimensions
    */
   double V_eff = this->eff[0] * this->eff[1] * this->eff[2];
   cout << "Symmetry volume " << box[0]*box[1]*box[2]
        << " * 1.4818e-04 nm^3, effectively " << V_eff
        << " * 1.4818e-04 nm^3.\n";
   double N_id = V_eff*rho;
   fl_units[1] = round(
                    pow(
                       (N_id * this->eff[1] * this->eff[1])
                            / (this->eff[0] * this->eff[2]), 1.0/3.0
                    )
                 );
   if(fl_units[1] == 0) fl_units[1] = 1;
   double bxbz_id = N_id / fl_units[1];
   fl_units[0] = round(sqrt(this->eff[0] * bxbz_id / this->eff[2]));
   if(fl_units[0] == 0) this->fl_units[0] = 1;
   fl_units[2] = ceil(bxbz_id / fl_units[0]);

   for(int i=0; i < 3; i++)
      this->fl_unit[i] = round(this->eff[i] / fl_units[i]);
   this->pfill = N_id / ((double)fl_units[0]*fl_units[1]*fl_units[2]);
}

void Domain::specifyNanotube(double rho, double m_per_n, unsigned N)
{
   cout << "Nanotubes are not yet implemented.\n";
   exit(16);
}
