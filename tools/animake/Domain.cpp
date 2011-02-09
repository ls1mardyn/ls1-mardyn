/*
 * GNU GPL version 2
 */

#include "Domain.h"
#include "Random.h"
#include <cmath>

#define DT 0.061240
#define PRECISION 8

Domain::Domain(
      int sp_fluid, double* sp_box, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      unsigned sp_N, double sp_T
) {
   this->fluid = sp_fluid;
   for(int d=0; d < 3; d++) this->box[d] = sp_box[d];
   this->SIG_REF = sp_SIG_REF;
   this->EPS_REF = sp_EPS_REF;
   this->REFMASS = sp_REFMASS;
   this->muVT = sp_muVT;
   this->N = sp_N;
   this->T = sp_T;
}

void Domain::write(char* prefix, int format, double mu)
{
   ofstream xdr, txt, buchholz;
   stringstream strstrm, txtstrstrm, buchholzstrstrm;
   if(format == FORMAT_BRANCH)
   {
      strstrm << prefix << ".xdr";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      strstrm << prefix << ".inp";
   }
   xdr.open(strstrm.str().c_str(), ios::trunc);
   if(format == FORMAT_BRANCH)
   {
      txtstrstrm << prefix << "_1R.txt";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      txtstrstrm << prefix << "_1R.cfg";
   }
   txt.open(txtstrstrm.str().c_str(), ios::trunc);
   if(format == FORMAT_BUCHHOLZ)
   {
      buchholzstrstrm << prefix << "_1R.xml";
      buchholz.open(buchholzstrstrm.str().c_str(), ios::trunc);

      /*
       * Gesamter Inhalt der Buchholz-Datei
       */
      buchholz << "<?xml version = \'1.0\' encoding = \'UTF-8\'?>\n<mardyn version=\""
               << TIME << "\">\n   <simulation type=\"MD\">\n      <input type=\"oldstyle\">"
               << prefix << "_1R.cfg</input>\n   </simulation>\n</mardyn>";
      buchholz.close();
   }

   Random* r = new Random();
   r->init(
      (int)(3162.3*box[0]) + (int)(31623.0*box[1]) - (int)(316.23*box[2])
   );
   double REFTIME = SIG_REF * sqrt(REFMASS / EPS_REF);
   double VEL_REF = SIG_REF / REFTIME;
   cout << "Velocity unit 1 = " << VEL_REF << " * 1620.34 m/s = "
        << 1620.34 * VEL_REF << " m/s.\n";
   double REFCARG = sqrt(EPS_REF * SIG_REF);
   cout << "Charge unit 1 = " << REFCARG << " e.\n";
   double DIP_REF = SIG_REF*REFCARG;
   double QDR_REF = SIG_REF*DIP_REF;
   double REFOMGA = 1.0 / REFTIME;

   unsigned fl_units[3];
   double fl_unit[3];
   double N_boxes = N / 3.0;
   fl_units[1] = round(
                    pow(
                       (N_boxes * box[1]*box[1])
                                / (this->box[0] * this->box[2]), 1.0/3.0
                    )
                 );
   if(fl_units[1] == 0) fl_units[1] = 1;
   double bxbz_id = N_boxes / fl_units[1];
   fl_units[0] = round(sqrt(this->box[0] * bxbz_id / this->box[2]));
   if(fl_units[0] == 0) fl_units[0] = 1;
   fl_units[2] = ceil(bxbz_id / fl_units[0]);
   for(int d=0; d < 3; d++) fl_unit[d] = box[d] / (double)fl_units[d];
   cout << "Unit cell dimensions: " << fl_unit[0] << " x " << fl_unit[1] << " x " << fl_unit[2] << ".\n";
   bool fill[fl_units[0]][fl_units[1]][fl_units[2]][3];
   unsigned slots = 3 * fl_units[0] * fl_units[1] * fl_units[2];
   /*
   double pfill = (double)N / slots;
   unsigned N1 = 0;
   for(unsigned i=0; i < fl_units[0]; i++)
      for(unsigned j=0; j < fl_units[1]; j++)
         for(unsigned k=0; k < fl_units[2]; k++)
            for(int d=0; d < 3; d++)
            {
               bool tfill = (pfill >= r->rnd());
               fill[i][j][k][d] = tfill;
               if(tfill) N1++;
            }
   */
   unsigned N1 = slots; 
   for(unsigned i=0; i < fl_units[0]; i++)
      for(unsigned j=0; j < fl_units[1]; j++)
         for(unsigned k=0; k < fl_units[2]; k++)
            for(int d=0; d < 3; d++) fill[i][j][k][d] = true;

   bool tswap;
   double pswap;
   for(int m=0; m < PRECISION; m++)
   {
      tswap = (N1 < N);
      pswap = ((double)N - (double)N1) / ((tswap? slots: 0) - (double)N1);
      // cout << "(N = " << N << ", N1 = " << N1 << ", tswap = " << tswap << ", pswap = " << pswap << ")\n";
      for(unsigned i=0; i < fl_units[0]; i++)
         for(unsigned j=0; j < fl_units[1]; j++)
            for(unsigned k=0; k < fl_units[2]; k++)
               for(int d=0; d < 3; d++)
                  if(pswap >= r->rnd())
                  {
                     if(fill[i][j][k][d]) N1--;
                     fill[i][j][k][d] = tswap;
                     if(tswap) N1++;
                  }
   }
   cout << "Filling " << N1 << " of 3*"
        << fl_units[0] << "*" << fl_units[1] << "*" << fl_units[2]
        << " = " << 3*fl_units[0]*fl_units[1]*fl_units[2]
        << " slots (ideally " << N << ").\n";

   double LJ_CUTOFF;
   double EL_CUTOFF;
   double FLUIDMASS, EPS_FLUID, SIG_FLUID, FLUIDLONG, QDR_FLUID;
   if(fluid == FLUID_AR)
   {
      FLUIDMASS = ARMASS;
      EPS_FLUID = EPS_AR;
      SIG_FLUID = SIG_AR;
      LJ_CUTOFF = CUT_AR;
   }
   else if(fluid == FLUID_CH4)
   {
      FLUIDMASS = CH4MASS;
      EPS_FLUID = EPS_CH4;
      SIG_FLUID = SIG_CH4;
      LJ_CUTOFF = CUT_CH4;
   }
   else if(fluid == FLUID_C2H6)
   {
      FLUIDMASS = C2H6MASS;
      EPS_FLUID = EPS_C2H6;
      SIG_FLUID = SIG_C2H6;
      FLUIDLONG = C2H6LONG;
      QDR_FLUID = QDR_C2H6;
      LJ_CUTOFF = CUT_C2H6;
   }
   else if(fluid == FLUID_N2)
   {
      FLUIDMASS = N2MASS;
      EPS_FLUID = EPS_N2;
      SIG_FLUID = SIG_N2;
      FLUIDLONG = N2LONG;
      QDR_FLUID = QDR_N2;
      LJ_CUTOFF = CUT_N2;
   }
   else if(fluid == FLUID_CO2)
   {
      FLUIDMASS = CO2MASS;
      EPS_FLUID = EPS_CO2;
      SIG_FLUID = SIG_CO2;
      FLUIDLONG = CO2LONG;
      QDR_FLUID = QDR_CO2;
      LJ_CUTOFF = CUT_CO2;
   }
   else if(fluid == FLUID_EOX)
   {
      FLUIDMASS = 2*CEOXMASS + OEOXMASS;
      LJ_CUTOFF = CUTLJEOX;
   }
   else
   {
      cout << "Unavailable fluid ID " << fluid << ".\n";
      exit(20);
   }

   if((fluid == FLUID_AR) || (fluid == FLUID_CH4))
   {
      EL_CUTOFF = LJ_CUTOFF;
   }
   else EL_CUTOFF = 1.2 * LJ_CUTOFF;

   xdr.precision(9);
   if(format == FORMAT_BRANCH)
   {
      xdr << "mardyn " << TIME << " tersoff\n"
          << "# mardyn input file, ls1 project\n"
          << "# written by animaker, the mesh generator\n";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      xdr << "mardyn trunk " << TIME << "\n";
   }

   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      xdr << "T\t" << T/EPS_REF << "\n";
      xdr << "t\t0.0\nL\t" << box[0]/SIG_REF << "\t"
          << box[1]/SIG_REF << "\t" << box[2]/SIG_REF
          << "\nC\t1\n";

      if((fluid == FLUID_AR) || (fluid == FLUID_CH4))
      {
         xdr << "1 0 0 0 0\n"  // LJ, C, Q, D, Tersoff
             << "0.0 0.0 0.0\t"
             << FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF << " "
             << SIG_FLUID/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\t0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_EOX)
      {
         xdr << "3 0 0 1 0\n";  // LJ, C, Q, D, Tersoff

         xdr << R0_C1EOX/SIG_REF << " " << R1_C1EOX/SIG_REF << " "
             << R2_C1EOX/SIG_REF << "\t" << CEOXMASS/REFMASS << " "
             << EPS_CEOX/EPS_REF << " " << SIG_CEOX/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << R0_C2EOX/SIG_REF << " " << R1_C2EOX/SIG_REF << " "
             << R2_C2EOX/SIG_REF << "\t" << CEOXMASS/REFMASS << " "
             << EPS_CEOX/EPS_REF << " " << SIG_CEOX/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << R0_O_EOX/SIG_REF << " " << R1_O_EOX/SIG_REF << " "
             << R2_O_EOX/SIG_REF << "\t" << OEOXMASS/REFMASS << " "
             << EPS_OEOX/EPS_REF << " " << SIG_OEOX/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n";

         xdr << R0DIPEOX/SIG_REF << " " << R1DIPEOX/SIG_REF << " "
             << R2DIPEOX/SIG_REF << "\t0.0 0.0 1.0 "
             << DIPOLEOX/DIP_REF;

         xdr << "\n0.0 0.0 0.0\n";
      }
      else
      {
         xdr << "2 0 1 0 0\n";  // LJ, C, Q, D, Tersoff

         xdr << "0.0 0.0 " << -0.5*FLUIDLONG/SIG_REF << "\t"
             << 0.5*FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF
             << " " << SIG_FLUID/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << "0.0 0.0 " << +0.5*FLUIDLONG/SIG_REF << "\t"
             << 0.5*FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF
             << " " << SIG_FLUID/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n";

         xdr << "0.0 0.0 0.0\t0.0 0.0 1.0\t" << QDR_FLUID/QDR_REF;

         xdr << "\n0.0 0.0 0.0\n";
      }
      xdr << "1.0e+10\n";

      xdr << "N" << "\t" << N1 << "\nM" << "\t" << "ICRVQD\n\n";

      txt.precision(6);
      txt << "mardynconfig\n# \ntimestepLength\t" << DT/REFTIME
          << "\ncutoffRadius\t" << EL_CUTOFF/SIG_REF
          << "\nLJCutoffRadius\t" << LJ_CUTOFF/SIG_REF;
   }
   if(format == FORMAT_BRANCH)
   {
      txt << "\ntersoffCutoffRadius\t"
          <<  LJ_CUTOFF/(3.0*SIG_REF);
   }
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt << "\ninitCanonical\t50000\n";
      if(muVT)
      {
         txt.precision(9);
         txt << "chemicalPotential " << mu/EPS_REF
             << " component 1 control 0.0 0.0 0.0 to "
             << this->box[0]/SIG_REF << " " << this->box[1]/SIG_REF
             << " " << this->box[2]/SIG_REF << " conduct "
             << 1 + (int)round(0.003 * N)
             << " tests every 3 steps\n";
	 txt << "planckConstant\t" 
	     << sqrt(6.28319 * T/EPS_REF) << "\n";  // sqrt(2 pi kT)
         txt << "initGrandCanonical\t100000\n";
      }
      txt.precision(5);
      txt << "initStatistics\t150000\n";
   }
   if(format == FORMAT_BRANCH)
   {
      txt << "phaseSpaceFile\t" << prefix << ".xdr\n";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      txt << "phaseSpaceFile\tOldStyle\t" << prefix << ".inp\n";
   }
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt << "datastructure\tLinkedCells\t1\noutput\t"
          << "ResultWriter\t100\t" << prefix
          << "_1R\noutput\tXyzWriter\t10000\t" << prefix
          << "_1R.buxyz\n";
   }
   txt.close();

   double I[3];
   for(int d=0; d < 3; d++) I[d] = 0.0;
   if(fluid == FLUID_EOX)
   {
      I[0] = I_XX_EOX;
      I[1] = I_YY_EOX;
      I[2] = I_ZZ_EOX;
   }
   else if(!((fluid == FLUID_AR) || (fluid == FLUID_CH4)))
   {
      I[0] = 0.25 * FLUIDMASS * FLUIDLONG * FLUIDLONG;
      I[1] = I[0];
   }

   unsigned id = 1;
   double tr[3];
   unsigned ii[3];
   for(ii[0]=0; ii[0] < fl_units[0]; (ii[0]) ++)
      for(ii[1]=0; ii[1] < fl_units[1]; (ii[1]) ++)
         for(ii[2]=0; ii[2] < fl_units[2]; (ii[2]) ++)
            for(int d=0; d < 3; d++)
            {
               if(fill[ ii[0] ][ ii[1] ][ ii[2] ][ d ])
               {
                  for(int k=0; k < 3; k++)
                  {
                     tr[k] = fl_unit[k] * (
                                ii[k] + 0.02*r->rnd() + ((k == d)? 0.24: 0.74)
                             );
                  }
                  for(int k=0; k < 3; k++)
                  {
                     if(tr[k] > box[k]) tr[k] -= box[k];
                     else if(tr[k] < 0.0) tr[k] += box[k];
                  }
                  double tv = sqrt(3.0*T / FLUIDMASS);
                  double phi = 6.283185 * r->rnd();
                  double omega = 6.283185 * r->rnd();
                  double w[3];
                  for(int d=0; d < 3; d++)
                     w[d] = (I[d] == 0)? 0.0: ((r->rnd() > 0.5)? 1: -1) * sqrt(2.0*r->rnd()*T / I[d]);
                  // xdr << "(" << ii[0] << "/" << ii[1] << "/" << ii[2] << "), j = " << j << ", d = " << d << ":\t";
                  if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
                  {
                     xdr << id << " " << 1 << "\t" << tr[0]/SIG_REF
                         << " " << tr[1]/SIG_REF << " " << tr[2]/SIG_REF
                         << "\t" << tv*cos(phi)*cos(omega)/VEL_REF << " "
                         << tv*cos(phi)*sin(omega)/VEL_REF << " "
                         << tv*sin(phi)/VEL_REF << "\t1.0 0.0 0.0 0.0\t"
                         << w[0]/REFOMGA << " " << w[1]/REFOMGA << " "
                         << w[2]/REFOMGA << "\n";
                  }
                  id++;
               }
               else xdr << "\n";
            }
   xdr << "\n";

   xdr.close();
}

