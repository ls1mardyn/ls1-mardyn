#include "Domain.h"
#include "Random.h"
#include <cmath>

#define DT 0.061240
#define PRECISION 8
#define BINS 1024

Domain::Domain(
      int sp_fluid, int sp_fluid2, double* sp_box, double sp_SIG_REF,
      double sp_EPS_REF, double sp_REFMASS, bool sp_muVT,
      unsigned sp_N, double sp_T, double sp_ETAF, double sp_XIF
) {
   this->fluid = sp_fluid;
   this->fluid2 = sp_fluid2;
   for(int d=0; d < 3; d++) this->box[d] = sp_box[d];
   this->SIG_REF = sp_SIG_REF;
   this->EPS_REF = sp_EPS_REF;
   this->REFMASS = sp_REFMASS;
   this->muVT = sp_muVT;
   this->N = sp_N;
   this->T = sp_T;
   this->ETAF = sp_ETAF;
   this->XIF = sp_XIF;
}

void Domain::write(char* prefix, int format, double mu, double x)
{
   ofstream xdr, txt, buchholz;
   stringstream strstrm, txtstrstrm, buchholzstrstrm;
   if(format == FORMAT_BRANCH)
   {
      strstrm << prefix << ".xdr";
   }
   if((format == FORMAT_BERNREUTHER) || (format == FORMAT_BUCHHOLZ))
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
   if(format == FORMAT_BERNREUTHER)
   {
      txtstrstrm << prefix << "_1R.xml";
   }
   txt.open(txtstrstrm.str().c_str(), ios::trunc);
   if(format == FORMAT_BERNREUTHER)
   {
      txt << "<?xml version=\'1.0\' encoding=\'UTF-8\'?>\n<mardyn version=\""
          << TIME << "\">\n";
   }
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
   if(format == FORMAT_BERNREUTHER)
   {
      txt << "<refunits type=\"SI\"><length unit=\"nm\">" << SIG_REF * 0.0529177 << "</length><mass unit=\"u\">" << REFMASS * 1000.0 << "</mass><energy unit=\"eV\">" << EPS_REF * 27.2126 << "</energy></refunits>\n";
   }

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
   else if(fluid == FLUID_JES)
   {
      FLUIDMASS = OJESMASS + 2.0*HJESMASS;
      EPS_FLUID = EPS_OJES;
      SIG_FLUID = SIG_OJES;
      FLUIDLONG = JES_LONG;
      LJ_CUTOFF = 4.0*SIG_FLUID + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_MER)
   {
      FLUIDMASS = CMERMASS + 2.0*OMERMASS;
      FLUIDLONG = MER_LONG;
      QDR_FLUID = QDR_CMER;
      LJ_CUTOFF = 4.0*SIG_OMER + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_TOL)
   {
      FLUIDMASS = CH3TOLMASS + CTRTOLMASS + 5.0*CH_TOLMASS;
      FLUIDLONG = TOL_LONG;
      LJ_CUTOFF = 4.0*SIG_CH_TOL + 0.5*FLUIDLONG;
   }
   else if(fluid == FLUID_VEG)
   {
      FLUIDMASS = OVEGMASS + 2.0*HVEGMASS;
      EPS_FLUID = EPS_OVEG;
      SIG_FLUID = SIG_OVEG;
      FLUIDLONG = VEG_LONG;
      LJ_CUTOFF = 4.0*SIG_FLUID + 0.5*FLUIDLONG;
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

   double LJ_CUTOFF2 = 0.0;
   double EL_CUTOFF2 = 0.0;
   double FLUIDMASS2, EPS_FLUID2, SIG_FLUID2, FLUIDLONG2, QDR_FLUID2;
   if(fluid2 == FLUID_AR)
   {
      FLUIDMASS2 = ARMASS;
      EPS_FLUID2 = EPS_AR;
      SIG_FLUID2 = SIG_AR;
      LJ_CUTOFF2 = 2.5*SIG_FLUID2;
   }
   else if(fluid2 == FLUID_CH4)
   {
      FLUIDMASS2 = CH4MASS;
      EPS_FLUID2 = EPS_CH4;
      SIG_FLUID2 = SIG_CH4;
      LJ_CUTOFF2 = 2.5*SIG_FLUID2;
   }
   else if(fluid2 == FLUID_C2H6)
   {
      FLUIDMASS2 = C2H6MASS;
      EPS_FLUID2 = EPS_C2H6;
      SIG_FLUID2 = SIG_C2H6;
      FLUIDLONG2 = C2H6LONG;
      QDR_FLUID2 = QDR_C2H6;
      LJ_CUTOFF2 = 4.0*SIG_FLUID2 + 0.5*FLUIDLONG2;
   }
   else if(fluid2 == FLUID_N2)
   {
      FLUIDMASS2 = N2MASS;
      EPS_FLUID2 = EPS_N2;
      SIG_FLUID2 = SIG_N2;
      FLUIDLONG2 = N2LONG;
      QDR_FLUID2 = QDR_N2;
      LJ_CUTOFF2 = 4.0*SIG_FLUID2 + 0.5*FLUIDLONG2;
   }
   else if(fluid2 == FLUID_CO2)
   {
      FLUIDMASS2 = CO2MASS;
      EPS_FLUID2 = EPS_CO2;
      SIG_FLUID2 = SIG_CO2;
      FLUIDLONG2 = CO2LONG;
      QDR_FLUID2 = QDR_CO2;
      LJ_CUTOFF2 = 4.0*SIG_FLUID2 + 0.5*FLUIDLONG2;
   }
   else if(fluid2 == FLUID_EOX)
   {
      FLUIDMASS2 = 2*CEOXMASS + OEOXMASS;
      LJ_CUTOFF2 = CUTLJEOX;
   }
   else if(fluid2 == FLUID_JES)
   {
      FLUIDMASS2 = OJESMASS + 2.0*HJESMASS;
      EPS_FLUID2 = EPS_OJES;
      SIG_FLUID2 = SIG_OJES;
      FLUIDLONG2 = JES_LONG;
      LJ_CUTOFF2 = 4.0*SIG_FLUID2 + 0.5*FLUIDLONG2;
   }
   else if(fluid2 == FLUID_MER)
   {
      FLUIDMASS2 = CMERMASS + 2.0*OMERMASS;
      FLUIDLONG2 = MER_LONG;
      QDR_FLUID2 = QDR_CMER;
      LJ_CUTOFF2 = 4.0*SIG_OMER + 0.5*FLUIDLONG2;
   }
   else if(fluid2 == FLUID_TOL)
   {
      FLUIDMASS2 = CH3TOLMASS + CTRTOLMASS + 5.0*CH_TOLMASS;
      FLUIDLONG2 = TOL_LONG;
      LJ_CUTOFF2 = 4.0*SIG_CH_TOL + 0.5*FLUIDLONG2;
   }
   else if(fluid2 == FLUID_VEG)
   {
      FLUIDMASS2 = OVEGMASS + 2.0*HVEGMASS;
      EPS_FLUID2 = EPS_OVEG;
      SIG_FLUID2 = SIG_OVEG;
      FLUIDLONG2 = VEG_LONG;
      LJ_CUTOFF2 = 4.0*SIG_FLUID2 + 0.5*FLUIDLONG2;
   }
   else fluid2 = FLUID_NIL;

   if((fluid2 == FLUID_AR) || (fluid2 == FLUID_CH4))
   {
      EL_CUTOFF2 = LJ_CUTOFF2;
   }
   else if(fluid2 != FLUID_NIL) EL_CUTOFF2 = 1.2*LJ_CUTOFF2;

   if(LJ_CUTOFF2 > LJ_CUTOFF) LJ_CUTOFF = LJ_CUTOFF2;
   if(EL_CUTOFF2 > EL_CUTOFF) EL_CUTOFF = EL_CUTOFF2;
   for(int dim=0; dim < 3; dim++)
   {
      if(0.5*box[dim] < EL_CUTOFF) EL_CUTOFF = 0.5*box[dim];
      if(0.5*box[dim] < LJ_CUTOFF) LJ_CUTOFF = 0.5*box[dim];
   }

   unsigned fluidcomp = (fluid2 != FLUID_NIL)? 2: 1;

   xdr.precision(9);
   if(format == FORMAT_BRANCH)
   {
      xdr << "mardyn " << TIME << " tersoff\n"
          << "# mardyn input file, ls1 project\n"
          << "# written by animake, the mesh generator\n";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      xdr << "mardyn trunk " << TIME << "\n";
   }

   if(format == FORMAT_BERNREUTHER)
   {
      txt << "<simulation type=\"MD\">\n";

      txt << "<integrator type=\"Leapfrog\"><timestep unit=\"reduced\">"
          << DT/REFTIME << "</timestep></integrator>\n"
          << "<run><currenttime>0</currenttime><production><steps>100000000</steps></production></run>\n"
          << "<algorithm><parallelisation type=\"DomainDecomposition\"></parallelisation><datastructure type=\"LinkedCells\"><cellsInCutoffRadius>1</cellsInCutoffRadius></datastructure><cutoffs type=\"CenterOfMass\"><radiusLJ unit=\"reduced\">"
          << LJ_CUTOFF/SIG_REF << "</radiusLJ></cutoffs><electrostatic type=\"ReactionField\"><epsilon>1.0e+10</epsilon></electrostatic></algorithm>\n"
          << "<output>\n"
          << "<outputplugin name=\"Resultwriter\"><writefrequency>100</writefrequency><outputprefix>"
          << prefix << "_1R</outputprefix></outputplugin>\n"
          << "<outputplugin name=\"CheckpointWriter\"><writefrequency>10000</writefrequency><outputprefix>"
          << prefix << "_1R</outputprefix></outputplugin>\n"
          << "</output>\n"
          << "<ensemble type=\"NVT\"><temperature unit=\"reduced\">"
          << T/EPS_REF << "</temperature>\n"
          << "<domain type=\"box\"><lx>" << box[0]/SIG_REF
          << "</lx><ly>" << box[1]/SIG_REF << "</ly><lz>"
          << box[2]/SIG_REF << "</lz></domain>\n"
          << "<components>\n";

      txt << "<moleculetype id=\"1\">\n";
      if((fluid == FLUID_AR) || (fluid == FLUID_CH4))
      {
         txt << "<site type=\"LJ126\" id=\"1\"><coord><x>0</x><y>0</y><z>0</z></coord><mass>"
             << FLUIDMASS/REFMASS << "</mass><sigma>"
             << SIG_FLUID/SIG_REF << "</sigma><epsilon>"
             << EPS_FLUID/EPS_REF
             << "</epsilon><shifted>false</shifted></site>\n";
      }
      else if(fluid == FLUID_EOX)
      {
         txt << "<site type=\"LJ126\" id=\"1\"><coord><x>"
             << R0_C1EOX/SIG_REF << "</x><y>" << R1_C1EOX/SIG_REF
             << "</y><z>" << R2_C1EOX/SIG_REF << "</z></coord><mass>"
             << CEOXMASS/REFMASS << "</mass><sigma>"
             << SIG_CEOX/SIG_REF << "</sigma><epsilon>"
             << EPS_CEOX/EPS_REF
             << "</epsilon><shifted>false</shifted></site>\n"
             << "<site type=\"LJ126\" id=\"2\"><coord><x>"
             << R0_C2EOX/SIG_REF << "</x><y>" << R1_C2EOX/SIG_REF
             << "</y><z>" << R2_C2EOX/SIG_REF << "</z></coord><mass>"
             << CEOXMASS/REFMASS << "</mass><sigma>"
             << SIG_CEOX/SIG_REF << "</sigma><epsilon>"
             << EPS_CEOX/EPS_REF
             << "</epsilon><shifted>false</shifted></site>\n"
             << "<site type=\"LJ126\" id=\"3\"><coord><x>"
             << R0_O_EOX/SIG_REF << "</x><y>" << R1_O_EOX/SIG_REF
             << "</y><z>" << R2_O_EOX/SIG_REF << "</z></coord><mass>"
             << OEOXMASS/REFMASS << "</mass><sigma>"
             << SIG_OEOX/SIG_REF << "</sigma><epsilon>"
             << EPS_OEOX/EPS_REF
             << "</epsilon><shifted>false</shifted></site>\n"
             << "<site type=\"Dipol\" id=\"4\"><coord><x>"
             << R0DIPEOX/SIG_REF << "</x><y>" << R1DIPEOX/SIG_REF
             << "</y><z>" << R2DIPEOX/SIG_REF
             << "</z></coord><dipolemoment><orientation><x>0</x><y>0</y><z>1</z></orientation><abs>"
             << DIPOLEOX/DIP_REF << "</abs></dipolemoment></site>\n";
      }
      txt << "</moleculetype>\n";

      txt << "</components>\n"
          << "<phasespacepoint><file type=\"ASCII\">" << prefix 
          << ".inp</file></phasespacepoint>\n"
          << "</ensemble>\n";
   }
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      xdr << "T\t" << T/EPS_REF << "\n";
      xdr << "t\t0.0\nL\t" << box[0]/SIG_REF << "\t"
          << box[1]/SIG_REF << "\t" << box[2]/SIG_REF
          << "\nC\t" << fluidcomp << "\n";
   }
   if(format == FORMAT_BRANCH)
   {
      if((fluid == FLUID_AR) || (fluid == FLUID_CH4))
      {
         xdr << "1 0 0 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid == FLUID_EOX)
      {
         xdr << "3 0 0 1 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid == FLUID_JES)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid == FLUID_MER)
      {
         xdr << "3 0 1 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid == FLUID_TOL)
      {
         xdr << "7 0 5 1 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid == FLUID_VEG)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else
      {
         xdr << "2 0 1 0 0\n";  // LJ, C, Q, D, Tersoff
      }
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      if((fluid == FLUID_AR) || (fluid == FLUID_CH4))
      {
         xdr << "1 0 0 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid == FLUID_EOX)
      {
         xdr << "3 0 1 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid == FLUID_JES)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid == FLUID_MER)
      {
         xdr << "3 0 0 1 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid == FLUID_TOL)
      {
         xdr << "7 0 1 5 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid == FLUID_VEG)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else
      {
         xdr << "2 0 0 1 0\n";  // LJ, C, D, Q, Tersoff
      }
   }
   
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      if((fluid == FLUID_AR) || (fluid == FLUID_CH4))
      {
         xdr << "0.0 0.0 0.0\t"
             << FLUIDMASS/REFMASS << " " << EPS_FLUID/EPS_REF << " "
             << SIG_FLUID/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\t0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_EOX)
      {
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
      else if(fluid == FLUID_JES)
      {
         xdr << R0_O_JES/SIG_REF << " " << R1_O_JES/SIG_REF << " " << R2_O_JES/SIG_REF << "\t"
             << OJESMASS/REFMASS << " " << EPS_OJES/EPS_REF << " " << SIG_OJES/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << R0_H1JES/SIG_REF << " " << R1_H1JES/SIG_REF << " " << R2_H1JES/SIG_REF << "\t"
             << HJESMASS/REFMASS << " " << CHG_HJES/REFCARG << "\n";
         xdr << R0_H2JES/SIG_REF << " " << R1_H2JES/SIG_REF << " " << R2_H2JES/SIG_REF << "\t"
             << HJESMASS/REFMASS << " " << CHG_HJES/REFCARG << "\n";
         xdr << R0_E_JES/SIG_REF << " " << R1_E_JES/SIG_REF << " " << R2_E_JES/SIG_REF << "\t"
             << "0.0 " << CHG_EJES/REFCARG << "\n";

         xdr << "0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_MER)
      {
         xdr << "0.0 0.0 " << -0.5*MER_LONG/SIG_REF << "\t"
             << OMERMASS/REFMASS << " " << EPS_OMER/EPS_REF
             << " " << SIG_OMER/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << "0.0 0.0 " << +0.5*MER_LONG/SIG_REF << "\t"
             << OMERMASS/REFMASS << " " << EPS_OMER/EPS_REF
             << " " << SIG_OMER/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << "0.0 0.0 0.0\t" << CMERMASS/REFMASS << " "
             << EPS_CMER/EPS_REF << " " << SIG_CMER/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n";

         xdr << "0.0 0.0 0.0\t0.0 0.0 1.0\t" << QDR_CMER/QDR_REF;

         xdr << "\n0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_TOL)
      {
         xdr << "0.0 " << R1_CH3_TOL/SIG_REF << " "
             << R2_CH3_TOL/SIG_REF << "\t" << CH3TOLMASS/REFMASS
             << " " << EPS_CH3TOL/EPS_REF << " "
             << SIG_CH3TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << "0.0 " << R1_CTR_TOL/SIG_REF << " "
             << R2_CTR_TOL/SIG_REF << "\t" << CTRTOLMASS/REFMASS
             << " " << EPS_CTRTOL/EPS_REF << " "
             << SIG_CTRTOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << "0.0 " << R1_CHA_TOL/SIG_REF << " "
             << R2_CHA_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHB_TOL/SIG_REF << " "
             << R2_CHB_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHC_TOL/SIG_REF << " "
             << R2_CHC_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHD_TOL/SIG_REF << " "
             << R2_CHD_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHE_TOL/SIG_REF << " "
             << R2_CHE_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         if(format == FORMAT_BUCHHOLZ)
         {
            xdr << "0.0 " << R1_CTR_TOL/SIG_REF << " "
                << R2_CTR_TOL/SIG_REF << "\t0.0 0.0 -1.0 "
                << DIP_CTRTOL/DIP_REF << "\n";
         }

         xdr << "0.0 " << R1_CHA_TOL/SIG_REF << " "
             << R2_CHA_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHB_TOL/SIG_REF << " "
             << R2_CHB_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHC_TOL/SIG_REF << " "
             << R2_CHC_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHD_TOL/SIG_REF << " "
             << R2_CHD_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHE_TOL/SIG_REF << " "
             << R2_CHE_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";

         if(format == FORMAT_BRANCH)
         {
            xdr << "0.0 " << R1_CTR_TOL/SIG_REF << " "
                << R2_CTR_TOL/SIG_REF << "\t0.0 0.0 -1.0 "
                << DIP_CTRTOL/DIP_REF << "\n";
         }

         xdr << "0.0 0.0 0.0\n";
      }
      else if(fluid == FLUID_VEG)
      {
         xdr << R0_O_VEG/SIG_REF << " " << R1_O_VEG/SIG_REF << " " << R2_O_VEG/SIG_REF << "\t"
             << OVEGMASS/REFMASS << " " << EPS_OVEG/EPS_REF << " " << SIG_OVEG/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << R0_H1VEG/SIG_REF << " " << R1_H1VEG/SIG_REF << " " << R2_H1VEG/SIG_REF << "\t"
             << HVEGMASS/REFMASS << " " << CHG_HVEG/REFCARG << "\n";
         xdr << R0_H2VEG/SIG_REF << " " << R1_H2VEG/SIG_REF << " " << R2_H2VEG/SIG_REF << "\t"
             << HVEGMASS/REFMASS << " " << CHG_HVEG/REFCARG << "\n";
         xdr << R0_E_VEG/SIG_REF << " " << R1_E_VEG/SIG_REF << " " << R2_E_VEG/SIG_REF << "\t"
             << "0.0 " << CHG_EVEG/REFCARG << "\n";

         xdr << "0.0 0.0 0.0\n";
      }
      else
      {
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
   }

   if(format == FORMAT_BRANCH)
   {
      if((fluid2 == FLUID_AR) || (fluid2 == FLUID_CH4))
      {
         xdr << "1 0 0 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid2 == FLUID_EOX)
      {
         xdr << "3 0 0 1 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid2 == FLUID_JES)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid2 == FLUID_MER)
      {
         xdr << "3 0 1 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid2 == FLUID_TOL)
      {
         xdr << "7 0 5 1 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid2 == FLUID_VEG)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, Q, D, Tersoff
      }
      else if(fluid2 != FLUID_NIL)
      {
         xdr << "2 0 1 0 0\n";  // LJ, C, Q, D, Tersoff
      }
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      if((fluid2 == FLUID_AR) || (fluid2 == FLUID_CH4))
      {
         xdr << "1 0 0 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid2 == FLUID_EOX)
      {
         xdr << "3 0 1 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid2 == FLUID_JES)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid2 == FLUID_MER)
      {
         xdr << "3 0 0 1 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid2 == FLUID_TOL)
      {
         xdr << "7 0 1 5 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid2 == FLUID_VEG)
      {
         xdr << "1 3 0 0 0\n";  // LJ, C, D, Q, Tersoff
      }
      else if(fluid2 != FLUID_NIL)
      {
         xdr << "2 0 0 1 0\n";  // LJ, C, D, Q, Tersoff
      }
   }

   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      if((fluid2 == FLUID_AR) || (fluid2 == FLUID_CH4))
      {
         xdr << "0.0 0.0 0.0\t"
             << FLUIDMASS2/REFMASS << " " << EPS_FLUID2/EPS_REF << " "
             << SIG_FLUID2/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\t0.0 0.0 0.0\n";
      }
      else if(fluid2 == FLUID_EOX)
      {
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
      else if(fluid2 == FLUID_JES)
      {
         xdr << R0_O_JES/SIG_REF << " " << R1_O_JES/SIG_REF << " " << R2_O_JES/SIG_REF << "\t"
             << OJESMASS/REFMASS << " " << EPS_OJES/EPS_REF << " " << SIG_OJES/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << R0_H1JES/SIG_REF << " " << R1_H1JES/SIG_REF << " " << R2_H1JES/SIG_REF << "\t"
             << HJESMASS/REFMASS << " " << CHG_HJES/REFCARG << "\n";
         xdr << R0_H2JES/SIG_REF << " " << R1_H2JES/SIG_REF << " " << R2_H2JES/SIG_REF << "\t"
             << HJESMASS/REFMASS << " " << CHG_HJES/REFCARG << "\n";
         xdr << R0_E_JES/SIG_REF << " " << R1_E_JES/SIG_REF << " " << R2_E_JES/SIG_REF << "\t"
             << "0.0 " << CHG_EJES/REFCARG << "\n";

         xdr << "0.0 0.0 0.0\n";
      }
      else if(fluid2 == FLUID_MER)
      {
         xdr << "0.0 0.0 " << -0.5*MER_LONG/SIG_REF << "\t"
             << OMERMASS/REFMASS << " " << EPS_OMER/EPS_REF
             << " " << SIG_OMER/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << "0.0 0.0 " << +0.5*MER_LONG/SIG_REF << "\t"
             << OMERMASS/REFMASS << " " << EPS_OMER/EPS_REF
             << " " << SIG_OMER/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << "0.0 0.0 0.0\t" << CMERMASS/REFMASS << " "
             << EPS_CMER/EPS_REF << " " << SIG_CMER/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n";

         xdr << "0.0 0.0 0.0\t0.0 0.0 1.0\t" << QDR_CMER/QDR_REF;

         xdr << "\n0.0 0.0 0.0\n";
      }
      else if(fluid2 == FLUID_TOL)
      {
         xdr << "0.0 " << R1_CH3_TOL/SIG_REF << " "
             << R2_CH3_TOL/SIG_REF << "\t" << CH3TOLMASS/REFMASS
             << " " << EPS_CH3TOL/EPS_REF << " "
             << SIG_CH3TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << "0.0 " << R1_CTR_TOL/SIG_REF << " "
             << R2_CTR_TOL/SIG_REF << "\t" << CTRTOLMASS/REFMASS
             << " " << EPS_CTRTOL/EPS_REF << " "
             << SIG_CTRTOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << "0.0 " << R1_CHA_TOL/SIG_REF << " "
             << R2_CHA_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHB_TOL/SIG_REF << " "
             << R2_CHB_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHC_TOL/SIG_REF << " "
             << R2_CHC_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHD_TOL/SIG_REF << " "
             << R2_CHD_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";
         xdr << "0.0 " << R1_CHE_TOL/SIG_REF << " "
             << R2_CHE_TOL/SIG_REF << "\t" << CH_TOLMASS/REFMASS
             << " " << EPS_CH_TOL/EPS_REF << " "
             << SIG_CH_TOL/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         if(format == FORMAT_BUCHHOLZ)
         {
            xdr << "0.0 " << R1_CTR_TOL/SIG_REF << " "
                << R2_CTR_TOL/SIG_REF << "\t0.0 0.0 -1.0 "
                << DIP_CTRTOL/DIP_REF << "\n";
         }

         xdr << "0.0 " << R1_CHA_TOL/SIG_REF << " "
             << R2_CHA_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHB_TOL/SIG_REF << " "
             << R2_CHB_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHC_TOL/SIG_REF << " "
             << R2_CHC_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHD_TOL/SIG_REF << " "
             << R2_CHD_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";
         xdr << "0.0 " << R1_CHE_TOL/SIG_REF << " "
             << R2_CHE_TOL/SIG_REF << "\t" << "1.0 0.0 0.0\t"
             << QDR_CH_TOL/QDR_REF << "\n";

         if(format == FORMAT_BRANCH)
         {
            xdr << "0.0 " << R1_CTR_TOL/SIG_REF << " "
                << R2_CTR_TOL/SIG_REF << "\t0.0 0.0 -1.0 "
                << DIP_CTRTOL/DIP_REF << "\n";
         }

         xdr << "0.0 0.0 0.0\n";
      }
      else if(fluid2 == FLUID_VEG)
      {
         xdr << R0_O_VEG/SIG_REF << " " << R1_O_VEG/SIG_REF << " " << R2_O_VEG/SIG_REF << "\t"
             << OVEGMASS/REFMASS << " " << EPS_OVEG/EPS_REF << " " << SIG_OVEG/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF << " 0";
         xdr << "\n";

         xdr << R0_H1VEG/SIG_REF << " " << R1_H1VEG/SIG_REF << " " << R2_H1VEG/SIG_REF << "\t"
             << HVEGMASS/REFMASS << " " << CHG_HVEG/REFCARG << "\n";
         xdr << R0_H2VEG/SIG_REF << " " << R1_H2VEG/SIG_REF << " " << R2_H2VEG/SIG_REF << "\t"
             << HVEGMASS/REFMASS << " " << CHG_HVEG/REFCARG << "\n";
         xdr << R0_E_VEG/SIG_REF << " " << R1_E_VEG/SIG_REF << " " << R2_E_VEG/SIG_REF << "\t"
             << "0.0 " << CHG_EVEG/REFCARG << "\n";

         xdr << "0.0 0.0 0.0\n";
      }
      else if(fluid2 != FLUID_NIL)
      {
         xdr << "0.0 0.0 " << -0.5*FLUIDLONG2/SIG_REF << "\t"
             << 0.5*FLUIDMASS2/REFMASS << " " << EPS_FLUID2/EPS_REF
             << " " << SIG_FLUID2/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n" << "0.0 0.0 " << +0.5*FLUIDLONG2/SIG_REF << "\t"
             << 0.5*FLUIDMASS2/REFMASS << " " << EPS_FLUID2/EPS_REF
             << " " << SIG_FLUID2/SIG_REF;
         if(format == FORMAT_BUCHHOLZ) xdr << "\t" << LJ_CUTOFF/SIG_REF << " 0";
         xdr << "\n";

         xdr << "0.0 0.0 0.0\t0.0 0.0 1.0\t" << QDR_FLUID2/QDR_REF;

         xdr << "\n0.0 0.0 0.0\n";
      }

      if(fluid2 != FLUID_NIL)
      {
         xdr << XIF << "\t" << ETAF << "\n";
      }
      xdr << "1.0e+10\n";

      txt.precision(6);
      txt << "mardynconfig\n# \ntimestepLength\t" << DT/REFTIME
          << "\ncutoffRadius\t" << EL_CUTOFF/SIG_REF
          << "\nLJCutoffRadius\t" << LJ_CUTOFF/SIG_REF;
   }
   if((format == FORMAT_BERNREUTHER) || (format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      xdr << "N" << "\t" << N1 << "\nM" << "\t" << "ICRVQD\n";
   }
   if(format == FORMAT_BRANCH)
   {
      txt << "\ntersoffCutoffRadius\t"
          <<  LJ_CUTOFF/(4.0*SIG_REF);
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
      txt << "phaseSpaceFile\tOldStyle\t" << prefix << ".inp\nparallelization\tDomainDecomposition\n";
   }
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt << "datastructure\tLinkedCells\t1\noutput\t"
          << "ResultWriter\t100\t" << prefix
          << "_1R\noutput\tXyzWriter\t10000\t" << prefix
          << "_1R.buxyz\n"
          << "RDF\t" << (EL_CUTOFF - 4.0*FLUIDLONG)/(SIG_REF * (double)BINS) << " "
          << BINS << "\nRDFOutputTimesteps\t150000\nRDFOutputPrefix\t" << prefix << "_1R\n";
   }
   if(format == FORMAT_BERNREUTHER)
   {
      txt << "</simulation></mardyn>\n";
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
   else if(fluid == FLUID_JES)
   {
      I[0] = I_XX_JES;
      I[1] = I_YY_JES;
      I[2] = I_ZZ_JES;
   }
   else if(fluid == FLUID_TOL)
   {
      I[0] = I_XX_TOL;
      I[1] = I_YY_TOL;
      I[2] = I_ZZ_TOL;
   }
   else if(fluid == FLUID_VEG)
   {
      I[0] = I_XX_VEG;
      I[1] = I_YY_VEG;
      I[2] = I_ZZ_VEG;
   }
   else if(!((fluid == FLUID_AR) || (fluid == FLUID_CH4)))
   {
      I[0] = 0.25 * FLUIDMASS * FLUIDLONG * FLUIDLONG;
      I[1] = I[0];
   }

   double I2[3];
   for(int k=0; k < 3; k++) I2[k] = 0.0;
   if(fluid2 == FLUID_EOX)
   {
      I2[0] = I_XX_EOX;
      I2[1] = I_YY_EOX;
      I2[2] = I_ZZ_EOX;
   }
   else if(fluid2 == FLUID_JES)
   {
      I2[0] = I_XX_JES;
      I2[1] = I_YY_JES;
      I2[2] = I_ZZ_JES;
   }
   else if(fluid2 == FLUID_TOL)
   {
      I2[0] = I_XX_TOL;
      I2[1] = I_YY_TOL;
      I2[2] = I_ZZ_TOL;
   }
   else if(fluid2 == FLUID_VEG)
   {
      I2[0] = I_XX_VEG;
      I2[1] = I_YY_VEG;
      I2[2] = I_ZZ_VEG;
   }
   else if(!((fluid2 == FLUID_AR) || (fluid2 == FLUID_CH4)))
   {
      I2[0] = 0.25 * FLUIDMASS2 * FLUIDLONG2 * FLUIDLONG2;
      I2[1] = I2[0];
   }

   unsigned id = 1;
   double tr[3];
   unsigned ii[3];
   double Nf[2];
   Nf[0] = (1.0 - x) * (double)N1 + 0.001;
   Nf[1] = x * (double)N1 + 0.001;
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
                  unsigned cid;
                  if(fluidcomp == 1) cid = 1;
                  else
                  {
                     cid = (r->rnd() > (Nf[0] / (Nf[0] + Nf[1])))? 2: 1;
                     --Nf[cid - 1];
                  }
                  double w[3];
                  for(int d=0; d < 3; d++)
                     if(cid == 1) w[d] = (I[d] == 0)? 0.0: ((r->rnd() > 0.5)? 1: -1) * sqrt(2.0*r->rnd()*T / I[d]);
                     else w[d] = (I2[d] == 0)? 0.0: ((r->rnd() > 0.5)? 1: -1) * sqrt(2.0*r->rnd()*T / I2[d]);
                  // xdr << "(" << ii[0] << "/" << ii[1] << "/" << ii[2] << "), j = " << j << ", d = " << d << ":\t";
                  if((format == FORMAT_BERNREUTHER) || (format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
                  {
                     xdr << id << " " << cid << "\t" << tr[0]/SIG_REF
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

