#include "Domain.h"
#include "Graphit.h"

#include <fstream>
#include <cstdio>
#include <math.h>
#include <string.h>


#define DEFAULT_TAU 1.0e+10

#define DEFAULT_CC_BONDLENGTH 2.6853
#define ORIGINAL_CC_BONDLENGTH 2.7609

// interlayer spacing graphite: 3.35 A
// interlayer spacing MWCNT: 3.40 A
//    (http://www.wag.caltech.edu/foresight/foresight_2.html)

int main(int argc, char** argv) 
{
   const char* usage = "usage: mkcp (-C|-P|-0) <prefix> [-a <initial acceleration>] [-A <C-C bond length>] -c <density> -d <layers> [-e] [-E] [-f <fluid>] [-g <second component>] -h <height> [-H <eta>] [-I <eta2>] [-J <fluid eta>] [-k] [-l] [-L] [-m <chemical potential>] [-M <CNT with m/n>] -N <N_fluid> [-p <polarity coefficient>] [-r] [-s <size unit [A]>] [-S] [-t <controller time parameter>] -T <temperature> [-u] -U <velocity> [-v <volume fraction without acceleration>] [-V <volume fraction without wall] [-w] [-W <energy and temperature unit [K]>] [-x <2nd comp. mole fract.>] [-Y <mass unit [u]>] [-3 <xi>] [-4 <xi2>] [-5 <fluid xi>] [-8]\n\n-A\treduced C-C bond length; default: 2.6853 a0 = 0.1421 nm, original Tersoff: 2.7609 a0\n-C\tCouette flow (flag followed by output prefix)\n-e\tuse B-e-rnreuther format\n-E\tgenerate an empty channel\n-f\tAr (default), Ave, CH4, C2H6, N2, CO2, H2O, CH3OH, or C6H14\n-H, -I\tdefault: eta according to Wang et al.\n-J\tdefault: eta = 1\n-k\tonly harmonic potentials active within the wall (default for polar walls)\n-l\tLJ interaction within the wall as well (default for unpolar walls)\n-L\tuse Lennard-Jones units instead of atomic units (cf. atomic_units.txt)\n-M\tgenerate a carbon nanotube (only for Poiseuille flow)\n-p\tdefault polarity coefficient: 0 (unpolar walls)\n-r\tuse b-r-anch format (active by default)\n-s\tgiven in units of Angstrom; default: 1 = 0.5291772 A 0\n-S\tsymmetric system (with two identical fluid components)\n-t\tdefault: tau extremely large (about 30 us)\n-u\tuse B-u-chholz format\n-w\tWidom\n-W\tgiven in units of K; default value: 1 = 315774.5 K\n-x\tdefault value: 0 (i.e. pure first comp. fluid)\n-Y\tgiven in units of g/mol; default value: 1 = 1000 g/mol\n-0\tstatic scenario without flow (followed by output prefix)\n-3, -4\tdefault: xi according to Wang et al.\n-5\tdefault: xi = 1\n-8\toriginal Tersoff potential as published in the 80s\n";
   if((argc < 15) || (argc > 69))
   {
      std::cout << "There are " << argc
           << " arguments where 15 to 61 should be given.\n\n";
      std::cout << usage;
      return 1;
   }

   unsigned d, N;
   double h, T, rho, U, XI, ETA, XI2, ETA2, XIF, ETAF, TAU, bondlength, a, m_per_n;
   double mu = 0.0;
   double x = 0.0;
   int fluid, flow;
   int fluid2 = FLUID_NIL;
   int format = FORMAT_BRANCH;
   char* prefix = (char*)0;

   double SIG_REF = 1.0;
   double EPS_REF = 1.0;
   double REFMASS = 1.0;

   double wo_acceleration = 0.0;
   double wo_wall = 0.0;
   double polarity = 0.0;

   bool in_a = false;
   bool in_bondlength = false;
   bool empty = false;
   bool in_ETA = false;
   bool in_ETA2 = false;
   bool in_ETAF = false;
   bool LJunits = false;
   bool in_TAU = false;
   bool in_N = false;
   bool in_XI = false;
   bool in_XI2 = false;
   bool in_XIF = false;
   bool in_rho = false;
   bool in_d = false;
   bool in_fluid = false;
   bool in_fluid2 = false;
   bool in_h = false;
   bool in_T = false;
   bool in_U = false;
   bool in_WLJ = false;
   bool in_x = false;
   bool original = false;
   bool muVT = false;
   bool nanotube = false;
   bool WLJ = true;
   bool symmetric = false;
   bool widom = false;

   for(int i=1; i < argc; i++)
   {
      if(*argv[i] != '-')
      {
         std::cout << "Flag expected where '" << argv[i]
              << "' was given.\n\n";
         std::cout << usage;
         return 2;
      }
      for(int j=1; argv[i][j]; j++)
      {
         // std::cout << argv[i][j] << "\n";
         if(argv[i][j] == 'a')
         {
            in_a = true;
            i++;
            a = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'A')
         {
            in_bondlength = true;
            i++;
            bondlength = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'c')
         {
            in_rho = true;
            i++;
            rho = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'C')
         {
            flow = FLOW_COUETTE;
            i++;
            prefix = argv[i];
            break;
         }
         else if(argv[i][j] == 'd')
         {
            in_d = true;
            i++;
            d = atoi(argv[i]);
            break;
         }
         else if(argv[i][j] == 'e') format = FORMAT_BERNREUTHER;
         else if(argv[i][j] == 'E') empty = true;
         else if(argv[i][j] == 'f')
         {
            in_fluid = true;
            i++;
            if(!strcmp(argv[i], "Ar")) fluid = FLUID_AR;
            else if(!strcmp(argv[i], "Ave")) fluid = FLUID_AVE;
            else if(!strcmp(argv[i], "CH4")) fluid = FLUID_CH4;
            else if(!strcmp(argv[i], "C2H6")) fluid = FLUID_C2H6;
            else if(!strcmp(argv[i], "N2")) fluid = FLUID_N2;
            else if(!strcmp(argv[i], "CO2")) fluid = FLUID_CO2;
            else if(!strcmp(argv[i], "H2O")) fluid = FLUID_H2O;
            else if(!strcmp(argv[i], "CH3OH")) fluid = FLUID_CH3OH;
            else if(!strcmp(argv[i], "C6H14")) fluid = FLUID_C6H14;
            else
            {
               std::cout << "Fluid '" << argv[i] 
                    << "' is not available.\n\n" << usage;
               return 9;
            }
            break;
         }
         else if(argv[i][j] == 'g')
         {
            in_fluid2 = true;
            i++;
            if(!strcmp(argv[i], "Ar")) fluid2 = FLUID_AR;
            else if(!strcmp(argv[i], "Ave")) fluid2 = FLUID_AVE;
            else if(!strcmp(argv[i], "CH4")) fluid2 = FLUID_CH4;
            else if(!strcmp(argv[i], "C2H6")) fluid2 = FLUID_C2H6;
            else if(!strcmp(argv[i], "N2")) fluid2 = FLUID_N2;
            else if(!strcmp(argv[i], "CO2")) fluid2 = FLUID_CO2;
            else if(!strcmp(argv[i], "H2O")) fluid2 = FLUID_H2O;
            else if(!strcmp(argv[i], "CH3OH")) fluid2 = FLUID_CH3OH;
            else if(!strcmp(argv[i], "C6H14")) fluid2 = FLUID_C6H14;
            else
            {
               std::cout << "(Secondary) fluid '" << argv[i] 
                    << "' is not available.\n\n" << usage;
               return 99;
            }
            break;
         }
         else if(argv[i][j] == 'h')
         {
            in_h = true;
            i++;
            h = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'H')
         {
            in_ETA = true;
            i++;
            ETA = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'I')
         {
            in_ETA2 = true;
            i++;
            ETA2 = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'J')
         {
            in_ETAF = true;
            i++;
            ETAF = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'k')
         {
            in_WLJ = true;
            WLJ = false;
         }
         else if(argv[i][j] == 'l')
         {
            in_WLJ = true;
            WLJ = true;
         }
         else if(argv[i][j] == 'L') LJunits = true;
         else if(argv[i][j] == 'm')
         {
            muVT = true;
            i++;
            mu = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'M')
         {
            nanotube = true;
            i++;
            m_per_n = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'N')
         {
            in_N = true;
            i++;
            N = atoi(argv[i]);
            break;
         }
         else if(argv[i][j] == 'p')
         {
            i++;
            polarity = atof(argv[i]);
            if(!in_WLJ)
            {
               WLJ = (polarity == 0.0);
            }
            break;
         }
         else if(argv[i][j] == 'P')
         {
            flow = FLOW_POISEUILLE;
            i++;
            prefix = argv[i];
            break;
         }
         else if(argv[i][j] == 'r') format = FORMAT_BRANCH;
         else if(argv[i][j] == 'S') symmetric = true;
         else if(argv[i][j] == 's')
         {
            i++;
            SIG_REF = atof(argv[i]) / 0.5291772;
            break;
         }
         else if(argv[i][j] == 't')
         {
            in_TAU = true;
            i++;
            TAU = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'T')
         {
            in_T = true;
            i++;
            T = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'u') format = FORMAT_BUCHHOLZ;
         else if(argv[i][j] == 'U')
         {
            in_U = true;
            i++;
            U = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'v')
         {
            i++;
            wo_acceleration = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'V')
         {
            i++;
            wo_wall = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'w') widom = true;
         else if(argv[i][j] == 'W')
         {
            i++;
            EPS_REF = atof(argv[i]) / 315774.5;
            break;
         }
         else if(argv[i][j] == 'x')
         {
            in_x = true;
            i++;
            x = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'Y')
         {
            i++;
            REFMASS = atof(argv[i]) / 1000.0;
            break;
         }
         else if(argv[i][j] == '0')
         {
            flow = FLOW_NONE;
            i++;
            prefix = argv[i];
            break;
         }
         else if(argv[i][j] == '3')
         {
            in_XI = true;
            i++;
            XI = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == '4')
         {
            in_XI2 = true;
            i++;
            XI2 = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == '5')
         {
            in_XIF = true;
            i++;
            XIF = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == '8') original = true;
         else
         {
            std::cout << "Invalid flag '-" << argv[i][j]
                 << "' was detected.\n\n" << usage;
            return 2; 
         }
      }
   }

   if(wo_wall < 0.003) wo_wall = 0.0;
   else if(wo_wall > 0.997) wo_wall = 1.0;

   if(nanotube)
   {
      if(flow == FLOW_COUETTE)
      {
         std::cout << "Couette flow is not well-defined for nanotubes.\n\n"
              << usage;
         return 18;
      }

      std::cout << "Carbon nanotubes are not yet implemented.\n\n"
           << usage;
      return 11;
   }
   if(WLJ && (polarity != 0.0))
   {
      std::cout << "Severe error: Non-zero polarity is incompatible with the LJ interlayer interaction.\n\n"
           << usage;
      return 181;
   }
   if(!prefix)
   {
      std::cout << "You have to specify an output prefix "
           << "via the -C, -P, or -0 option.\n\n" << usage;
      return 5;
   }
   if(!in_rho)
   {
      std::cout << "Fatal error: no fluid density was specified.\n\n";
      std::cout << usage;
      return 6;
   }
   if(!in_d)
   {
      std::cout << "Wall thickness unknown. Reconsider parameter set.\n\n";
      std::cout << usage;
      return 7;
   }
   if(format == FORMAT_BERNREUTHER)
   {
      std::cout << "B-e-rnreuther format (flag -e) "
           << "is unavailable at present.\n\n" << usage;
      return 14;
   }

   if(fluid == fluid2)
   {
      symmetric = true;
      fluid2 = FLUID_NIL;
   }
   if((x < 0.0) || (x > 1.0))
   {
      std::cout << "Invalid mole fraction x = " << x << ".\n\n" << usage;
      return 15;
   }
   if((in_fluid2 == false) || (fluid2 == FLUID_NIL))
   {
      x = 0.0;
   }
   if(x == 1.0)
   {
      fluid = fluid2;
      x = 0.0;
   }
   if(x == 0.0)
   {
      in_fluid2 = false;
      fluid2 = FLUID_NIL;
   }

   if(symmetric && (fluid2 != FLUID_NIL))
   {
      std::cout << "The symmetric scenario cannot contain an (actual) fluid mixture.\n\n" << usage;
      return 16;
   }
   if(symmetric) x = 0.5;

   if(!in_fluid) fluid = FLUID_AR;
   double SIG_FLUID, EPS_FLUID, FLUIDMASS;
   double SIG_FLUID2, EPS_FLUID2, FLUIDMASS2;
   if(fluid == FLUID_CH4)
   {
      SIG_FLUID = SIG_CH4;
      EPS_FLUID = EPS_CH4;
      FLUIDMASS = CH4MASS;
   }
   else if(fluid == FLUID_C2H6)
   {
      SIG_FLUID = SIG_C2H6;
      EPS_FLUID = EPS_C2H6;
      FLUIDMASS = C2H6MASS;
   }
   else if(fluid == FLUID_N2)
   {
      SIG_FLUID = SIG_N2;
      EPS_FLUID = EPS_N2;
      FLUIDMASS = N2MASS;
   }
   else if(fluid == FLUID_CO2)
   {
      SIG_FLUID = SIG_CO2;
      EPS_FLUID = EPS_CO2;
      FLUIDMASS = CO2MASS;
   }
   else if(fluid == FLUID_AVE) // Avendaño mode
   {
      SIG_FLUID = SIG_AVE;
      EPS_FLUID = EPS_AVE;
      FLUIDMASS = AVEMASS;
   }
   else if(fluid == FLUID_AR)
   {
      SIG_FLUID = SIG_AR;
      EPS_FLUID = EPS_AR;
      FLUIDMASS = ARMASS;
   }
   else // effectively sets ETA = XI = 1 in all other cases
   {
      SIG_FLUID = SIG_WANG;
      EPS_FLUID = EPS_WANG;
      if(fluid == FLUID_H2O) FLUIDMASS = OH2OMASS + 2.0*HH2OMASS;
      else if(fluid == FLUID_CH3OH) FLUIDMASS = CCH3OHMASS + OCH3OHMASS + HCH3OHMASS;
      else FLUIDMASS = 2.0*FC6H14MASS + 4.0*MC6H14MASS;
   }
   if(fluid2 == FLUID_CH4)
   {
      SIG_FLUID2 = SIG_CH4;
      EPS_FLUID2 = EPS_CH4;
      FLUIDMASS2 = CH4MASS;
   }
   else if(fluid2 == FLUID_C2H6)
   {
      SIG_FLUID2 = SIG_C2H6;
      EPS_FLUID2 = EPS_C2H6;
      FLUIDMASS2 = C2H6MASS;
   }
   else if(fluid2 == FLUID_N2)
   {
      SIG_FLUID2 = SIG_N2;
      EPS_FLUID2 = EPS_N2;
      FLUIDMASS2 = N2MASS;
   }
   else if(fluid2 == FLUID_CO2)
   {
      SIG_FLUID2 = SIG_CO2;
      EPS_FLUID2 = EPS_CO2;
      FLUIDMASS2 = CO2MASS;
   }
   else if(fluid2 == FLUID_AVE) // Avendaño mode
   {
      SIG_FLUID2 = SIG_AVE;
      EPS_FLUID2 = EPS_AVE;
      FLUIDMASS2 = AVEMASS;
   }
   else if(fluid2 == FLUID_AR)
   {
      SIG_FLUID2 = SIG_AR;
      EPS_FLUID2 = EPS_AR;
      FLUIDMASS2 = ARMASS;
   }
   else // effectively sets ETA2 = XI2 = 1 in all other cases
   {
      SIG_FLUID2 = SIG_WANG;
      EPS_FLUID2 = EPS_WANG;
      if(fluid2 == FLUID_H2O) FLUIDMASS2 = OH2OMASS + 2.0*HH2OMASS;
      else if(fluid2 == FLUID_CH3OH) FLUIDMASS2 = CCH3OHMASS + OCH3OHMASS + HCH3OHMASS;
      else FLUIDMASS2 = 2.0*FC6H14MASS + 4.0*MC6H14MASS;
   }
   if(!in_N)
   {
      std::cout << "Missing essential input parameter "
           << "N (number of fluid molecules).\n\n" << usage;
      return 12;
   }
   if(!in_T)
   {
      std::cout << "No temperature specified: aborting.\n\n";
      std::cout << usage;
      return 13;
   }
   if(!in_U) U = 0.0;
   else if(FLOW_NONE && (U != 0.0))
   {
      std::cout << "Specified flow velocity U = " << U 
           << " is incompatible with -0 option (static scenario).\n\n"
           << usage;
      return 15;
   }
   if(muVT && symmetric)
   {
      std::cout << "The combination between the uVT (-m) and symmetry (-S) options remains unsupported.\n\n"
           << usage;
      return 115;
   }

   if(!in_h) h = pow((double)N/rho, 1.0/3.0);
   if(!in_ETA) ETA = 0.5*(1.0 + SIG_WANG/SIG_FLUID);
   if(!in_XI) XI = sqrt(EPS_WANG / EPS_FLUID);
   if(!in_ETA2) ETA2 = (SIG_WANG + SIG_FLUID2) / (SIG_FLUID + SIG_FLUID2);
   if(!in_XI2) XI2 = sqrt(EPS_WANG / EPS_FLUID);
   if(!in_ETAF) ETAF = 1.0;
   if(!in_XIF) XIF = 1.0;

   SIG_REF = LJunits? SIG_FLUID: 1.0;
   EPS_REF = LJunits? EPS_FLUID: 1.0;
   REFMASS = LJunits? FLUIDMASS: 1.0;
   if(LJunits)
   {
      SIG_REF = SIG_FLUID;
      EPS_REF = EPS_FLUID;
      REFMASS = FLUIDMASS;
   }
   double REFTIME = SIG_REF * sqrt(REFMASS / EPS_REF);

   h *= SIG_REF;
   T *= EPS_REF;
   rho /= (SIG_REF * SIG_REF * SIG_REF);
   U *= SIG_REF / REFTIME;
   TAU *= REFTIME;
   bondlength *= SIG_REF;
   a *= SIG_REF / (REFTIME * REFTIME);
   mu *= EPS_REF;

   if(!in_bondlength)
      bondlength = original? ORIGINAL_CC_BONDLENGTH
                           : DEFAULT_CC_BONDLENGTH;
   if(!in_TAU) TAU = DEFAULT_TAU;
   if(!in_a) a = 0.01 * U / TAU;

   Domain* dalet = new Domain(
      flow, bondlength, rho, d, fluid, fluid2, h, ETA, ETA2, ETAF, SIG_REF, EPS_REF,
      REFMASS, muVT, nanotube, m_per_n, N, T, XI, XI2, XIF, wo_wall
   );
   if(nanotube)
   {
      /*
       * not yet implemented
       */
   }
   else
   {
      dalet->write(
         prefix, a, empty, format, mu, TAU, U, original,
         wo_acceleration, polarity, WLJ, symmetric, widom, x
      );
   }

   return 0;
}

