#include "Domain.h"

#include <fstream>
#include <cstdio>
#include <math.h>
#include <string.h>


int main(int argc, char** argv)
{
   const char* usage = "usage: animake <prefix> -c <density> [-e] [-f <fluid>] [-g <second component>] [-J <eta>] [-m <chemical potential>] -N <N_fluid> [-r] [-s <size unit [A]>] -T <temperature> [-u] [-W <energy and temperature unit [K]>] [-x <2nd comp. mole fract.>] [-Y <mass unit [u]>] [-y <y dim> [-z <z dim>]] [-5 <xi>]\n\n-e\tuse B-e-rnreuther format\n-f\tCH4 (default), Ar, C2H6, N2, CO2, EOX, JES, MER, TOL or VEG\n-r\tuse b-r-anch format\n-s\tgiven in units of Angstrom; default: 1 = 0.5291772 A\n-u\tuse B-u-chholz format (active by default)\n-W\tgiven in units of K; default value: 1 = 315774.5 K\n-Y\tgiven in units of g/mol; default value: 1 = 1000 g/mol\n\n";
   if((argc < 8) || (argc > 23))
   {
      std::cout << "There are " << argc
           << " arguments where 8 to 23 should be given.\n\n";
      std::cout << usage;
      return 1;
   }

   unsigned N;
   double x, dimx, dimy, dimz, T, rho, XIF, ETAF;
   double mu = 0.0;
   int fluid;
   int fluid2 = FLUID_NIL;
   int format = FORMAT_BUCHHOLZ;
   char* prefix = argv[1];
   bool muVT = false;

   double SIG_REF = 1.0;
   double EPS_REF = 1.0;
   double REFMASS = 1.0;

   bool in_N = false;
   bool in_rho = false;
   bool in_fluid = false;
   bool in_fluid2 = false;
   bool in_x = false;
   bool in_dimy = false;
   bool in_dimz = false;
   bool in_T = false;
   bool in_ETAF = false;
   bool in_XIF = false;

   for(int i=2; i < argc; i++)
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
         if(argv[i][j] == 'c')
         {
            in_rho = true;
            i++;
            rho = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'e') format = FORMAT_BERNREUTHER;
         else if(argv[i][j] == 'f')
         {
            in_fluid = true;
            i++;
            if(!strcmp(argv[i], "CH4")) fluid = FLUID_CH4;
            else if(!strcmp(argv[i], "Ar")) fluid = FLUID_AR;
            else if(!strcmp(argv[i], "C2H6")) fluid = FLUID_C2H6;
            else if(!strcmp(argv[i], "N2")) fluid = FLUID_N2;
            else if(!strcmp(argv[i], "CO2")) fluid = FLUID_CO2;
            else if(!strcmp(argv[i], "EOX")) fluid = FLUID_EOX;
            else if(!strcmp(argv[i], "JES")) fluid = FLUID_JES;
            else if(!strcmp(argv[i], "MER")) fluid = FLUID_MER;
            else if(!strcmp(argv[i], "TOL")) fluid = FLUID_TOL;
            else if(!strcmp(argv[i], "VEG")) fluid = FLUID_VEG;
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
            else if(!strcmp(argv[i], "CH4")) fluid2 = FLUID_CH4;
            else if(!strcmp(argv[i], "C2H6")) fluid2 = FLUID_C2H6;
            else if(!strcmp(argv[i], "N2")) fluid2 = FLUID_N2;
            else if(!strcmp(argv[i], "CO2")) fluid2 = FLUID_CO2;
            else if(!strcmp(argv[i], "EOX")) fluid2 = FLUID_EOX;
            else if(!strcmp(argv[i], "JES")) fluid2 = FLUID_JES;
            else if(!strcmp(argv[i], "MER")) fluid2 = FLUID_MER;
            else if(!strcmp(argv[i], "TOL")) fluid2 = FLUID_TOL;
            else if(!strcmp(argv[i], "VEG")) fluid2 = FLUID_VEG;
            else
            {
               std::cout << "(Secondary) fluid '" << argv[i]
                    << "' is not available.\n\n" << usage;
               return 99;
            }
            break;
         }
         else if(argv[i][j] == 'J')
         {
            in_ETAF = true;
            i++;
            ETAF = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'm')
         {
            muVT = true;
            i++;
            mu = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'N')
         {
            in_N = true;
            i++;
            N = atoi(argv[i]);
            break;
         }
         else if(argv[i][j] == 'r') format = FORMAT_BRANCH;
         else if(argv[i][j] == 's')
         {
            i++;
            SIG_REF = atof(argv[i]) / 0.5291772;
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
         else if(argv[i][j] == 'y')
         {
            in_dimy = true;
            i++;
            dimy = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'z')
         {
            in_dimz = true;
            i++;
            dimz = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == '5')
         {
            in_XIF = true;
            i++;
            XIF = atof(argv[i]);
            break;
         }
      }
   }

   if(!in_ETAF) ETAF = 1.0;
   if(!in_XIF) XIF = 1.0;
   if(!in_XIF && in_fluid2 && (((fluid == FLUID_MER) && (fluid2 == FLUID_TOL)) || ((fluid == FLUID_TOL) && (fluid2 == FLUID_MER))))
   {
      std::cerr << "Default binary parameter for the modified Berthelot mixing rule with CO2 (Merker) and toluene xi = 0.95.\n";
      XIF = 0.95;
   }
   if(!in_rho)
   {
      std::cout << "Fatal error: no fluid density was specified.\n\n";
      std::cout << usage;
      return 16;
   }
   if(in_dimz && !in_dimy)
   {
      std::cout << "Fatal error: z is specified while y is unknown. Please replace the -z option by -y.\n\n";
      std::cout << usage;
      return 17;
   }
   if(!in_fluid) fluid = FLUID_CH4;
   if(!in_N)
   {
      std::cout << "The essential input parameter "
           << "N (number of fluid molecules) is missing.\n\n" << usage;
      return 20;
   }
   if(!in_T)
   {
      std::cout << "Unspecified temperature: aborting.\n\n";
      std::cout << usage;
      return 21;
   }
   if(in_x && ((x < 0.0) || (x > 1.0)))
   {
      std::cout << "Invalid mole fraction x = " << x << ".\n\n" << usage;
      return 15;
   }
   if((fluid2 != FLUID_NIL) && !in_x)
   {
      std::cout << "Unspecified composition: aborting.\n\n";
      std::cout << usage;
      return 24;
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

   double V = (double)N/rho;
   if(!in_dimy) dimy = pow(V, 1.0/3.0);
   if(!in_dimz) dimz = sqrt(V/dimy);
   dimx = V / (dimy*dimz);

   double box[3];
   box[0] = dimx * SIG_REF;
   box[1] = dimy * SIG_REF;
   box[2] = dimz * SIG_REF;
   T *= EPS_REF;
   mu *= EPS_REF;

   Domain* delta = new Domain(
      fluid, fluid2, box, SIG_REF, EPS_REF, REFMASS, muVT, N, T, ETAF, XIF
   );
   delta->write(
      prefix, format, mu, x
   );

   return 0;
}

