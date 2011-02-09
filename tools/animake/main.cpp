/*
 * GNU GPL version 2
 */

#include "Domain.h"

#include <fstream>
#include <cstdio>
#include <math.h>
#include <string.h>

using namespace std;

int main(int argc, char** argv) 
{
   const char* usage = "usage: animake <prefix> -c <density> [-e] [-f <fluid>] [-m <chemical potential>] -N <N_fluid> [-r] [-s <size unit [A]>] -T <temperature> [-u] [-W <energy and temperature unit [K]>] [-Y <mass unit [u]>] [-y <y dim> [-z <z dim>]]\n\n-e\tuse B-e-rnreuther format\n-f\tCH4 (default), Ar, C2H6, N2, CO2, or EOX\n-r\tuse b-r-anch format (active by default)\n-s\tgiven in units of Angstrom; default: 1 = 0.5291772 A\n-u\tuse B-u-chholz format\n-W\tgiven in units of K; default value: 1 = 315774.5 K\n-Y\tgiven in units of g/mol; default value: 1 = 1000 g/mol\n\n";
   if((argc < 8) || (argc > 23))
   {
      cout << "There are " << argc
           << " arguments where 8 to 23 should be given.\n\n";
      cout << usage;
      return 1;
   }

   unsigned N;
   double x, y, z, T, rho;
   double mu = 0.0;
   int fluid;
   int format = FORMAT_BRANCH;
   char* prefix = argv[1];
   bool muVT = false;

   double SIG_REF = 1.0;
   double EPS_REF = 1.0;
   double REFMASS = 1.0;

   bool in_N = false;
   bool in_rho = false;
   bool in_fluid = false;
   bool in_y = false;
   bool in_z = false;
   bool in_T = false;

   for(int i=2; i < argc; i++)
   {
      if(*argv[i] != '-')
      {
         cout << "Flag expected where '" << argv[i]
              << "' was given.\n\n";
         cout << usage;
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
            else
            {
               cout << "Fluid '" << argv[i] 
                    << "' is not available.\n\n" << usage;
               return 9;
            }
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
         else if(argv[i][j] == 'Y')
         {
            i++;
            REFMASS = atof(argv[i]) / 1000.0;
            break;
         }
         else if(argv[i][j] == 'y')
         {
            in_y = true;
            i++;
            y = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'z')
         {
            in_z = true;
            i++;
            z = atof(argv[i]);
            break;
         }
      }
   }

   if(!in_rho)
   {
      cout << "Fatal error: no fluid density was specified.\n\n";
      cout << usage;
      return 16;
   }
   if(in_z && !in_y)
   {
      cout << "Fatal error: z is specified while y is unknown. Please replace the -z option by -y.\n\n";
      cout << usage;
      return 17;
   }
   if(format == FORMAT_BERNREUTHER)
   {
      cout << "B-e-rnreuther format (flag -e) "
           << "is unavailable at present.\n\n" << usage;
      return 18;
   }
   if(!in_fluid) fluid = FLUID_CH4;
   if(!in_N)
   {
      cout << "The essential input parameter "
           << "N (number of fluid molecules) is missing.\n\n" << usage;
      return 20;
   }
   if(!in_T)
   {
      cout << "Unspecified temperature: aborting.\n\n";
      cout << usage;
      return 21;
   }

   double V = (double)N/rho;
   if(!in_y) y = pow(V, 1.0/3.0);
   if(!in_z) z = sqrt(V/y);
   x = V / (y*z);

   double box[3];
   box[0] = x * SIG_REF;
   box[1] = y * SIG_REF;
   box[2] = z * SIG_REF;
   T *= EPS_REF;
   mu *= EPS_REF;

   Domain* delta = new Domain(
      fluid, box, SIG_REF, EPS_REF, REFMASS, muVT, N, T
   );
   delta->write(
      prefix, format, mu
   );

   return 0;
}

