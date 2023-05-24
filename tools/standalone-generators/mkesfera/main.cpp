#include "Domain.h"

#include <fstream>
#include <cstdio>
#include <math.h>
#include <string.h>


int main(int argc, char** argv) 
{
   const char* usage = "usage: mkesfera <prefix> [-e] -I <radius inside> -i <density inside> -O <radius outside> -o <density outside> [-R <cutoff>] [-r] [-S] -T <temperature> [-U] [-u]\n\n-e\tuse B-e-rnreuther format\n-r\tuse b-r-anch format (active by default)\n-S\tshift (active by default)\n-U\tunshift\n-u\tuse B-u-chholz format\n";
   if((argc < 12) || (argc > 16))
   {
      std::cout << "There are " << argc
           << " arguments where 12 to 16 should be given.\n\n";
      std::cout << usage;
      return 1;
   }

   bool do_shift = true;
   unsigned format = FORMAT_BRANCH;

   double cutoff = 2.5;
   double rho = 0.319;
   double rho2 = 0.319;
   double R = 7.0;
   double R2 = 14.0;
   double T = 1.0;

   char* prefix = argv[1];

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
         if(argv[i][j] == 'e') format = FORMAT_BERNREUTHER;
         else if(argv[i][j] == 'I')
         {
            i++;
            R = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'i')
         {
            i++;
            rho = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'O')
         {
            i++;
            R2 = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'o')
         {
            i++;
            rho2 = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'R')
         {
            i++;
            cutoff = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'r') format = FORMAT_BRANCH;
         else if(argv[i][j] == 'S') do_shift = true;
         else if(argv[i][j] == 'T')
         {
            i++;
            T = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'U') do_shift = false;
         else if(argv[i][j] == 'u') format = FORMAT_BUCHHOLZ;
         else
         {
            std::cout << "Invalid flag '-" << argv[i][j]
                 << "' was detected.\n\n" << usage;
            return 2; 
         }
      }
   }

   if(format == FORMAT_BERNREUTHER)
   {
      std::cout << "B-e-rnreuther format (flag -e) "
           << "is unavailable at present.\n\n" << usage;
      return 3;
   }

   Domain* dalet;
   dalet = new Domain(R, R2, rho, rho2);
   dalet->write(prefix, cutoff, T, do_shift, format);

   return 0;
}

