#include "Domain.h"

#include <fstream>
#include <cstdio>
#include <math.h>
#include <string.h>


int main(int argc, char** argv) 
{
   const char* usage = "usage: mkTcTS <prefix> -c <density> [-a] [-d <second density>] [-e] [-H <pressure1> <pressure2> <mu_low> <mu_high>] [-h <height>] [-m <chemical potential>] -N <particles> [-p <pair correlation cutoff>] [-R <cutoff>] [-r] [-S] -T <temperature> [-U] [-u]\n\n-a\tcompute velocity autocorrelation\n-e\tuse B-e-rnreuther format\n-P\tcompute pseudomomenta\n-r\tuse b-r-anch format\n-S\tshift (active by default)\n-U\tunshift\n-u\tuse B-u-chholz format (active by default)\n";
   if((argc < 8) || (argc > 28))
   {
      std::cout << "There are " << argc
           << " arguments where 8 to 28 should be given.\n\n";
      std::cout << usage;
      return 1;
   }

   bool compute_autocorr = false;
   bool do_shift = true;
   bool in_h = false;
   bool use_mu = false;
   bool gradient = false;
   unsigned format = FORMAT_BUCHHOLZ;

   double cutoff = 2.5;
   double RDF = 12.0;
   double h = 0.0;
   double mu = 0.0;
   unsigned N = 131072;
   double rho = 0.319;
   double rho2 = 0.319;
   double T = 1.0779;
   
   bool use_hato = false;
   double p1 = 0.0;
   double p2 = 0.0;
   double mu_low = 0.0;
   double mu_high = 0.0;

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
         if(argv[i][j] == 'a') compute_autocorr = true;
         else if(argv[i][j] == 'c')
         {
            i++;
            rho = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'd')
         {
            gradient = true;
            i++;
            rho2 = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'e') format = FORMAT_BERNREUTHER;
         else if(argv[i][j] == 'H')
         {
            use_hato = true;
            i++;
            p1 = atof(argv[i]);
            i++;
            p2 = atof(argv[i]);
            i++;
            mu_low = atof(argv[i]);
            i++;
            mu_high = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'h')
         {
            in_h = true;
            i++;
            h = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'm')
         {
            use_mu = true;
            i++;
            mu = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'N')
         {
            i++;
            N = atoi(argv[i]);
            break;
         }
         else if(argv[i][j] == 'p')
         {
            i++;
            RDF = atof(argv[i]);
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

   if(in_h && !gradient)
   {
      std::cout << "The box dimension can only be specified for "
           << "systems with a density gradient.\n\n" << usage;
      return 5;
   }

   if(format == FORMAT_BERNREUTHER)
   {
      std::cout << "B-e-rnreuther format (flag -e) "
           << "is unavailable at present.\n\n" << usage;
      return 3;
   }

   if(!in_h) h = pow((double)N/rho, 1.0/3.0);

   Domain* dalet;
   if(gradient) dalet = new Domain(h, N, rho, rho2);
   else dalet = new Domain(N, rho, RDF);
   if(use_hato) dalet->hato(p1, p2, mu_low, mu_high);
   dalet->write(prefix, cutoff, mu, T, do_shift, use_mu, compute_autocorr, format);

   return 0;
}

