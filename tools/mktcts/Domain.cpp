/*
 * GNU GPL version 2
 */

#include "Domain.h"
#include "Random.h"
#include <cmath>

#define BINS 1024
#define DT 0.0015
#define PRECISION 6
#define TIME 20100210
#define VARFRACTION 0.125

Domain::Domain(double t_h, unsigned t_N, double t_rho, double t_rho2)
{
   this->gradient = true;
   this->rho = t_rho;
   this->rho2 = t_rho2;
   this->RDF = 0.0;

   double V = 2.0 * t_N / (rho + rho2);
   this->box[0] = sqrt(V / t_h);
   this->box[1] = t_h;
   this->box[2] = sqrt(V / t_h);
}

Domain::Domain(unsigned t_N, double t_rho, double t_RDF)
{
   this->gradient = false;
   this->rho = t_rho;
   this->rho2 = t_rho;
   this->RDF = t_RDF;

   this->box[0] = pow((double)t_N/rho, 1.0/3.0);
   this->box[1] = this->box[0];
   this->box[2] = this->box[0];
}

void Domain::write(char* prefix, double cutoff, double mu, double T, bool do_shift, bool use_mu, int format)
{
   ofstream xdr, txt;
   stringstream strstrm, txtstrstrm;
   strstrm << prefix << ".xdr";
   xdr.open(strstrm.str().c_str(), ios::trunc);
   txtstrstrm << prefix << "_1R.txt";
   txt.open(txtstrstrm.str().c_str(), ios::trunc);

   unsigned fl_units[3][2];
   double fl_unit[3][2];
   double N_id[2];
   for(int i=0; i < 2; i++)
   {
      N_id[i] = box[0]*(0.5*box[1])*box[2] * ((i == 0)? rho: rho2);
      double N_boxes = N_id[i] / 3.0;
      fl_units[1][i] = round(
                          pow(
                             (N_boxes * (0.5*box[1]) * (0.5*box[1]))
                                      / (this->box[0] * this->box[2]), 1.0/3.0
                          )
                       );
      if(fl_units[1][i] == 0) fl_units[1][i] = 1;
      double bxbz_id = N_boxes / fl_units[1][i];
      fl_units[0][i] = round(sqrt(this->box[0] * bxbz_id / this->box[2]));
      if(fl_units[0][i] == 0) fl_units[0][i] = 1;
      fl_units[2][i] = ceil(bxbz_id / fl_units[0][i]);
      for(int d=0; d < 3; d++) fl_unit[d][i] = ((d == 1)? 0.5: 1.0) * box[d] / (double)fl_units[d][i];
      cout << "Elementary cell " << i << ": " << fl_unit[0][i] << " x " << fl_unit[1][i] << " x " << fl_unit[2][i] << ".\n";
   }

   Random* r = new Random();
   r->init(
      (int)(10000.0*box[0]) - (int)(3162.3*cutoff)
        + (int)(1000.0*T) - (int)(316.23*mu)
        + (int)(100.0*box[1])
   );

   bool fill0[fl_units[0][0]][fl_units[1][0]][fl_units[2][0]][3];
   bool fill1[fl_units[0][1]][fl_units[1][1]][fl_units[2][1]][3];
   unsigned N[2];
   unsigned slots[2];
   for(int l=0; l < 2; l++)
   {
      for(int i=0; i < fl_units[0][l]; i++)
         for(int j=0; j < fl_units[1][l]; j++)
            for(int k=0; k < fl_units[2][l]; k++)
               for(int d=0; d < 3; d++)
               {
                  if(l == 0) fill0[i][j][k][d] = true;
                  else fill1[i][j][k][d] = true;
               }
      slots[l] = 3 * fl_units[0][l] * fl_units[1][l] * fl_units[2][l];
      N[l] = slots[l];
   }
   bool tswap;
   double pswap;
   for(int l=0; l < 2; l++)
   {
      for(int m=0; m < PRECISION; m++)
      {
         tswap = (N[l] < N_id[l]);
         pswap = (N_id[l] - (double)N[l]) / ((tswap? slots[l]: 0) - (double)N[l]);
         // cout << "N = " << N[l] << ", N_id = " << N_id[l] << " => tswap = " << tswap << ", pswap = " << pswap << "\n";
         for(int i=0; i < fl_units[0][l]; i++)
            for(int j=0; j < fl_units[1][l]; j++)
               for(int k=0; k < fl_units[2][l]; k++)
                  for(int d=0; d < 3; d++)
                     if(pswap >= r->rnd())
                     {
                        if(((l == 0) && fill0[i][j][k][d]) || ((l == 1) && fill1[i][j][k][d])) N[l] --;
                        if(l == 0) fill0[i][j][k][d] = tswap;
                        else fill1[i][j][k][d] = tswap;
                        if(tswap) N[l] ++;
                     }
      }
      cout << "Filling " << N[l] << " of 3*"
           << fl_units[0][l] << "*" << fl_units[1][l] << "*" << fl_units[2][l]
           << " = " << slots[l] << " slots (ideally " << N_id[l] << ").\n";
   }

   if(format == FORMAT_BRANCH)
   {
      txt.precision(6);
      txt << "mardynconfig\n# \ntimestepLength\t" << DT << "\ncutoffRadius\t" << ((RDF > cutoff)? 1.031623*RDF: cutoff)  << "\nLJCutoffRadius\t" << cutoff << "\ntersoffCutoffRadius\t0.5\n";

      txt << "phaseSpaceFile\t" << prefix << ".xdr\n";

      txt << "datastructure\tLinkedCells 1\n";

      txt << "output\tResultWriter 200\t" << prefix << "_1R\nresultOutputTimesteps\t200\noutput\tXyzWriter 50000\t" << prefix << "_1R\ninitCanonical\t10\ninitStatistics\t120000\n";
      if(gradient)
      {
         txt << "profile\t1 " << BINS << " 1\nprofileRecordingTimesteps\t1\nprofileOutputTimesteps\t40000\nprofiledComponent\t1\nprofileOutputPrefix\t" << prefix << "_1R\n";
      }
      else
      {
         txt << "RDF\t" << RDF/(double)BINS << " " << BINS << "\nRDFOutputTimesteps\t40000\nRDFOutputPrefix\t" << prefix << "_1R\n";
      }

      if(use_mu)
      {
         if(gradient)
         {
            txt << "chemicalPotential\t" << mu << " component 1\tcontrol 0 0 0 to " << box[0] << " " << 0.1*box[1] << " " << box[2] << "\tconduct " << (int)round(N[2] / 5000.0) << " tests every 1 steps\n";
            txt << "chemicalPotential\t" << mu << " component 1\tcontrol 0 " << 0.5*box[1] << " 0 to " << box[0] << " " << 0.6*box[1] << " " << box[2] << "\tconduct " << (int)round(N[1] / 5000.0) << " tests every 1 steps\n";
         }
         else
         {
            txt << "chemicalPotential\t" << mu << " component 1\tconduct " << (int)round((N[1]+N[2]) / 1000.0) << " tests every 1 steps\n";
         }
         txt << "planckConstant\t" << sqrt(6.2831853 * T) << "\ninitGrandCanonical\t60000\n";
      }
   }

   if(format == FORMAT_BRANCH)
   {
      xdr.precision(8);
      xdr << "mardyn " << TIME << " tersoff\n"
          << "# mardyn input file, ls1 project\n"
          << "# generated by the mkTcTS tool\n# \n";

      xdr << "t\t0\nL\t" << box[0] << " " << box[1]
          << " " << box[2] << "\nC\t1\t1 0 0 0 0\t0 0 0 1 1 1\t0 0 0 1e+10\nT\t" << T << "\nN\t"
          << N[0]+N[1] << "\nM\tICRVQD\n\n";
   }

   double v = sqrt(3.0 * T);
   double loffset[3][2];
   loffset[0][0] = 0.1; loffset[1][0] = 0.3; loffset[2][0] = 0.1;
   loffset[0][1] = 0.1; loffset[1][1] = 0.8; loffset[2][1] = 0.1;
   double goffset[3][3];
   goffset[0][0] = 0.0; goffset[1][0] = 0.5; goffset[2][0] = 0.5;
   goffset[0][1] = 0.5; goffset[1][1] = 0.0; goffset[2][1] = 0.5;
   goffset[0][2] = 0.5; goffset[1][2] = 0.5; goffset[2][2] = 0.0;

   unsigned ID = 1;
   for(int l=0; l < 2; l++)
      for(int i=0; i < fl_units[0][l]; i++)
         for(int j=0; j < fl_units[1][l]; j++)
            for(int k=0; k < fl_units[2][l]; k++)
               for(int d=0; d < 3; d++)
               {
                  if(((l == 0) && fill0[i][j][k][d]) || ((l == 1) && fill1[i][j][k][d]))
                  {
                     double q[3];
                     q[0] = i * fl_unit[0][l];
                     q[1] = j * fl_unit[1][l];
                     q[2] = k * fl_unit[2][l];
                     for(int m=0; m < 3; m++)
                     {
                        q[m] += box[m]*loffset[m][l] + fl_unit[m][l]*goffset[m][d];
                        q[m] += VARFRACTION * fl_unit[m][l] * (r->rnd() - 0.5);
                        if(q[m] > box[m]) q[m] -= box[m];
                     }
                     double phi = 6.283185 * r->rnd();
                     double omega = 6.283185 * r->rnd();

                     xdr << ID << " " << 1 << "\t" << q[0]
                         << " " << q[1] << " " << q[2]
                         << "\t" << v*cos(phi)*cos(omega) << " "
                         << v*cos(phi)*sin(omega) << " "
                         << v*sin(phi) << "\t1 0 0 0 0 0 0\n";

                     ID++;
                  }
                  else
                  {
                     if(format == FORMAT_BRANCH)
                     {
                        xdr << "\n";
                     }
                  }
               }

   xdr.close();
   txt.close();
}

