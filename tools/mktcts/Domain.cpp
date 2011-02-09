/*
 * GNU GPL version 2
 */

#include "Domain.h"
#include "Random.h"
#include <cmath>

#define BINS 2048
#define DT 0.002
#define PRECISION 5
#define TIME 20110209
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
   }

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
   for(unsigned l=0; l < 2; l++)
   {
      for(unsigned i=0; i < fl_units[0][l]; i++)
         for(unsigned j=0; j < fl_units[1][l]; j++)
            for(unsigned k=0; k < fl_units[2][l]; k++)
               for(unsigned d=0; d < 3; d++)
               {
                  if(l == 0) fill0[i][j][k][d] = true;
                  else fill1[i][j][k][d] = true;
               }
      slots[l] = 3 * fl_units[0][l] * fl_units[1][l] * fl_units[2][l];
      N[l] = slots[l];
   }
   bool tswap;
   double pswap;
   for(unsigned l=0; l < 2; l++)
   {
      for(unsigned m=0; m < PRECISION; m++)
      {
         tswap = (N[l] < N_id[l]);
         pswap = (N_id[l] - (double)N[l]) / ((tswap? slots[l]: 0) - (double)N[l]);
         // cout << "N = " << N[l] << ", N_id = " << N_id[l] << " => tswap = " << tswap << ", pswap = " << pswap << "\n";
         for(unsigned i=0; i < fl_units[0][l]; i++)
            for(unsigned j=0; j < fl_units[1][l]; j++)
               for(unsigned k=0; k < fl_units[2][l]; k++)
                  for(unsigned d=0; d < 3; d++)
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

   txt.precision(6);
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt << "mardynconfig\ntimestepLength\t" << DT << "\ncutoffRadius\t"
          << ((RDF > cutoff)? 1.031623*RDF: cutoff) << "\nLJCutoffRadius\t" << cutoff << "\n";
   }
   if(format == FORMAT_BRANCH)
   {
      txt << "tersoffCutoffRadius\t0.5\nphaseSpaceFile\t" << prefix << ".xdr\n";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      txt << "phaseSpaceFile\tOldStyle\t" << prefix << ".inp\n";
   }
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt << "datastructure\tLinkedCells\t1\n";

      txt << "output\tResultWriter\t1500\t" << prefix << "_1R\n";
   }
   if(format == FORMAT_BRANCH) txt << "resultOutputTimesteps\t1500\n";
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt << "output\tXyzWriter 60000\t" << prefix << "_1R\ninitCanonical\t10\ninitStatistics\t3003003\n";
      if(gradient)
      {
         txt << "profile\t1 " << BINS << " 1\nprofileRecordingTimesteps\t1\nprofileOutputTimesteps\t3000000\nprofiledComponent\t1\nprofileOutputPrefix\t" << prefix << "_1R\nAlignCentre\t100 0.01\n";
      }
      else
      {
         txt << "RDF\t" << RDF/(double)BINS << " " << BINS << "\nRDFOutputTimesteps\t3000000\nRDFOutputPrefix\t" << prefix << "_1R\n";
      }
      txt << "nomomentum\t1024\n";

      if(use_mu)
      {
         if(gradient)
         {
            txt << "chemicalPotential\t" << mu << " component 1\tcontrol 0 0 0 to " << box[0] << " " << 0.1*box[1] << " " << box[2] << "\tconduct " << (int)round(N[0] / 2500.0) << " tests every 2 steps\n";
            txt << "chemicalPotential\t" << mu << " component 1\tcontrol 0 " << 0.5*box[1] << " 0 to " << box[0] << " " << 0.6*box[1] << " " << box[2] << "\tconduct " << (int)round(N[1] / 2500.0) << " tests every 2 steps\n";
         }
         else
         {
            cout << "N[0] = " << N[0] << ", N[1] = " << N[1] << ", conducting " << round((double)(N[0]+N[1]) / 5000.0) << " test actions.\n";
            txt << "chemicalPotential\t" << mu << " component 1\tconduct " << (int)round((double)(N[0]+N[1]) / 5000.0) << " tests every 2 steps\n";
         }
         txt << "planckConstant\t" << sqrt(6.2831853 * T) << "\ninitGrandCanonical\t60000\n";
      }
   }

   xdr.precision(8);
   if(format == FORMAT_BRANCH)
   {
      xdr << "mardyn " << TIME << " tersoff\n"
          << "# mardyn input file, ls1 project\n"
          << "# generated by the mkTcTS tool\n# \n"
          << "t\t0\nL\t" << box[0] << " " << box[1]
          << " " << box[2] << "\nC\t1\t1 0 0 0 0\t0 0 0 1 1 1\t0 0 0 1e+10\nT\t" << T << "\nN\t"
          << N[0]+N[1] << "\nM\tICRVQD\n\n";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      xdr << "mardyn trunk " << TIME << "\n"
          << "t\t0\nL\t" << box[0] << " " << box[1]
          << " " << box[2] << "\nC\t1\t1 0 0 0 0\t0 0 0 1 1 1\t2.5 1\t0 0 0 1e+10\nT\t" << T << "\nN\t"
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
   for(unsigned l=0; l < 2; l++)
      for(unsigned i=0; i < fl_units[0][l]; i++)
         for(unsigned j=0; j < fl_units[1][l]; j++)
            for(unsigned k=0; k < fl_units[2][l]; k++)
               for(unsigned d=0; d < 3; d++)
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

                     if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
                     {
                        xdr << ID << " " << 1 << "\t" << q[0]
                            << " " << q[1] << " " << q[2]
                            << "\t" << v*cos(phi)*cos(omega) << " "
                            << v*cos(phi)*sin(omega) << " "
                            << v*sin(phi) << "\t1 0 0 0 0 0 0\n";
                     }
                     ID++;
                  }
                  else
                  {
                     if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
                     {
                        xdr << "\n";
                     }
                  }
               }

   xdr.close();
   txt.close();
   if(format == FORMAT_BUCHHOLZ)
   {
      buchholz.close();
   }
}

