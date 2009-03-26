#include <iomanip>
#include <fstream>
#include <sstream>
#include "md_io/VISWriter.h"
#include "md_io/Common.h"
#include "datastructures/ParticleContainer.h"
#include "molecules/Molecule.h"

VISWriter::VISWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
   _filename = filename;
   _writeFrequency = writeFrequency;
   _incremental = incremental;
   _numberOfTimesteps = numberOfTimesteps;

   if (filename == "default")
      _filenameisdate = true;
   else
      _filenameisdate = false;
}

VISWriter::~VISWriter(){}

void VISWriter::initOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain) {
}

void VISWriter::doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
   if(simstep%_writeFrequency == 0) {
      stringstream filenamestream;
      if(_filenameisdate) {
         filenamestream << gettimestring();
      } else {
         filenamestream << _filename;
      }
      if(_incremental) {
         unsigned long temp = simstep/_writeFrequency;
	 filenamestream << "-";
         while(temp < floor(_numberOfTimesteps/_writeFrequency)){
            filenamestream << "0";
            temp = temp*10;
         }
         filenamestream << simstep/_writeFrequency << ".vis";
      } else {
	 filenamestream << ".vis";
      }

      ofstream ostrm(filenamestream.str().c_str());

      // originally VIS files had a fixed width of 8 (and no t), here I use 12 (with 2 for t)
      //ostrm << "t           x           y           z          q0          q1          q2          q3" << endl;
      for(Molecule* pos=particleContainer->begin() ; pos != particleContainer->end(); pos=particleContainer->next())
      {
         //printf("%2d %6.2f %6.2f %6.2f %6.3f %6.3f %6.3f %6.3f",pos->componentid(),pos->r(0),pos->r(1),pos->r(2),(pos->q()).qw(),(pos->q()).qx(),(pos->q()).qy(),(pos->q()).qz());
         //ostrm << setw(2) << pos->componentid()
         //      << setw(12) << pos->r(0) << setw(12) << pos->r(1) << setw(12) << pos->r(2)
         //      << setw(12) << (pos->q()).qw() << setw(12) << (pos->q()).qx() << setw(12) << (pos->q()).qy() << setw(12) << (pos->q()).qz()
         //      << endl;
         ostrm << setiosflags(ios::fixed) //<< setiosflags(ios::floatfield) << << setiosflags(ios::showpoint)
               << setw(1) << pos->componentid()
               << setprecision(1) << setw(6) << pos->r(0) << setw(6) << pos->r(1)
               << setw(6) << pos->r(2) << setprecision(3) << setw(7)
               << (pos->q()).qw() << setw(7) << (pos->q()).qx() << setw(7)
               << (pos->q()).qy() << setw(7) << (pos->q()).qz()
               << setw(6);
         
         // RK
	 //if ((pos->clusterid() != -1) && (MINCLUSTERSIZE < m_clusters.find(pos->clusterid())->second))
         //{
         //  ostrm << setiosflags(ios::fixed) << setw(11) << pos->clusterid()+1 << endl;
         //}
         //else
         //{
           ostrm << setiosflags(ios::fixed) << setw(11) << "0" << endl;
         //}
	 // /RK
      }
      ostrm.close();
   }
}

void VISWriter::finishOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain){
}
