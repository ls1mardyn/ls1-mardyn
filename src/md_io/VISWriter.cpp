#include <iomanip>
#include <fstream>
#include <sstream>
#include "md_io/VISWriter.h"
#include "md_io/Common.h"
#include "datastructures/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"

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
                           DomainDecompBase* domainDecomp, Domain* domain)
{
   this->_wroteVIS = false;
}

void VISWriter::doOutput( ParticleContainer* particleContainer,
                          DomainDecompBase* domainDecomp, Domain* domain,
			  unsigned long simstep, list<ChemicalPotential>* lmu) 
{
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
         while(temp < floor((double) _numberOfTimesteps / (double) _writeFrequency)){
            filenamestream << "0";
            temp = temp*10;
         }
         filenamestream << simstep/_writeFrequency << ".vis_";
      } else {
	 filenamestream << ".vis_";
      }

      ofstream visittfstrm(filenamestream.str().c_str());

      if(!domain->getlocalRank() && !this->_wroteVIS)
      {
         visittfstrm << "      id t          x          y          z     q0     q1     q2     q3        c\n";
         this->_wroteVIS = true;
      }
      else if(!domain->getlocalRank()) visittfstrm << "#\n";

      // originally VIS files had a fixed width of 8 (and no t), here I use 12 (with 2 for t)
      //ostrm << "t           x           y           z          q0          q1          q2          q3" << endl;
      for(Molecule* pos=particleContainer->begin() ; pos != particleContainer->end(); pos=particleContainer->next())
      {
         bool halo = false;
         for(unsigned short d = 0; d < 3; d++)
         {
            if((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d)))
            {
               halo = true;
               break;
            }
         }
         if(!halo)
	 {
            visittfstrm << setiosflags(ios::fixed) << setw(8) << pos->id() << setw(2)
                       << pos->componentid() << setprecision(3);
            for(unsigned short d = 0; d < 3; d++) visittfstrm << setw(11) << pos->r(d);
            visittfstrm << setprecision(3) << setw(7) << pos->q().qw() << setw(7) << pos->q().qx()
                        << setw(7) << pos->q().qy()<< setw(7) << pos->q().qz()
                        << setw(9) << right << 0 << "\n";
	 }
      }
      visittfstrm.close();
   }
}

void VISWriter::finishOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain){
}
