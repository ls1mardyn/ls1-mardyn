#include "md_io/PovWriter.h"
#include "md_io/Common.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "molecules/Molecule.h"

md_io::PovWriter::PovWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
   _filename = filename;
   _writeFrequency = writeFrequency;
   _incremental = incremental;
   _numberOfTimesteps = numberOfTimesteps;

   if (filename == "default")
      _filenameisdate = true;
   else
      _filenameisdate = false;
}

md_io::PovWriter::~PovWriter(){}

void md_io::PovWriter::initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain){
}

void md_io::PovWriter::doOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep){
   if(simstep%_writeFrequency == 0) {
      
      stringstream filenamestream;
      filenamestream << _filename;
      if(_filenameisdate) {
         filenamestream << gettimestring();
      }
      if(_incremental) {
         unsigned long temp = simstep/_writeFrequency;
         filenamestream << "-";
         while(temp < floor(_numberOfTimesteps/_writeFrequency)){
            filenamestream << "0";
            temp = temp*10;
         }
         filenamestream << simstep/_writeFrequency << ".pov";
      } else {
         filenamestream << ".pov";
      }

      ofstream ostrm(filenamestream.str().c_str());
   
      ostrm << "// " << filenamestream.str() << endl;
      ostrm << "// moldy" << endl;
      time_t now;
      now=time(NULL);
      ostrm << "// " << ctime(&now) << endl;
   
      ostrm << "// bb: [0," << domain->getGlobalLength(0) << "]^3" << endl;
      ostrm << "//*PMRawBegin" << endl;
      ostrm << "background {rgb <1,1,1>}" << endl;
      ostrm << "//*PMRawEnd" << endl;
      vector<Component> dcomponents = domain->getComponents();
      for(unsigned int i=0;i<dcomponents.size();++i)
      {
              ostringstream osstrm;
              osstrm.clear(); osstrm.str("");
              osstrm << " pigment {color rgb <" << (i+1)%2 << "," << (i+1)/2%2 << "," << (i+1)/4%2 << ">}";
              osstrm << " finish{ambient 0.5 diffuse 0.4 phong 0.3 phong_size 3}";
              ostrm << "#declare T" << i << " = ";
              //ostrm << "sphere {<0,0,0>,0.5 pigment {color rgb<1,0,0>} finish{ambient 0.5 diffuse 0.4 phong 0.3 phong_size 3} scale 1.}";
   	   dcomponents.at(i).writePOVobjs(ostrm,osstrm.str());
              ostrm << endl;
      }
      ostrm << endl;
      ostrm << "camera { perspective" << endl;
      float xloc=-.1*domain->getGlobalLength(0);
      float yloc=1.1*domain->getGlobalLength(1);
      float zloc=-1.5*domain->getGlobalLength(2);
      ostrm << " location <" << xloc << ", " << yloc << ", " << zloc << ">" << endl;
      //ostrm << " direction <0, 0, 1>" << endl;
      //ostrm << " right <1.33333, 0, 0>" << endl;
      //ostrm << " up <0, 1, 0>" << endl;
      //ostrm << " sky <0, 1, 0>" << endl;
      ostrm << " look_at <" << .5*domain->getGlobalLength(0) << ", " << .5*domain->getGlobalLength(1) << ", " << .5*domain->getGlobalLength(2) << ">" << endl;
      ostrm << "}" << endl;
      ostrm << endl;
      ostrm << "light_source { <" << xloc << ", " << yloc << ", " << zloc << ">, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <0,0,0>, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <0,0," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <0," << domain->getGlobalLength(1) << ",0>, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <0," << domain->getGlobalLength(1) << "," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <" << domain->getGlobalLength(0) << ",0,0>, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <" << domain->getGlobalLength(0) << ",0," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <" << domain->getGlobalLength(0) << "," << domain->getGlobalLength(1) << ",0>, color rgb <1,1,1> }" << endl;
      ostrm << "light_source { <" << domain->getGlobalLength(0) << "," << domain->getGlobalLength(1) << "," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
      ostrm << endl;
      ostrm << "// " << dcomponents.size() << " objects for the atoms following..." << endl;
      double mrot[3][3];
      for(Molecule* pos=particleContainer->begin() ; pos != particleContainer->end(); pos=particleContainer->next())
      {
              (pos->q()).getRotinvMatrix(mrot);
              //cout << "object { T0 rotate <0,0,0> translate <0,0,0>}" << endl;
              ostrm << "object { T" << pos->componentid();
              ostrm << " matrix <"
                    << mrot[0][0] << "," << mrot[0][1] << "," << mrot[0][2] << ","
                    << mrot[1][0] << "," << mrot[1][1] << "," << mrot[1][2] << ","
                    << mrot[2][0] << "," << mrot[2][1] << "," << mrot[2][2] << ","
                    << pos->r(0) << "," << pos->r(1) << "," << pos->r(2)
                    << ">";
              ostrm << "}" << endl;
      }

      // RK
      /* map cluster ID to color */
      //if ((pos->clusterid() != -1) && (m_clusters.find(pos->clusterid())->second > MINCLUSTERSIZE))
      //{
      //  ostrm << "  pigment {color rgb<" <<
      //      (0.25*(pos->clusterid()%5)) << "," <<
      //      (0.25*((pos->clusterid()/5)%5)) << "," <<
      //      (0.25*((pos->clusterid()/125)%5)) << "," <<
      //      ">}";
      //} else
      //{
      //  ostrm << "  pigment {color rgb<0.9,0.9,0.9>}";
      //}
      // /RK
   
      ostrm.close();
   }
}

void md_io::PovWriter::finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain){
}
