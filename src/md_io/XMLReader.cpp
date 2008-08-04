#include "md_io/XMLReader.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "XMLReader_main.h"

#include <string>
#include <fstream>
#include <iostream>

md_io::XMLReader::XMLReader() {
}

md_io::XMLReader::~XMLReader(){}

void md_io::XMLReader::setPhaseSpaceFile(string filename) {
  _phaseSpaceFileName = filename;
}

void md_io::XMLReader::setPhaseSpaceHeaderFile(string filename) {
  _phaseSpaceHeaderFileName = filename;
}

void md_io::XMLReader::readPhaseSpaceHeader(Domain* domain) {
  md_io::XMLReader_main * _xmlreader = new md_io::XMLReader_main();
  TiXmlDocument XMLdoc = _xmlreader->XMLReader_main_get_doc_chara(_phaseSpaceHeaderFileName.c_str());
  TiXmlDocument *XMLdoc_p = &XMLdoc;

  vector<Component>& dcomponents = domain->getComponents();

  domain->setinpversion(0);
  domain->setinpversion(_xmlreader->Eval_ul(XMLdoc_p,
     (string)"/mardyncfg/header/version/text()"));
  if(domain->getinpversion() < 20070801 && domain->getlocalRank() == 0)
  {
    cerr << "Input: old version: " << domain->getinpversion() << endl;
    exit(1);
  }

  domain->setCurrentTime(_xmlreader->Eval_d(XMLdoc_p,
     (string)"/mardyncfg/experiment/current-time/text()"));

  domain->setGlobalTemperature(_xmlreader->Eval_d(XMLdoc_p,
     (string)"/mardyncfg/experiment/temperature/text()"));

  domain->setGlobalLength(0, _xmlreader->Eval_d(XMLdoc_p,
     (string)"/mardyncfg/experiment/length/x/text()"));
  domain->setGlobalLength(1, _xmlreader->Eval_d(XMLdoc_p,
     (string)"/mardyncfg/experiment/length/y/text()"));
  domain->setGlobalLength(2, _xmlreader->Eval_d(XMLdoc_p,
     (string)"/mardyncfg/experiment/length/z/text()"));

  int numcomponents = _xmlreader->Eval_i(XMLdoc_p,
     (string)"/mardyncfg/experiment/components/data/components/@amount");
  dcomponents.resize(numcomponents);
  
  
  // TODO: parse the rest
  double x,y,z;
  unsigned int i,j;
  double xi,eta;
  std::ostringstream itoastr1,itoastr2;

  for(i=0;i<numcomponents;++i)
  {
    dcomponents[i].setID(i);
    unsigned int numljcenters=0;
    unsigned int numdipoles=0;
    unsigned int numquadrupoles=0;
    itoastr1 << (i+1);
    numljcenters = _xmlreader->Eval_i(XMLdoc_p,
       (string)"count(/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/ljcenter)");
    numdipoles = _xmlreader->Eval_i(XMLdoc_p,
       (string)"count(/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole)");
    numquadrupoles = _xmlreader->Eval_i(XMLdoc_p,
       (string)"count(/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole)");
    itoastr1.str("");

    for(j=0;j<numljcenters;++j)
    {
      double m,sigma,eps;
      itoastr1 << (i+1);
      itoastr2 << (j+1);
      x = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/ljcenter["+itoastr2.str()+"]/x/text()");
      y = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/ljcenter["+itoastr2.str()+"]/y/text()");
      z = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/ljcenter["+itoastr2.str()+"]/z/text()");
      m = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/ljcenter["+itoastr2.str()+"]/m/text()");
      eps = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/ljcenter["+itoastr2.str()+"]/eps/text()");
      sigma = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/ljcenter["+itoastr2.str()+"]/sigma/text()");
      dcomponents[i].addLJcenter(x,y,z,m,eps,sigma);

      itoastr1.str("");
      itoastr2.str("");
    }
    for(j=0;j<numdipoles;++j)
    {
      double eMyx,eMyy,eMyz,absMy;
      itoastr1 << (i+1);
      itoastr2 << (j+1);
      x = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole["+itoastr2.str()+"]/x/text()");
      y = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole["+itoastr2.str()+"]/y/text()");
      z = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole["+itoastr2.str()+"]/z/text()");
      eMyx = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole["+itoastr2.str()+"]/eMyx/text()");
      eMyy = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole["+itoastr2.str()+"]/eMyy/text()");
      eMyz = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole["+itoastr2.str()+"]/eMyz/text()");
      absMy = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dipole["+itoastr2.str()+"]/absMy/text()");
      dcomponents[i].addDipole(x,y,z,eMyx,eMyy,eMyz,absMy);
      itoastr1.str("");
      itoastr2.str("");
    }
    for(j=0;j<numquadrupoles;++j)
    {
      double eQx,eQy,eQz,absQ;
      itoastr1 << (i+1);
      itoastr2 << (j+1);
      x = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole["+itoastr2.str()+"]/x/text()");
      y = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole["+itoastr2.str()+"]/y/text()");
      z = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole["+itoastr2.str()+"]/z/text()");
      eQx = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole["+itoastr2.str()+"]/eQx/text()");
      eQy = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole["+itoastr2.str()+"]/eQy/text()");
      eQz = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole["+itoastr2.str()+"]/eQz/text()");
      absQ = _xmlreader->Eval_d(XMLdoc_p,
         (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/quadrupole["+itoastr2.str()+"]/absQ/text()");
      dcomponents[i].addQuadrupole(x,y,z,eQx,eQy,eQz,absQ);
      itoastr1.str("");
      itoastr2.str("");
    }
    double IDummy1,IDummy2,IDummy3;
    itoastr1 << (i+1);
    IDummy1 = _xmlreader->Eval_d(XMLdoc_p,
      (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dummy1/text()");
    IDummy2 = _xmlreader->Eval_d(XMLdoc_p,
      (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dummy2/text()");
    IDummy3 = _xmlreader->Eval_d(XMLdoc_p,
      (string)"/mardyncfg/experiment/components/data/components/comp["+itoastr1.str()+"]/dummy3/text()");
    if(IDummy1>0.) dcomponents[i].setI11(IDummy1);
    if(IDummy2>0.) dcomponents[i].setI22(IDummy2);
    if(IDummy3>0.) dcomponents[i].setI33(IDummy3);
    itoastr1.str("");
  }
  vector<double>& dmixcoeff = domain->getmixcoeff();
  dmixcoeff.clear();
  unsigned int imax = ((numcomponents*numcomponents)/2) - (numcomponents/2) ;
  for(i=1;i<=imax;++i)
  {
     itoastr1 << i;
     xi = _xmlreader->Eval_d(XMLdoc_p,
        (string)"/mardyncfg/experiment/components/data/components/mixcoeff/xi["+itoastr1.str()+"]/text()");
     eta = _xmlreader->Eval_d(XMLdoc_p,
        (string)"/mardyncfg/experiment/components/data/components/mixcoeff/eta["+itoastr1.str()+"]/text()");
     dmixcoeff.push_back(xi);
     dmixcoeff.push_back(eta);
     itoastr1.str("");
  }
  domain->setepsilonRF(_xmlreader->Eval_d(XMLdoc_p,
     (string)"/mardyncfg/experiment/components/data/components/epsilon-rf/text()"));
  
}

void md_io::XMLReader::readPhaseSpace(datastructures::ParticleContainer<Molecule>* particleContainer, Domain* domain) {
  // won't be implemented (doesn't make sense)
}
