#include "md_io/InputOldstyle.h"
#include "datastructures/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include <climits>

InputOldstyle::InputOldstyle() {
}

InputOldstyle::~InputOldstyle(){}

void InputOldstyle::setPhaseSpaceFile(string filename) {
 _phaseSpaceFileStream.open(filename.c_str());
}

void InputOldstyle::setPhaseSpaceHeaderFile(string filename) {

}

void InputOldstyle::readPhaseSpaceHeader(Domain* domain) {
  vector<Component>& dcomponents = domain->getComponents();
  string token;
  _phaseSpaceFileStream >> token;
  domain->setinpversion(0);
  if(token != "MDProject" && domain->getlocalRank() == 0)
  {
    cerr << "Input: NOT A MOLDY INPUT! (starts with " << token << ")" << endl;
  }
  else
  {
    _phaseSpaceFileStream >> token;
    domain->setinpversion(strtoul(token.c_str(),NULL,0));
    if(domain->getinpversion() < 20070111 && domain->getlocalRank() == 0)
    {
      cerr << "Input: OLD VERSION (" << domain->getinpversion() << ")" << endl;
    }
  }

  char c;
  double x,y,z;
  unsigned int numcomponents=0;
  unsigned int j;
  double m,sigma,eps;
  double xi,eta;
  unsigned long i;

  // When the last header element is reached, "header" is set to false
  bool header = true;

  while(header) {
    _phaseSpaceFileStream >> c;
    if(c=='#') {
      _phaseSpaceFileStream.ignore(INT_MAX,'\n');
      continue;
    }
    else {
      _phaseSpaceFileStream.putback(c);
    }
    token.clear();
    _phaseSpaceFileStream >> token;

    if(token=="currentTime"){
      _phaseSpaceFileStream >> token;
      domain->setCurrentTime(strtod(token.c_str(), NULL));
    }

    if(token=="Temperature"){
      _phaseSpaceFileStream >> token;
      domain->setGlobalTemperature(strtod(token.c_str(), NULL));
    }
    else if(token=="Length"){
      _phaseSpaceFileStream >> token;
      domain->setGlobalLength(0, strtod(token.c_str(),NULL));
      _phaseSpaceFileStream >> token;
      domain->setGlobalLength(1, strtod(token.c_str(),NULL));
      _phaseSpaceFileStream >> token;
      domain->setGlobalLength(2, strtod(token.c_str(),NULL));
    }
    if(token=="NumberOfComponents") {

      _phaseSpaceFileStream >> numcomponents;
      dcomponents.resize(numcomponents);
      for(i=0;i<numcomponents;++i)
      {
        dcomponents[i].setID(i);
        unsigned int numljcenters=0;
        unsigned int numdipoles=0;
        unsigned int numquadrupoles=0;
        _phaseSpaceFileStream >> numljcenters >> numdipoles >> numquadrupoles;
        for(j=0;j<numljcenters;++j)
        {
          _phaseSpaceFileStream >> x >> y >> z >> m >> eps >> sigma;
          dcomponents[i].addLJcenter(x,y,z,m,eps,sigma,2.5*sigma,true);
        }

        for(j=0;j<numdipoles;++j)
        {
          double eMyx,eMyy,eMyz,absMy;
          _phaseSpaceFileStream >> x >> y >> z >> eMyx >> eMyy >> eMyz >> absMy;
          dcomponents[i].addDipole(x,y,z,eMyx,eMyy,eMyz,absMy);
        }
        for(j=0;j<numquadrupoles;++j)
        {
          double eQx,eQy,eQz,absQ;
          _phaseSpaceFileStream >> x >> y >> z >> eQx >> eQy >> eQz >> absQ;
          dcomponents[i].addQuadrupole(x,y,z,eQx,eQy,eQz,absQ);
        }
        double IDummy1,IDummy2,IDummy3;
        _phaseSpaceFileStream >> IDummy1 >> IDummy2 >> IDummy3;
        if(IDummy1>0.) dcomponents[i].setI11(IDummy1);
        if(IDummy2>0.) dcomponents[i].setI22(IDummy2);
        if(IDummy3>0.) dcomponents[i].setI33(IDummy3);
      }
      vector<double>& dmixcoeff = domain->getmixcoeff();
      dmixcoeff.clear();
      for(i=0;i<numcomponents-1;++i)
      {
        for(j=i+1;j<numcomponents;++j)
        {
          _phaseSpaceFileStream >> xi >> eta;
          dmixcoeff.push_back(xi);
          dmixcoeff.push_back(eta);
        }
      }
      _phaseSpaceFileStream >> token;
      domain->setepsilonRF(strtod(token.c_str(),NULL));
    }
    else if(token=="NumberOfMolecules"){
      _phaseSpaceFileStream >> token;
      domain->setglobalNumMolecules(strtoul(token.c_str(),NULL,0));
      header = false;
    }
  }
}

void InputOldstyle::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {

  string token;

  vector<Component>& dcomponents = domain->getComponents();

  double x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz;
  double Fx,Fy,Fz,Mx,My,Mz;
  unsigned int numcomponents=dcomponents.size();
  unsigned long i,id;
  int componentid;

  x=y=z=vx=vy=vz=q1=q2=q3=Dx=Dy=Dz=0.;
  q0=1.;
  Fx=Fy=Fz=Mx=My=Mz=0.;
  _phaseSpaceFileStream >> token;
  if(token=="NumberOfMolecules") {
    string nummolecules;
    _phaseSpaceFileStream >> nummolecules;
    domain->setglobalNumMolecules(strtoul(nummolecules.c_str(),NULL,0));
    _phaseSpaceFileStream >> token;
  }

  if(token=="MoleculeFormat") {
    string ntypestring("ICRVQD");
    enum Ndatatype { ICRVQD, IRV, ICRV, ICRVFQDM } ntype=ICRVQD;

    if(domain->getlocalRank()==0) cout << "reading " << domain->getglobalNumMolecules() << " molecules" << flush;
    if(domain->getinpversion() >= 51129) _phaseSpaceFileStream >> ntypestring;
    ntypestring.erase(ntypestring.find_last_not_of(" \t\n")+1);
    ntypestring.erase(0,ntypestring.find_first_not_of(" \t\n"));
    if (ntypestring=="ICRVFQDM")
      ntype=ICRVFQDM;
    else if (ntypestring=="ICRV")
      ntype=ICRV;
    else if (ntypestring=="IRV")
      ntype=IRV;
    if(domain->getlocalRank()==0) cout << " (" << ntypestring << ")" << flush;
    if(!numcomponents)
    {
      if(domain->getlocalRank()==0) cout << endl << "No components defined! Setting up single one-centered LJ" << endl;
      numcomponents=1;
      dcomponents.resize(numcomponents);
      dcomponents[0].setID(0);
      dcomponents[0].addLJcenter(0.,0.,0.,1.,1.,1.,2.5,true);
    }
    //m_molecules.clear();
    for(i=0;i<domain->getglobalNumMolecules();++i)
    {
      if(ntype==ICRVFQDM)
        _phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz >> Fx >> Fy >> Fz
              >> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz >> Mx >> My >> Mz;
      else if(ntype==ICRV)
        _phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz;
      else if(ntype==IRV)
        _phaseSpaceFileStream >> id >> x >> y >> z >> vx >> vy >> vz;
      else
        _phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
              >> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz;
      if((x<0.0 || x>=domain->getGlobalLength(0) || y<0.0 || y>=domain->getGlobalLength(1) || z<0.0 || z>=domain->getGlobalLength(2)) && domain->getlocalRank() == 0) {
        cerr << endl << id << ": Molecule " << x << ";" << y << ";" << z << " out of box! " << flush;
      }

      if(((int)componentid > (int)numcomponents) && domain->getlocalRank() == 0) {
        cerr << "Molecule id " << id << " has wrong componentid: " << componentid << ">" << numcomponents << endl;
      }
      // @todo why do componetids start with 0?
      --componentid;

      //  store only those molecules within the domain of this process

      // @todo Pointer!!! new!!!
      Molecule m1 = Molecule(id,componentid,x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz,&dcomponents);
      particleContainer->addParticle(m1);
      //(_molecules.back()).setFM(Fx,Fy,Fz,Mx,My,Mz);
      dcomponents[componentid].incrnumMolecules();
      domain->setglobalRotDOF(dcomponents[componentid].rot_dof()+domain->getglobalRotDOF());

      if(!(i%1000)) if(domain->getlocalRank()==0) cout << '.' << flush;
    }
    if(domain->getlocalRank()==0) cout << " done" << endl;

    if(!domain->getglobalRho()){
      domain->setglobalRho(domain->getglobalNumMolecules()/(domain->getGlobalLength(0)*domain->getGlobalLength(1)*domain->getGlobalLength(2)));
      if(domain->getlocalRank()==0) cout << "calculated global Rho:\t" << domain->getglobalRho() << endl;
    }

    _phaseSpaceFileStream.close();

  }
  else {
    cerr << "ERROR in AsciiReader::readPhaseSpace: Error in the PhaseSpace File" << endl;
  }

}
