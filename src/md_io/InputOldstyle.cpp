// InputOldstyle.cpp

#include "md_io/InputOldstyle.h"
#include "Domain.h"
#include "datastructures/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"
#include "Logger.h"
#include <climits>

using Log::global_log;


InputOldstyle::InputOldstyle() {
}

InputOldstyle::~InputOldstyle(){}

void InputOldstyle::setPhaseSpaceFile(string filename) {
 _phaseSpaceFileStream.open(filename.c_str());
}

void InputOldstyle::setPhaseSpaceHeaderFile(string filename) {

}

void InputOldstyle::readPhaseSpaceHeader(Domain* domain, double timestep, double cutoff)
{
  vector<Component>& dcomponents = domain->getComponents();
  string token;
  _phaseSpaceFileStream >> token;
  domain->setinpversion(0);

  if((token != "mardyn") && (token != "MOLDY") && (token != "ls1r1") && (token != "mrdyn") && (token != "MDProject") && (token != "Mardyn") && (token != "MARDYN"))
  {
    if(domain->ownrank() == 0)
       global_log->error() << "Input: NOT AN OldStyle LS1 MARDYN INPUT! (starts with "
            << token << ")" << endl;
    exit(1);
  }
  _phaseSpaceFileStream >> token;
  if(token != "trunk")
  {
     global_log->error() << "Wrong input file version (\'"
	               << token << "\' instead of \'trunk\'). Aborting.\n";
     exit(787);
  }
  _phaseSpaceFileStream >> token;
  domain->setinpversion(strtoul(token.c_str(),NULL,0));
  if(domain->getinpversion() < 20080701 && domain->getlocalRank() == 0)
  {
	  global_log->error() << "Input: OLD VERSION (" << domain->getinpversion() << ")" << endl;
	  exit(1);
  }

  char c;
  double x,y,z;
  unsigned int numcomponents=0;
  unsigned int j;
  double m,sigma,eps;
  double xi,eta;
  unsigned long i;
  double do_shift;

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

    if((token == "currentTime") || (token == "t"))
    {
      _phaseSpaceFileStream >> token;
      domain->setCurrentTime(strtod(token.c_str(), NULL));
    }
    else if((token == "Temperature") || (token == "T"))
    {
       domain->disableCT();
       double targetT;
       _phaseSpaceFileStream >> targetT;
       domain->setGlobalTemperature(targetT);
    }
    else if((token == "ThermostatTemperature") || (token == "ThT") || (token == "h"))
    {
       int i;
       double targetT;
       _phaseSpaceFileStream >> i;
       _phaseSpaceFileStream >> targetT;
       domain->setTargetT(i, targetT);
    }
    else if((token == "ComponentThermostat") || (token == "CT") || (token == "o"))
    {
       if(!domain->severalThermostats())
	  domain->enableCT();
       int cid;
       _phaseSpaceFileStream >> cid;
       cid--;
       int th;
       _phaseSpaceFileStream >> th;
       if(0 >= th) continue;
       domain->setComponentThermostat(cid, th);
    }
    else if((token == "Undirected") || (token == "U"))
    {
       int tst;
       _phaseSpaceFileStream >> tst;
       domain->enableUndirectedThermostat(tst);
    }
    else if((token == "Length") || (token == "L"))
    {
       double globalLength[3];
       _phaseSpaceFileStream >> globalLength[0] >> globalLength[1] >> globalLength[2];
       for(int i=0; i < 3; i++)
          domain->setGlobalLength(i, 1.0000000333 * globalLength[i]);
    }
    else if((token == "NumberOfComponents") || (token == "C"))
    {
      _phaseSpaceFileStream >> numcomponents;
      dcomponents.resize(numcomponents);
      for(i=0;i<numcomponents;++i)
      {
        dcomponents[i].setID(i);
        unsigned int numljcenters = 0;
	unsigned int numcharges = 0;
        unsigned int numdipoles = 0;
        unsigned int numquadrupoles = 0;
	unsigned int numtersoff = 0;
        _phaseSpaceFileStream >> numljcenters >> numcharges
	                      >> numdipoles >> numquadrupoles
			      >> numtersoff;
        for(j=0;j<numljcenters;++j)
        {
          _phaseSpaceFileStream >> x >> y >> z >> m >> eps >> sigma >> cutoff >> do_shift;
          dcomponents[i].addLJcenter(x, y, z, m, eps, sigma, cutoff, (do_shift != 0));
        }
        for(j = 0; j < numcharges; j++)
        {
          double q;
          _phaseSpaceFileStream >> x >> y >> z >> m >> q;
          dcomponents[i].addCharge(x, y, z, m, q);
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
        for(j = 0; j < numtersoff; j++)
        {
          double x, y, z, m, A, B, lambda, mu, R, S, c, d, h, n, beta;
          _phaseSpaceFileStream >> x >> y >> z;
          _phaseSpaceFileStream >> m >> A >> B;
          _phaseSpaceFileStream >> lambda >> mu >> R >> S;
          _phaseSpaceFileStream >> c >> d >> h >> n >> beta;
          dcomponents[i].addTersoff(
            x, y, z,
            m, A, B,
            lambda, mu, R, S,
            c, d, h, n, beta
          );
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
      header = false;
    }
    else if((token == "NumberOfMolecules") || (token == "N"))
    {
      _phaseSpaceFileStream >> token;
      domain->setglobalNumMolecules(strtoul(token.c_str(),NULL,0));
    }
    else if((token == "AssignCoset") || (token == "S"))
    {
      unsigned cid, cosetid;
      _phaseSpaceFileStream >> cid >> cosetid;
      cid--;
      domain->assignCoset(cid, cosetid);
    }
    else if((token == "Accelerate") || (token == "A"))
    {
       unsigned cosetid;
       _phaseSpaceFileStream >> cosetid;
       double v[3];
       for(unsigned d = 0; d < 3; d++) _phaseSpaceFileStream >> v[d];
       double tau;
       _phaseSpaceFileStream >> tau;
       double ainit[3];
       for(unsigned d = 0; d < 3; d++) _phaseSpaceFileStream >> ainit[d];
       domain->specifyComponentSet(cosetid, v, tau, ainit, timestep);
    }
    else
    {
       global_log->error() << "Invalid token \'" << token << "\' found. Skipping rest of the header." << endl;
       header = false;
    }
  }
}

unsigned long InputOldstyle::readPhaseSpace(ParticleContainer* particleContainer, double cutoffRadius, list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {
  global_log->set_mpi_output_root(0);

  string token;

  vector<Component>& dcomponents = domain->getComponents();

  double x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz;
  double Fx,Fy,Fz,Mx,My,Mz;
  unsigned int numcomponents=dcomponents.size();
  unsigned long i,id;
  int componentid;
  
  unsigned long maxid = 0;

  x=y=z=vx=vy=vz=q1=q2=q3=Dx=Dy=Dz=0.;
  q0=1.;
  Fx=Fy=Fz=Mx=My=Mz=0.;
  _phaseSpaceFileStream >> token;
  if((token == "NumberOfMolecules") || (token == "N"))
  {
    string nummolecules;
    _phaseSpaceFileStream >> nummolecules;
    domain->setglobalNumMolecules(strtoul(nummolecules.c_str(),NULL,0));
    _phaseSpaceFileStream >> token;
  }

  if((token=="MoleculeFormat") || (token == "M")) {
    string ntypestring("ICRVQD");
    enum Ndatatype { ICRVQD, IRV, ICRV, ICRVFQDM } ntype=ICRVQD;

    global_log->info() << "Reading " << domain->getglobalNumMolecules() << " molecules" << endl;
    if(domain->getinpversion() >= 51129) _phaseSpaceFileStream >> ntypestring;
    ntypestring.erase(ntypestring.find_last_not_of(" \t\n")+1);
    ntypestring.erase(0,ntypestring.find_first_not_of(" \t\n"));
    if (ntypestring=="ICRVFQDM")
      ntype=ICRVFQDM;
    else if (ntypestring=="ICRV")
      ntype=ICRV;
    else if (ntypestring=="IRV")
      ntype=IRV;
    global_log->info() << "Molecul format: " << ntypestring << endl;
    if(!numcomponents)
    {
      global_log->warning() << "No components defined! Setting up single one-centered LJ" << endl;
      numcomponents=1;
      dcomponents.resize(numcomponents);
      dcomponents[0].setID(0);
      dcomponents[0].addLJcenter(0.,0.,0.,1.,1.,1.,cutoffRadius,false);
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

      if(componentid > numcomponents) {
        global_log->error() << "Molecule id " << id << " has wrong componentid: " << componentid << ">" << numcomponents << endl;
	exit(1);
      }
      componentid --;

      //  store only those molecules within the domain of this process

      // @todo Pointer!!! new!!!
      Molecule m1 = Molecule(id,componentid,x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz,&dcomponents);
      particleContainer->addParticle(m1);
      //(_molecules.back()).setFM(Fx,Fy,Fz,Mx,My,Mz);
      dcomponents[componentid].incrnumMolecules();
      domain->setglobalRotDOF(dcomponents[componentid].rot_dof()+domain->getglobalRotDOF());
      if(id > maxid) maxid = id;
      std::list<ChemicalPotential>::iterator cpit;
      for(cpit = lmu->begin(); cpit != lmu->end(); cpit++)
      {
         if( !cpit->hasSample() &&
	     (componentid == cpit->getComponentID()) )
	 {
	    cpit->storeMolecule(m1);
	 }
      }
      unsigned long iph = domain->getglobalNumMolecules() / 100;
      if(i%iph == 0)
        global_log->info() << "Finished reading molecules: " << i/iph << "%" << endl;
    }
    global_log->info() << "Finished reading molecules: 100%" << endl;
    global_log->info() << "Reading Molecules done" << endl;

    if(!domain->getglobalRho()){
      domain->setglobalRho(domain->getglobalNumMolecules()/(domain->getGlobalLength(0)*domain->getGlobalLength(1)*domain->getGlobalLength(2)));
      global_log->info() << "Calculated Rho_global = " << domain->getglobalRho() << endl;
    }

    _phaseSpaceFileStream.close();

  }
  else {
    global_log->error() << "Error in the PhaseSpace File. expected 'MoleculeFormat' or 'M' but found token " << token << endl;
    exit(1);
  }

  return maxid;
}
