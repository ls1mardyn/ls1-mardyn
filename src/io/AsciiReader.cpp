#include "io/AsciiReader.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"
#include "ensemble/PressureGradient.h"
#include "Domain.h"

#include <climits>

using namespace std;

AsciiReader::AsciiReader() {
}

AsciiReader::~AsciiReader(){}

void AsciiReader::setPhaseSpaceFile(string filename) {
#ifdef ENABLE_MPI
 _phaseSpaceFileStream.str(filename);
#else
 _phaseSpaceFileName = filename;
 _phaseSpaceFileStream.open(filename.c_str());
#endif
}

void AsciiReader::setPhaseSpaceHeaderFile(string filename) {
 // if the components are internal there's no need to open another ifstream
 if (_phaseSpaceFileName != filename)
		_phaseSpaceHeaderFileStream.open(filename.c_str());
 _phaseSpaceHeaderFileName = filename;
}

void AsciiReader::readPhaseSpaceHeader(Domain* domain, double timestep)
{
	vector<Component>& dcomponents = domain->getComponents();
	string token;
	_phaseSpaceFileStream >> token;
	domain->setinpversion(0);

	if((token != "mardyn") && (token != "MOLDY") && (token != "ls1r1") && (token != "mrdyn") && (token != "MDProject") && (token != "Mardyn") && (token != "MARDYN"))
	{
		if(domain->ownrank() == 0) cerr << "Input: NOT AN Ascii LS1 MARDYN INPUT! (starts with " << token << ")" << endl;
		exit(1);
	}
	_phaseSpaceFileStream >> token;
	if(token != "trunk")
	{
		 if(domain->ownrank() == 0) cout << "Wrong input file version (\'"
	               << token << "\' instead of \'trunk\'). Aborting.\n";
		 exit(787);
	}
	_phaseSpaceFileStream >> token;
	domain->setinpversion(strtoul(token.c_str(),NULL,0));
	if(domain->getinpversion() < 20080701 && domain->ownrank() == 0)
	{
	    cerr << "Input: OLD VERSION (" << domain->getinpversion() << ")" << endl;
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
			 _phaseSpaceFileStream >> token;
			 domain->disableComponentwiseThermostat();
			 domain->setGlobalTemperature(strtod(token.c_str(), NULL));
		}
		else if((token == "ThermostatTemperature") || (token == "ThT") || (token == "h"))
		{
			 int i;
			 double targetT;
			 _phaseSpaceFileStream >> i;
			 _phaseSpaceFileStream >> targetT;
			 domain->setTargetTemperature(i, targetT);
		}
		else if((token == "ComponentThermostat") || (token == "CT") || (token == "o"))
		{
			 if(!domain->severalThermostats())
	  domain->enableComponentwiseThermostat();
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
			    domain->setGlobalLength(i, globalLength[i]);
		}
		else if((token == "HeatCapacity") || (token == "cv") || (token == "I"))
		{
			 unsigned N;
			 double U, UU;
			 _phaseSpaceFileStream >> N >> U >> UU;
			 domain->init_cv(N, U, UU);
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
				double tcutoff = 0.0;
				_phaseSpaceFileStream >> numljcenters >> numcharges
	                      >> numdipoles >> numquadrupoles
			      >> numtersoff;
				for(j=0;j<numljcenters;++j)
				{
					_phaseSpaceFileStream >> x >> y >> z >> m >> eps >> sigma >> tcutoff >> do_shift;
					dcomponents[i].addLJcenter(x, y, z, m, eps, sigma, tcutoff, (do_shift != 0));
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
		}
		else if((token == "NumberOfMolecules") || (token == "N"))
		{
			_phaseSpaceFileStream >> token;
			domain->setglobalNumMolecules(strtoul(token.c_str(),NULL,0));
			header = false;
		}
		else if((token == "AssignCoset") || (token == "S"))
		{
			unsigned cid, cosetid;
			_phaseSpaceFileStream >> cid >> cosetid;
			cid--;
			domain->getPG()->assignCoset(cid, cosetid);
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
			 domain->getPG()->specifyComponentSet(cosetid, v, tau, ainit, timestep);
		}
	}
}

unsigned long AsciiReader::readPhaseSpace(ParticleContainer* particleContainer, list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {
	
	string token;

	vector<Component>& dcomponents = domain->getComponents();
	
	double x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz;
	double Fx,Fy,Fz,Mx,My,Mz;
	unsigned int numcomponents=dcomponents.size();
	unsigned long i,id;
	unsigned int componentid;
	
	unsigned long maxid = 0;
	
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
	
	if((token=="MoleculeFormat") || (token == "M")) {
		string ntypestring("ICRVQD");
		enum Ndatatype { ICRVQD, IRV, ICRV, ICRVFQDM } ntype=ICRVQD;

		if(domain->ownrank() == 0) cout << "reading " << domain->getglobalNumMolecules() << " molecules" << flush;
		if(domain->getinpversion() >= 51129) _phaseSpaceFileStream >> ntypestring;
		ntypestring.erase(ntypestring.find_last_not_of(" \t\n")+1);
		ntypestring.erase(0,ntypestring.find_first_not_of(" \t\n"));
		if (ntypestring=="ICRVFQDM")
			ntype=ICRVFQDM;
		else if (ntypestring=="ICRV")
			ntype=ICRV;
		else if (ntypestring=="IRV")
			ntype=IRV;
		if(domain->ownrank() == 0) cout << " (" << ntypestring << ")" << flush;
		if(!numcomponents)
		{
			if(domain->ownrank() == 0) cout << endl << "No components defined! Setting up single one-centered LJ" << endl;
			numcomponents=1;
			dcomponents.resize(numcomponents);
			dcomponents[0].setID(0);
			dcomponents[0].addLJcenter(0.,0.,0.,1.,1.,1.,6.,false);
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
			if((x<0.0 || x>=domain->getGlobalLength(0) || y<0.0 || y>=domain->getGlobalLength(1) || z<0.0 || z>=domain->getGlobalLength(2)) && (domain->ownrank() == 0)) {
				cerr << endl << id << ": Molecule " << x << ";" << y << ";" << z << " out of box! " << flush;
			}

			if(((int)componentid > (int)numcomponents) && (domain->ownrank() == 0)) {
				cerr << "Molecule id " << id << " has wrong componentid: " << componentid << ">" << numcomponents << endl;
			}
			--componentid;

			//  store only those molecules within the domain of this process

			// @todo Pointer!!! new!!!  
			Molecule m1 = Molecule(id,componentid,x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz,&dcomponents);
			particleContainer->addParticle(m1);
			dcomponents[componentid].incNumMolecules();
			domain->setglobalRotDOF(dcomponents[componentid].getRotationalDegreesOfFreedom() + domain->getglobalRotDOF());
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

			if(!(i%4096)) if(domain->ownrank() == 0) cout << '.' << flush;
		}
		if(domain->ownrank() == 0) cout << " done" << endl;

		if(!domain->getglobalRho()){
			domain->setglobalRho(domain->getglobalNumMolecules()/(domain->getGlobalLength(0)*domain->getGlobalLength(1)*domain->getGlobalLength(2)));
			if(domain->ownrank() == 0) cout << "calculated global Rho:\t" << domain->getglobalRho() << endl;
		}
#ifndef ENABLE_MPI
		_phaseSpaceFileStream.close();
#endif
	}
	else {
		cerr << "ERROR in AsciiReader::readPhaseSpace: Error in the PhaseSpace File" << endl;
	}

	return maxid;
}
