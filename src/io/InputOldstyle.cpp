#include "io/InputOldstyle.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"
#include "ensemble/PressureGradient.h"
#include "utils/Logger.h"
#include "utils/Timer.h"
#include <climits>

using Log::global_log;
using namespace std;


InputOldstyle::InputOldstyle() {
}

InputOldstyle::~InputOldstyle(){}

void InputOldstyle::setPhaseSpaceFile(string filename) {
	_phaseSpaceFile = filename;
	_phaseSpaceFileStream.open( filename.c_str() );

	if (!_phaseSpaceFileStream.is_open()) {
		global_log->error() << " Reader for old-style input file: " << endl;
		global_log->error() << "Could not open phaseSpaceFile " << _phaseSpaceFile << endl;
		exit(1);
	}
}

void InputOldstyle::setPhaseSpaceHeaderFile(string filename) {
	_phaseSpaceHeaderFile = filename;
	_phaseSpaceHeaderFileStream.open( filename.c_str() );
}

void InputOldstyle::readPhaseSpaceHeader(Domain* domain, double timestep)
{
	string token, token2;

	_phaseSpaceHeaderFileStream >> token;
	domain->setinpversion(0);

	if(        
			(token != "mardyn") 
			&& (token != "MOLDY") 
			&& (token != "ls1r1") 
			&& (token != "mrdyn") 
			&& (token != "MDProject") 
			&& (token != "Mardyn") 
			&& (token != "MARDYN") 
	  ) 
	{
		global_log->error() << _phaseSpaceHeaderFile << " not a valid OldStyle LS1 input file." << endl;
		exit(1);
	}

	string inputversion;
	_phaseSpaceHeaderFileStream >> token >> inputversion;
	// FIXME: remove tag trunk from file specification?
	if(token != "trunk") {
		global_log->error() << "Wrong input file specifier (\'" << token << "\' instead of \'trunk\')." << endl;
		exit(1);
	}

	domain->setinpversion( strtoul(inputversion.c_str(), NULL, 0) );
	if( domain->getinpversion() < 20080701 ) {
		global_log->error() << "Input version tool old (" << domain->getinpversion() << ")" << endl;
		exit(1);
	}

	global_log->info() << "Reading phase space header from file " << _phaseSpaceHeaderFile << endl;


	vector<Component>& dcomponents = domain->getComponents();
	bool header = true; // When the last header element is reached, "header" is set to false

	while(header) {
		char c;
		_phaseSpaceHeaderFileStream >> c;
		if(c == '#') {
			// comment line
			_phaseSpaceHeaderFileStream.ignore( INT_MAX,'\n' );
			continue;
		}
		_phaseSpaceHeaderFileStream.putback(c);

		token.clear();
		_phaseSpaceHeaderFileStream >> token;
		global_log->debug() << "[[" << token << "]]" << endl;

		if((token == "currentTime") || (token == "t")) {
			// set current simulation time
			_phaseSpaceHeaderFileStream >> token;
			domain->setCurrentTime( strtod(token.c_str(), NULL) );
		}
		else if((token == "Temperature") || (token == "T")) {
			// set global thermostat temperature
			domain->disableComponentwiseThermostat(); //disable component wise thermostats
			double targetT;
			_phaseSpaceHeaderFileStream >> targetT;
			domain->setGlobalTemperature( targetT );
		}
		else if((token == "ThermostatTemperature") || (token == "ThT") || (token == "h")) {
			// set up a new thermostat
			int thermostat_id;
			double targetT;
			_phaseSpaceHeaderFileStream >> thermostat_id;
			_phaseSpaceHeaderFileStream >> targetT;
			domain->setTargetTemperature( thermostat_id, targetT );
		}
		else if((token == "ComponentThermostat") || (token == "CT") || (token == "o")) {
			// specify a thermostat for a component
			if( !domain->severalThermostats() )
				domain->enableComponentwiseThermostat();
			int component_id;
			int thermostat_id;
			_phaseSpaceHeaderFileStream >> component_id >> thermostat_id;
			component_id--; // FIXME thermostat IDs start with 0 in the program but not in the config file?!
			if( thermostat_id < 0 ) // thermostat IDs start with 0
				continue;
			domain->setComponentThermostat( component_id, thermostat_id );
		}
		else if((token == "Undirected") || (token == "U")) {
			// set undirected thermostat
			int thermostat_id;
			_phaseSpaceHeaderFileStream >> thermostat_id;
			domain->enableUndirectedThermostat( thermostat_id );
		}
		else if((token == "Length") || (token == "L")) {
			// simulation box dimensions
			double globalLength[3];
			_phaseSpaceHeaderFileStream >> globalLength[0] >> globalLength[1] >> globalLength[2];
			for( int d = 0; d < 3; d++ )
				domain->setGlobalLength( d, globalLength[d] );
		}
		else if((token == "HeatCapacity") || (token == "cv") || (token == "I"))
		{
			unsigned N;
			double U, UU;
			_phaseSpaceFileStream >> N >> U >> UU;
			domain->init_cv(N, U, UU);
		}
		else if((token == "NumberOfComponents") || (token == "C")) {
			// read in component definitions and
			// read in mixing coefficients

			// components:
			unsigned int numcomponents = 0;
			_phaseSpaceHeaderFileStream >> numcomponents;
			global_log->debug() << "Reading " << numcomponents << " components" << endl;
			dcomponents.resize(numcomponents);
			for( unsigned int i = 0; i < numcomponents; i++ ) {
				dcomponents[i].setID(i);
				unsigned int numljcenters = 0;
				unsigned int numcharges = 0;
				unsigned int numdipoles = 0;
				unsigned int numquadrupoles = 0;
				unsigned int numtersoff = 0;
				_phaseSpaceHeaderFileStream >> numljcenters >> numcharges >> numdipoles 
					>> numquadrupoles >> numtersoff;

				double x, y, z, m;
				for( unsigned int j = 0; j < numljcenters; j++ ) {
					double eps, sigma, tcutoff, do_shift;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> m >> eps >> sigma >> tcutoff >> do_shift;
					dcomponents[i].addLJcenter( x, y, z, m, eps, sigma, tcutoff, (do_shift != 0) );
				}
				for( unsigned int j = 0; j < numcharges; j++ ) {
					double q;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> m >> q;
					dcomponents[i].addCharge( x, y, z, m, q );
				}
				for( unsigned int j = 0; j < numdipoles; j++ ) {
					double eMyx,eMyy,eMyz,absMy;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> eMyx >> eMyy >> eMyz >> absMy;
					dcomponents[i].addDipole( x, y, z, eMyx, eMyy, eMyz, absMy );
				}
				for( unsigned int j = 0; j < numquadrupoles; j++ ) {
					double eQx,eQy,eQz,absQ;
					_phaseSpaceHeaderFileStream >> x >> y >> z >> eQx >> eQy >> eQz >> absQ;
					dcomponents[i].addQuadrupole(x,y,z,eQx,eQy,eQz,absQ);
				}
				for( unsigned int j = 0; j < numtersoff; j++ ) {
					double x, y, z, m, A, B, lambda, mu, R, S, c, d, h, n, beta;
					_phaseSpaceHeaderFileStream >> x >> y >> z;
					_phaseSpaceHeaderFileStream >> m >> A >> B;
					_phaseSpaceHeaderFileStream >> lambda >> mu >> R >> S;
					_phaseSpaceHeaderFileStream >> c >> d >> h >> n >> beta;
					dcomponents[i].addTersoff( x, y, z, m, A, B, lambda, mu, R, S, c, d, h, n, beta );
				}
				double IDummy1,IDummy2,IDummy3;
				// FIXME! Was soll das hier? Was ist mit der Initialisierung im Fall I <= 0.
				_phaseSpaceHeaderFileStream >> IDummy1 >> IDummy2 >> IDummy3;
				if( IDummy1 > 0. ) dcomponents[i].setI11(IDummy1);
				if( IDummy2 > 0. ) dcomponents[i].setI22(IDummy2);
				if( IDummy3 > 0. ) dcomponents[i].setI33(IDummy3);
			}

#ifndef NDEBUG
			for (int i = 0; i < numcomponents; i++) {
				global_log->debug() << "Component " << (i+1) << " of " << numcomponents << endl;
				global_log->debug() << dcomponents[i] << endl;
			}
#endif

			// Mixing coefficients
			vector<double>& dmixcoeff = domain->getmixcoeff();
			dmixcoeff.clear();
			for( unsigned int i = 1; i < numcomponents; i++ ) {
				for( unsigned int j = i + 1; j <= numcomponents; j++ ) {
					double xi, eta;
					_phaseSpaceHeaderFileStream >> xi >> eta;
					dmixcoeff.push_back( xi );
					dmixcoeff.push_back( eta );
				}
			}
			// read in global factor \epsilon_{RF}
			// FIXME: Maybe this should go better to a seperate token?!
			_phaseSpaceHeaderFileStream >> token;
			domain->setepsilonRF( strtod(token.c_str(),NULL) );
			long int fpos;
			if( _phaseSpaceFile == _phaseSpaceHeaderFile ) {
				// in the case of a single phase space header + phase space file
				// find out the actual position, because the phase space definition will follow
				// FIXME: is there a more elegant way?
				fpos = _phaseSpaceHeaderFileStream.tellg();
				_phaseSpaceFileStream.seekg( fpos, ios::beg );
			}
			// FIXME: Is there a better solution than skipping the rest of the file?
			header = false;
		}
		else if((token == "NumberOfMolecules") || (token == "N")) {
			// set number of Molecules 
			// FIXME: Is this part called in any case as the token is handled in the readPhaseSpace method?
			_phaseSpaceHeaderFileStream >> token;
			domain->setglobalNumMolecules( strtoul(token.c_str(),NULL,0) );
		}
		else if((token == "AssignCoset") || (token == "S")) {
			unsigned component_id, cosetid;
			_phaseSpaceHeaderFileStream >> component_id >> cosetid;
			component_id--; // FIXME: Component ID starting with 0 in program ...
			domain->getPG()->assignCoset( component_id, cosetid );
		}
		else if((token == "Accelerate") || (token == "A")) {
			unsigned cosetid;
			_phaseSpaceHeaderFileStream >> cosetid;
			double v[3];
			for(unsigned d = 0; d < 3; d++) 
				_phaseSpaceHeaderFileStream >> v[d];
			double tau;
			_phaseSpaceHeaderFileStream >> tau;
			double ainit[3];
			for(unsigned d = 0; d < 3; d++) 
				_phaseSpaceHeaderFileStream >> ainit[d];
			domain->getPG()->specifyComponentSet(cosetid, v, tau, ainit, timestep);
		}
		else {
			global_log->error() << "Invalid token \'" << token << "\' found. Skipping rest of the header." << endl;
			header = false;
		}
	}
}

unsigned long InputOldstyle::readPhaseSpace(ParticleContainer* particleContainer, list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {

	Timer inputTimer;
	inputTimer.start();
	global_log->info() << "Reading phase space file " << _phaseSpaceFile << endl;
	
	string token;
	vector<Component>& dcomponents = domain->getComponents();
	unsigned int numcomponents = dcomponents.size();
	unsigned long nummolecules;
	unsigned long maxid = 0; // stores the highest molecule ID found in the phase space file
	string ntypestring("ICRVQD");
	enum Ndatatype { ICRVQD, IRV, ICRV, ICRVFQDM } ntype = ICRVQD;

	_phaseSpaceFileStream >> token;
	if((token != "NumberOfMolecules") && (token != "N")) {
		global_log->error() << "Expected token 'NumberOfMolecules (N)' instead of '" << token << "'" << endl;
		exit(1);
	}
	_phaseSpaceFileStream >> nummolecules;
	global_log->info() << " number of molecules: " << nummolecules << endl;
	domain->setglobalNumMolecules( nummolecules );
	
	streampos spos = _phaseSpaceFileStream.tellg();
	_phaseSpaceFileStream >> token;
	if((token=="MoleculeFormat") || (token == "M")) {

		if( domain->getinpversion() >= 51129 )  // TODO: Is this a magic number?
			_phaseSpaceFileStream >> ntypestring;
		ntypestring.erase( ntypestring.find_last_not_of( " \t\n") + 1 );
		ntypestring.erase( 0, ntypestring.find_first_not_of( " \t\n" ) );
		
		if (ntypestring == "ICRVQD") ntype = ICRVQD;
		else if (ntypestring == "ICRVFQDM")  ntype = ICRVFQDM;
		else if (ntypestring == "ICRV") ntype = ICRV;
		else if (ntypestring == "IRV")  ntype = IRV;
		else {
			global_log->error() << "Unknown molecule format '" << ntypestring << "'" << endl;
			exit(1);
		}
	} else {
		_phaseSpaceFileStream.seekg(spos);
	}
	global_log->info() << " molecule format: " << ntypestring << endl;
	
	if( numcomponents < 1 ) {
		global_log->warning() << "No components defined! Setting up single one-centered LJ" << endl;
		numcomponents = 1;
		dcomponents.resize( numcomponents );
		dcomponents[0].setID(0);
		dcomponents[0].addLJcenter(0., 0., 0., 1., 1., 1., 6., false);
	}


	double x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz;
	double Fx,Fy,Fz,Mx,My,Mz;
	unsigned long id;
	unsigned int componentid;

	x=y=z=vx=vy=vz=q1=q2=q3=Dx=Dy=Dz=0.;
	q0=1.;
	Fx=Fy=Fz=Mx=My=Mz=0.;

	for( unsigned long i = 0; i < domain->getglobalNumMolecules(); i++ ) {

		switch ( ntype ) {
			case ICRVQD:
				_phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
					>> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz;
				break;
			case ICRVFQDM :
				_phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz >> Fx >> Fy >> Fz
					>> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz >> Mx >> My >> Mz;
				break;
			case ICRV :
				_phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz;
				break;
			case IRV :
				_phaseSpaceFileStream >> id >> x >> y >> z >> vx >> vy >> vz;
				break;
		}
		if(        ( x < 0.0 || x >= domain->getGlobalLength(0) )
				|| ( y < 0.0 || y >= domain->getGlobalLength(1) ) 
				|| ( z < 0.0 || z >= domain->getGlobalLength(2) ) ) 
		{
			global_log->warning() << "Molecule " << id << " out of box: " << x << ";" << y << ";" << z << endl;
		}

		if( componentid > numcomponents ) {
			global_log->error() << "Molecule id " << id << " has wrong componentid: " << componentid << ">" << numcomponents << endl;
			exit(1);
		}
		componentid --; // TODO: Component IDs start with 0 in the program.

		// store only those molecules within the domain of this process
		// The neccessary check is performed in the particleContainer addPartice method
		// FIXME: Datastructures? Pass pointer instead of object, so that we do not need to copy?!
		Molecule m1 = Molecule(id,componentid,x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz,&dcomponents);
		particleContainer->addParticle(m1);
		
		 // TODO: The following should be done by the addPartice method.
		dcomponents[componentid].incNumMolecules();
		domain->setglobalRotDOF(dcomponents[componentid].getRotationalDegreesOfFreedom() + domain->getglobalRotDOF());
		
		if(id > maxid) maxid = id;

		std::list<ChemicalPotential>::iterator cpit;
		for(cpit = lmu->begin(); cpit != lmu->end(); cpit++) {
			if( !cpit->hasSample() && (componentid == cpit->getComponentID()) ) {
				cpit->storeMolecule(m1);
			}
		}

		// Print status message
		unsigned long iph = domain->getglobalNumMolecules() / 100;
		if( iph != 0 && (i % iph) == 0 )
			global_log->info() << "Finished reading molecules: " << i/iph << "%\r" << flush;
	}

	global_log->info() << "Finished reading molecules: 100%" << endl;
	global_log->info() << "Reading Molecules done" << endl;

	// TODO: Shouldn't we always calculate this?
	if( !domain->getglobalRho() ){
		domain->setglobalRho( domain->getglobalNumMolecules() / domain->getGlobalVolume() );
		global_log->info() << "Calculated Rho_global = " << domain->getglobalRho() << endl;
	}

	_phaseSpaceFileStream.close();

	inputTimer.stop();
	global_log->info() << "Initial IO took:                 " << inputTimer.get_etime() << " sec" << endl;
	return maxid;
}
