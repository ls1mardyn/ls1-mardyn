#include <iostream>
#include <string>
#include "utils/xmlfileUnits.h"
#include "molecules/Component.h"
#include "integrators/Leapfrog.h"

#include "io/xmlreader.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;


XMLReader::XMLReader(const std::string& filename) : _inp(filename) {
	ROOT = "/mardyn";
	if( _inp.changecurrentnode( ROOT ) <= 0 ) {
		global_log->error() << filename << " invalid MarDyn XML file!" << endl;
	}
}



std::string XMLReader::getVersion() {
	std::string version("unknown");
	_inp.getNodeValue("@version", version);
	return version;
}


bool XMLReader::getIntegrator( Integrator* &integrator ) {
	if( _inp.changecurrentnode( ROOT + "/simulation/integrator" ) ) {
		double timestepLength = 0;
		_inp.getNodeValueReduced( "timestep", timestepLength );
		
		std::string integratorType;
		_inp.getNodeValue( "@type", integratorType );
		if( integratorType == "Leapfrog" ) {
			integrator = new Leapfrog(timestepLength);
			return true;
		} 
		else if ( integratorType == "" ) {
			global_log->error() << "No integrator type specified." << endl;
		} else {
			global_log->error() << "Unknown integrator type '" << integratorType << "'." << endl;
		}
		
		global_log->info() << "Integrator: " << integratorType << endl;
		global_log->info() << "timestep: " << timestepLength << endl;
	}
	return false;
}


double XMLReader::getTimestepLength() {
	/* TODO: Do we need this interface? We should go through the integrator. */
	double timestepLength = 0;
	_inp.getNodeValueReduced( ROOT + "/simulation/integrator/timestep", timestepLength );
	return timestepLength;
}


bool XMLReader::getSimBoxSize( double simBoxLength[3]) {
	bool isBox = false;
	if (_inp.changecurrentnode( ROOT + "/simulation/ensemble/volume" )) {
		std::string type;
		_inp.getNodeValue( "@type", type );
		if( type == "box" ) {
			isBox = true;
			_inp.getNodeValueReduced( "lx", simBoxLength[0] );
			_inp.getNodeValueReduced( "ly", simBoxLength[1] );
			_inp.getNodeValueReduced( "lz", simBoxLength[2] );
		}
		_inp.changecurrentnode( ROOT );
	}
	return isBox;
}


long XMLReader::getComponents( std::vector<Component>& components ) {
	long numComponents = 0;
	XMLfile::Query query = _inp.query( ROOT + "/simulation/ensemble/components/moleculetype" );
	numComponents = query.card();
	XMLfile::Query::const_iterator componentIter;
	unsigned int i = 0;
	for( componentIter = query.begin(); componentIter; componentIter++, i++ ) {
		_inp.changecurrentnode( componentIter );
		std::cerr << "PATH: " << _inp.getcurrentnodepath() << endl;
		components.push_back( readComponent( i ) );
	}
	return numComponents;
}


Component XMLReader::readComponent( unsigned int i ) {
	Component component(i);

	std::string cid("");
	_inp.getNodeValue( "@id", cid );
	global_log->info() << "Adding component " << cid << " with id " << i << endl;
// 			std::cerr << "PATH: " << _inp.getcurrentnodepath() << endl;

	/* TODO: read in sites */
	XMLfile::Query query = _inp.query( "site" );
	XMLfile::Query::const_iterator siteIter;
	for( siteIter = query.begin(); siteIter; siteIter++ ) {
		_inp.changecurrentnode( siteIter );
		
		std::string siteType;
		double x, y, z, m;
		x = y = z = m = 0;
		
		_inp.getNodeValue( "@type", siteType );
		_inp.getNodeValueReduced( "coord/x", x );
		_inp.getNodeValueReduced( "coord/y", y );
		_inp.getNodeValueReduced( "coord/z", z );
		_inp.getNodeValueReduced( "mass", m );
		std::cout << "Adding new site of type " << siteType
			<< " (" << x << ", " << y << ", " << z << ")" << endl;
		
		if ( siteType == "LJ126" ) {
			double epsilon = 0;
			double sigma = 0;
			_inp.getNodeValueReduced( "epsilon", epsilon );
			_inp.getNodeValueReduced( "sigma", sigma );
			component.addLJcenter( x, y, z, m, epsilon, sigma );
		} else
		if ( siteType == "Charge" ) {
			double q = 0;
			_inp.getNodeValueReduced( "charge", q );
			component.addCharge( x, y, z, m, q );
		} else
		if ( siteType == "Dipole" ) {
			// TODO: Units?
			// double eMyx, eMyy, eMyz, eMyabs;
			double p[3] = {0, 0, 0};
			double pAbs = 0;
			_inp.getNodeValueReduced( "dipolemoment/x", p[0] );
			_inp.getNodeValueReduced( "dipolemoment/y", p[1] );
			_inp.getNodeValueReduced( "dipolemoment/z", p[2] );
			_inp.getNodeValueReduced( "dipolemoment/abs", pAbs );
			/* In case of unsqecified absolute value assume an absolute vector. */
			if( pAbs == 0 ) {
				pAbs = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
				for( int d = 0; d < 3; d++) {
					p[d] /= pAbs;
				}
			}
			component.addDipole( x, y, z, p[0], p[1], p[2], pAbs );
		} else
		if ( siteType == "Quadrupole" ) {
			// TODO: Units?
			// double eQx, eQy, eQz, eQabs
			double q[3] = {0, 0, 0};
			double qAbs = 0;	
			_inp.getNodeValueReduced( "quadrupolemoment/x", q[0] );
			_inp.getNodeValueReduced( "quadrupolemoment/y", q[1] );
			_inp.getNodeValueReduced( "quadrupolemoment/z", q[2] );
			_inp.getNodeValueReduced( "quadrupolemoment/abs", qAbs );
			/* In case of unsqecified absolute value assume an absolute vector. */
			if( qAbs == 0 ) {
				qAbs = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] );
				for( int d = 0; d < 3; d++) {
					q[d] /= qAbs;
				}
			}
			component.addQuadrupole( x, y, z, q[0], q[1], q[2], qAbs );
		} else
		if ( siteType == "Tersoff" ) {
			// TODO 
// 					double A, B, lambda, my, R, S, c, d, h, n, beta;
// 					A = B = lambda = my = R = S = c = d = h = n = beta = 0;
// 					component.addTersoff( x, y, z, m, A, B, lambda, my, R, S, c, d, h, n, beta );
		}
	}
	
	return component;
}


DomainDecompType XMLReader::getDomainDecompositionType() {
	if ( _inp.changecurrentnode( ROOT + "/simulation/algorithm/parallelization" ) ) {
		std::string domainDecompositionType("DummyDecomposition");
		_inp.getNodeValue( "@type", domainDecompositionType );
		global_log->debug() << "Domain decomposition type is " << domainDecompositionType << endl;
		
		if ( domainDecompositionType == "DummyDecomposition" ) { return DUMMY_DECOMPOSITION; }
#ifdef ENABLE_MPI
		else if ( domainDecompositionType == "DomainDecomposition" ) { return DOMAIN_DECOMPOSITION; }
		else if ( domainDecompositionType == "KDDecomposition" ) { return KD_DECOMPOSITION; }
#endif
		else {
			global_log->error() << "Unknown domain decomposition type '" << domainDecompositionType << "'." << endl;
			return UNKNOWN_DECOMPOSITION;
		}
	} 
	
	global_log->warning() << "No parallelization section found. Using dummy domain decomposition." << endl;
	return DUMMY_DECOMPOSITION;
}

double XMLReader::getTemperature() {
	double temperature = 0;
	_inp.getNodeValue( ROOT + "/simulation/ensemble/temperature", temperature);
	return 	temperature;

}
