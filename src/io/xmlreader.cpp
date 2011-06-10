#include <iostream>
#include <string>
#include "utils/xmlfileUnits.h"
#include "molecules/Component.h"

#include "io/xmlreader.h"

XMLReader::XMLReader(const std::string& filename) : _inp(filename) {
			ROOT = "/mardyn";
			_inp.changecurrentnode( ROOT ); /* TODO: checks */
			std::cerr << "PATH: " << _inp.getcurrentnodepath() << std::endl;
		}

		
		std::string XMLReader::getVersion() {
			std::string version("unknown");
			_inp.getNodeValue("@version", version);
			return version;
		}
		
		double XMLReader::getTimestepLength() {
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
				std::cerr << "PATH: " << _inp.getcurrentnodepath() << std::endl;
				components.push_back( readComponent( i ) );
			}
			return numComponents;
		}
		
		Component XMLReader::readComponent( unsigned int i ) {
			Component component(i);

			std::string cid("");
			_inp.getNodeValue( "@id", cid );
			std::cout << "Adding new component " << cid << " with id " << i << std::endl;
// 			std::cerr << "PATH: " << _inp.getcurrentnodepath() << std::endl;

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
				_inp.getNodeValueReduced( "coord/mass", m );
				std::cout << "Adding new site of type " << siteType
					<< " (" << x << ", " << y << ", " << z << ")" << std::endl;
				
				if ( siteType == "LJ126" ) {
					double epsilon = 0;
					double sigma = 0;
					_inp.getNodeValueReduced( "coord/epsilon", epsilon );
					_inp.getNodeValueReduced( "coord/sigma", sigma );
					component.addLJcenter( x, y, z, m, epsilon, sigma );
				} else
				if ( siteType == "Charge" ) {
					double q = 0;
					_inp.getNodeValueReduced( "coord/charge", q );
					component.addCharge( x, y, z, m, q );
				} else
				if ( siteType == "Dipole" ) {
					// TODO
					// double eMyx, eMyy, eMyz, eMyabs;
					// component.addDipole( x, y, z, eMyx, eMyy, eMyz, eMyabs );
				} else
				if ( siteType == "Quadrupole" ) {
					// TODO
				} else
				if ( siteType == "Tersoff" ) {
					// TODO
				}
			}
			
			return component;
		}
