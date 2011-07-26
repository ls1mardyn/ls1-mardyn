#ifndef XMLREADER_H
#define XMLREADER_H

#include <string>

#include "parallel/DomainDecompTypes.h"

class Component;
class Integrator;
class XMLfileUnits;


/** Reader for the MarDyn XML input file format
 * 
 * The XMLReader encapsulates the data format of the MarDyn XML input file 
 * format and provides access to the inout data using a generic interface.
 */
class XMLReader {
	
	private:
  
		std::string ROOT; /**< root node of the XML file */
		XMLfileUnits _inp;
		
	public:
		XMLReader(const std::string& filename); 
		~XMLReader() {}
		
		/** get input file format version */
		std::string getVersion();
		
		/** get timestep length */
		double getTimestepLength();
		
		/** get simulation box length 
		 *  @param[out] simBoxLength array to hold the box dimensions
		 *  @return false if simulation volume is not a box
		 */
		bool getSimBoxSize( double simBoxLength[3]);
		
		/** get components 
		 *  @param[out] components vector to which components will be appended
		 *  @return number of read components
		 */
		long getComponents( std::vector<Component>& components );
		
		/** get and initialize integrator
		 *  @param[out] integrator
		 *  @return false on error
		 */
		bool getIntegrator( Integrator *integrator );
		
		/** get the domain decomposition type
		 *  @return type of the domain decomposition as defined in DomainDecompTypes.h
		 */
		DomainDecompType getDomainDecompositionType();
		
	private:
		/** read in a single component */
		Component readComponent( unsigned int i );
};


#endif /* XMLREADER_H */
