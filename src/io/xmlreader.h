#ifndef XMLREADER_H
#define XMLREADER_H

#include <string>

class Component;
class XMLfileUnits;

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
		
	private:
		/** read in a single component */
		Component readComponent( unsigned int i );
};


#endif /* XMLREADER_H */
