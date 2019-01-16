#ifndef MKESFERA_H
#define MKESFERA_H

#include "io/InputBase.h"

/** @brief Single droplet scenario generator.
 *
 * Creates a droplet with a given center position and radius inside the simulation box.
 * The density of the droplet and its surrounding can be chosen seperately.
 * The droplet and its surrounding consist both out of molecules from the component
 * with ID=1 in the xml input file.
 */
class MkesferaGenerator : public InputBase {

public:
	MkesferaGenerator() {}

	~MkesferaGenerator() {}

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {}

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);

	/** @brief Read in XML configuration for MkesferaGenerator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <generator name="mkesfera">
	     <outer-density>DOUBLE</outer-density>
	     <droplet>
	       <radius>DOUBLE</radius>
	       <density>DOUBLE</density>
	       <center>
	         <x>DOUBLE</x>
	         <y>DOUBLE</y>
	         <z>DOUBLE</z>
	       </center>
	     </droplet>
	   </generator>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

private:
	double R_i, R_o;
	double rho_i, rho_o;
	double center[3]; /**< droplet center */
};

#endif /* MKESFERA_H */

