#ifndef TCTS_H_
#define TCTS_H_

#include "io/InputBase.h"

/** @brief Generates a scenario of two layers of the same component but with different densities.
 *
 * Creates two layers of the same height in the simulation box. The density of
 * the two layers can be chosen separately. The layers consist both out of
 * molecules from the component with ID=1 in the xml input file.
 */
class MkTcTSGenerator : public InputBase {

public:
	MkTcTSGenerator() {}

	~MkTcTSGenerator() {}

	void readPhaseSpaceHeader(Domain* domain, double timestep);

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);


	/** @brief Read in XML configuration for MkTcTSGenerator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <generator name="mkTcTS">
	     <layer1>
	       <density>DOUBLE</density>
	     </layer1>
	     <layer2>
	       <density>DOUBLE</density>
	     </layer2>
	   </generator>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

private:
	double l1_ratio;
	double l1_offset;

	double rho1;
	double rho2;

	void _msimulation(int arg1);
};

#endif /* TCTS_H_ */
