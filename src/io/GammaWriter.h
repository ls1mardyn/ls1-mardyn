#ifndef SRC_IO_GAMMAWRITER_H_
#define SRC_IO_GAMMAWRITER_H_

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <fstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

/** @brief The GammaWriter plugin writes the surface tension to a file.
 *
 * @todo What is the actual surface? Implementation is in the Domain.
 */
class XMLfileUnits;
class GammaWriter : public OutputBase {
public:
	GammaWriter() {}
	~GammaWriter() {}

	/** @brief Read in XML configuration for GammaWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="GammaWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	//! @todo comment
	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	//! @todo comment
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu,
			std::map<unsigned, CavityEnsemble>* mcav
	);
	//! @todo comment
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("GammaWriter");
	}
	static OutputBase* createInstance() { return new GammaWriter(); }

private:
	void calculateGamma(ParticleContainer* particleContainer, DomainDecompBase* domainDecom);
	double getGamma(unsigned id, double globalLength[3]);
	void resetGamma();

	//! prefix for the names of all output files
	std::ofstream _gammaStream;
	unsigned long _writeFrequency;
	std::string _outputPrefix;
	std::map<unsigned,double> _Gamma; //!< Surface tension component wise
};

#endif  // SRC_IO_GAMMAWRITER_H_
