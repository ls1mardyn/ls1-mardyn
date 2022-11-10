#ifndef SRC_IO_GAMMAWRITER_H_
#define SRC_IO_GAMMAWRITER_H_

#include "plugins/PluginBase.h"
#include <string>
#include <fstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

/** @brief The GammaWriter plugin writes the surface tension to a file.
 *
 * @todo What is the actual surface? y-plane?
 */
class XMLfileUnits;
class GammaWriter : public PluginBase {
public:
	GammaWriter() : _gammaStream(), _writeFrequency(1), _outputPrefix("mardyn"), _Gamma() {}
	~GammaWriter() {}

	/** @brief Read in XML configuration for GammaWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="GammaWriter">
	     <writefrequency>INTEGER</writefrequency>  <!-- frequency of writing output file (default: 5000) -->
	     <outputprefix>STRING</outputprefix>       <!-- prefix of output file (default: "gamma") -->
		 <numInterfaces>INT</numInterfaces>        <!-- number of interfaces within sampling range (default: 2) -->
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
	//! @todo comment
	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);
	//! @todo comment
	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    );
	//! @todo comment
	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);
	
	std::string getPluginName() {
		return std::string("GammaWriter");
	}
	static PluginBase* createInstance() { return new GammaWriter(); }

private:
	void calculateGamma(ParticleContainer* particleContainer, DomainDecompBase* domainDecom);
	double getGamma(unsigned id, double globalLength[3]);
	void resetGamma();

	std::ofstream _gammaStream;
	unsigned long _writeFrequency {5000ul};
	std::string _outputPrefix {"gamma"};     // prefix the output file
	unsigned short _numInterfaces {2};     // prefix the output file
	std::map<unsigned,double> _Gamma;        // component-wise surface tension
};

#endif  // SRC_IO_GAMMAWRITER_H_
