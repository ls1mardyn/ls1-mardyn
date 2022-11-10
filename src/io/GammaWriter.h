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
    GammaWriter() :
        _gammaStream(),
        _writeFrequency(5000ul),
        _outputPrefix("gamma"),
		_numInterfaces(2),
        _gamma(),
        _range{0,20000}, // Reasonable filling value for ymax as domain length is not available here
		_globalLength{0.0,0.0,0.0},
		_numComp(0)
    {}

    ~GammaWriter() {}

    /** @brief Read in XML configuration for GammaWriter.
     *
     * The following xml object structure is handled by this method:
     * \code{.xml}
        <outputplugin name="GammaWriter">
            <writefrequency>INTEGER</writefrequency>   <!-- Frequency to write out data to file; default: 1000 -->
            <outputprefix>STRING</outputprefix>        <!-- Prefix string of output file name; default: gamma -->
			<numInterfaces>INT</numInterfaces>         <!-- Number of interfaces within sampling range; default: 2 -->
            <range>                                    <!-- Range in which to sample data -->
                <ymin>DOUBLE</ymin>                    <!-- Sampled in y-direction between ymin and ymax; defaults: 0 and domain length -->
                <ymax>DOUBLE</ymax>
            </range>
        </outputplugin>
       \endcode
     */
    void readXML(XMLfileUnits& xmlconfig) override;

    // Initialization of output file
    void init(ParticleContainer */* particleContainer */,
              DomainDecompBase *domainDecomp, Domain *domain) override;

    // Calculation of surface tension and, if applicable, output to file
    void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    ) override;

    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {}
    
    std::string getPluginName() override { return "GammaWriter"; }
    static PluginBase* createInstance() { return new GammaWriter(); }

private:
    void calculateGamma(ParticleContainer *particleContainer, DomainDecompBase *domainDecom);
    inline double getGamma(unsigned id);
    inline void resetGamma(Domain *domain);

    std::ofstream _gammaStream;
	unsigned long _writeFrequency;
	std::string _outputPrefix;      // prefix the output file
	unsigned int _numInterfaces;    // number of interfaces within range
	std::vector<double> _gamma;     // component-wise surface tension; 0 is all components
	// Range within particles are considered for calculation of surface tension
    struct Range {
        double ymin, ymax;
    } _range;

	std::array<double, 3> _globalLength;
	unsigned int _numComp;
};

#endif  // SRC_IO_GAMMAWRITER_H_
