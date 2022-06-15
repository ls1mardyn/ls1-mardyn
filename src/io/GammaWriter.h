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
        _writeFrequency(1000),
        _outputPrefix("gamma"),
        _Gamma(),
        _range{0,2000} // Resonable filling value for ymax as domain length is not available here
    {}

    ~GammaWriter() {}

    /** @brief Read in XML configuration for GammaWriter.
     *
     * The following xml object structure is handled by this method:
     * \code{.xml}
        <outputplugin name="GammaWriter">
            <writefrequency>INTEGER</writefrequency>   <!-- Frequency to write out data to file; default: 1000 -->
            <outputprefix>STRING</outputprefix>        <!-- Prefix string of output file name; default: gamma -->
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

    // Calculation of surface tension and, if applicable, output to file takes place in endStep()
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
    double getGamma(unsigned id, double globalLength[3]);
    void resetGamma(Domain *domain);

    std::ofstream _gammaStream;
    unsigned long _writeFrequency;
    std::string _outputPrefix;  //!< prefix the output file
    std::map<unsigned,double> _Gamma;  //!< Surface tension component wise; 0 represents all components

    // Range within particles are considered for calculation of surface tension
    struct Range {
        double ymin, ymax;
    } _range;
};

#endif  // SRC_IO_GAMMAWRITER_H_