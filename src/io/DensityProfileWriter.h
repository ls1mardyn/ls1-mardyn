#ifndef SRC_IO_DENSITYPROFILEWRITER_H_
#define SRC_IO_DENSITYPROFILEWRITER_H_

#include "plugins/PluginBase.h"

/** @brief The DensityProfileWriter obtains a planar interface density profile.
 */
class DensityProfileWriter : public PluginBase {
public:
	DensityProfileWriter() : _writeFrequency(1), _initStatistics(0), _profileRecordingTimesteps(1), _outputPrefix("profile2"), _doRecordVirialProfile(false) {}
	~DensityProfileWriter() {}

	/** @brief Read in XML configuration for DensityProfileWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="PovWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits & xmlconfig);

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);

	void endStep(ParticleContainer *particleContainer,
                 DomainDecompBase *domainDecomp, Domain *domain,
                 unsigned long simstep, std::list<ChemicalPotential> *lmu,
                 std::map<unsigned, CavityEnsemble> *mcav);

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() { return std::string("DensityProfileWriter"); }
	static PluginBase * createInstance() {
		return new DensityProfileWriter();
	}
private:
	unsigned long _writeFrequency;  //!< was Simulation::_profileOutputTimesteps
	unsigned long _initStatistics;  //!< @todo needs documentation
	unsigned long _profileRecordingTimesteps;  //!< @todo needs documentation
	std::string _outputPrefix;  //!< output file prefix
	bool _doRecordVirialProfile;  //!< @todo needs documentation
};

#endif  // SRC_IO_DENSITYPROFILEWRITER_H_
