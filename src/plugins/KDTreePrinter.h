#include "PluginBase.h"

/**
 * Class to print the kd tree.
 */
class KDTreePrinter : public PluginBase {
public:
	/** @brief Read in XML configuration for KDTreePrinter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="KDTreePrinter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <incremental>INTEGER</incremental>
	     <appendTimestamp>INTEGER</appendTimestamp>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits &xmlconfig) override;

	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

	std::string getPluginName() override;

	static PluginBase* createInstance() { return new KDTreePrinter(); }

private:
	std::string _outputPrefix {"kd-tree"};
	unsigned long _writeFrequency {1000};
	bool	_incremental {true};
	bool	_appendTimestamp {false};

};
