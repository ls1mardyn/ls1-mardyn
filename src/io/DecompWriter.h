#ifndef SRC_IO_DECOMPWRITER_H_
#define SRC_IO_DECOMPWRITER_H_

#include <string>

#include "plugins/PluginBase.h"


/** @brief writes out information about decomposition of the simulation domain.
 *
 * Writes out decomposition information. The data written to the file depend
 * on the used domain decomposition.
 */
class DecompWriter : public PluginBase {
public:
	DecompWriter();
	~DecompWriter() {}


	/** @brief Read in XML configuration for DecompWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="DecompWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <incremental>INTEGER</incremental>
	     <appendTimestamp>INTEGER</appendTimestamp>
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
		return std::string("DecompWriter");
	}
	static PluginBase* createInstance() { return new DecompWriter(); }
private:
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	bool _incremental;
	std::string _outputPrefix;
};

#endif  // SRC_IO_DECOMPWRITER_H_
