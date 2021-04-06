#ifndef SRC_IO_XYZWRITER_H_
#define SRC_IO_XYZWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

/** @brief Writes out an ASCII file in the *.xyz-format containing coordinates of each molecule.
 *
 * The Plugin writes out a file in the xyz file format. The specifications of the xyz format
 * can be found under  http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html.
 */

class XyzWriter : public PluginBase {
public:
	XyzWriter() = default;
	~XyzWriter() = default;

	/** @brief Read in XML configuration for DecompWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="XyzWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <incremental>INTEGER</incremental>
	     <appendTimestamp>INTEGER</appendTimestamp>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);

	void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    );

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);
	
	std::string getPluginName() {
		return std::string("XyzWriter");
	}
	static PluginBase* createInstance() { return new XyzWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency {0ul};
	bool _appendTimestamp {false};
	bool _incremental {false};
};

#endif  // SRC_IO_XYZWRITER_H_
