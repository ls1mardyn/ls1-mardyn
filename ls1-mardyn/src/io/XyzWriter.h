#ifndef SRC_IO_XYZWRITER_H_
#define SRC_IO_XYZWRITER_H_

#include <string>

#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"

/** @brief Writes out an ASCII file in the *.xyz-format containing coordinates of each molecule.
 *
 * The Plugin writes out a file in the xyz file format. The specifications of the xyz format
 * can be found under  http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html.
 */

class XyzWriter : public OutputBase {
public:
	XyzWriter() {}
	~XyzWriter() {}

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

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu,
			std::map<unsigned, CavityEnsemble>* mcav
	);

	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("XyzWriter");
	}
	static OutputBase* createInstance() { return new XyzWriter(); }
private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	bool _incremental;
};

#endif  // SRC_IO_XYZWRITER_H_
