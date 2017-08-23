#ifndef SRC_IO_POVWRITER_H_
#define SRC_IO_POVWRITER_H_


#include "io/OutputBase.h"


class PovWriter : public OutputBase {
public:
	PovWriter() {}
	~PovWriter() {}

	/** @brief Read in XML configuration for PovWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="PovWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
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
		return std::string("PovWriter");
	}
	static OutputBase* createInstance() { return new PovWriter(); }

private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool  _incremental;
	bool  _appendTimestamp;
};

#endif  // SRC_IO_POVWRITER_H_
