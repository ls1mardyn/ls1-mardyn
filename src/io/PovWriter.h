#ifndef SRC_IO_POVWRITER_H_
#define SRC_IO_POVWRITER_H_


#include "plugins/PluginBase.h"


class PovWriter : public PluginBase {
public:
	PovWriter() {}
	~PovWriter() {}

	/** @brief Read in XML configuration for PovWriter and all its included objects.
	 *
	 * The povray output plugin writes povray input files. It is intended for smaller
	 * (serial) runs. For visualisation of larger systems please use the MmpldWriter.
	 *
	 * @note The PovWriter works only in serial execution, it is not parallel!
	 * @todo Implement parallel output
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="PovWriter">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <incremental>BOOL</incremental>
	     <appendTimestamp>BOOL</appendTimestamp>
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
		return std::string("PovWriter");
	}
	static PluginBase* createInstance() { return new PovWriter(); }

private:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool  _incremental;
	bool  _appendTimestamp;
};

#endif  // SRC_IO_POVWRITER_H_
