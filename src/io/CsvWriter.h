#ifndef SRC_IO_CSVWRITER_H_
#define SRC_IO_CSVWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

/** @brief Writes out an ASCII file in the *.csv-format containing coordinates of each molecule.
 *
 * The Plugin writes out a file in the csv file format for easy post-processing or visualization with ParaView.
 * Values are separated by ",".
 */

class CsvWriter : public PluginBase {
public:
	CsvWriter() = default;
	~CsvWriter() override = default;

	/** @brief Read in XML configuration for DecompWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="CsvWriter">
		 <writefrequency>INTEGER</writefrequency>
		 <outputprefix>STRING</outputprefix>
		 <incremental>0|1</incremental>
		 <appendTimestamp>0|1</appendTimestamp>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits &xmlconfig) override;

	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

	std::string getPluginName() override { return std::string("CsvWriter"); }

	static PluginBase *createInstance() { return new CsvWriter(); }

private:
	std::string _outputPrefix;
	unsigned long _writeFrequency{0ul};
	bool _appendTimestamp{false};
	bool _incremental{false};
	const char *const _separator = ",";
};

#endif  // SRC_IO_CSVWRITER_H_
