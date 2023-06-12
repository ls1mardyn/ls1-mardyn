#ifndef MAXWRITER_H_
#define MAXWRITER_H_

#include <fstream>
#include <string>
#include <cstdint>
#include <limits>
#include <vector>

#include "plugins/PluginBase.h"

#define FORMAT_SCI_MAX_DIGITS std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)

/** @brief Writes max values of velocity angular momentum and force to a file.
 *
 */
class MaxWriter : public PluginBase
{
public:
	MaxWriter();
	~MaxWriter();

	/** @brief Read in XML configuration for MaxWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

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
		return std::string("MaxWriter");
	}
	static PluginBase* createInstance() { return new MaxWriter(); }

private:
	void initDataStructures();
	void doSampling(ParticleContainer* particleContainer);
	void calculateGlobalValues(DomainDecompBase *domainDecomp);
	void resetLocalValues();
	void writeData(DomainDecompBase* domainDecomp);

private:
	uint64_t _writeFrequency;
	std::string _outputPrefix;
	std::vector<double> _dMaxValuesLocal;
	std::vector<double> _dMaxValuesGlobal;
	uint32_t _numQuantities;
	uint32_t _numValsPerQuantity;
	uint32_t _numValsPerComponent;
	uint32_t _numComponents;
	uint32_t _numVals;
};

#endif /*MAXWRITER_H_*/
