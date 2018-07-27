#ifndef MAXWRITER_H_
#define MAXWRITER_H_

#include <fstream>
#include <string>
#include <cstdint>
#include <limits>

#include "io/OutputBase.h"

#define FORMAT_SCI_MAX_DIGITS std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)

/** @brief Writes max values of velocity angular momentum and force to a file.
 *
 */
class MaxWriter : public OutputBase
{
public:
	MaxWriter();
	~MaxWriter();

	/** @brief Read in XML configuration for MaxWriter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

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
		return std::string("MaxWriter");
	}
	static OutputBase* createInstance() { return new MaxWriter(); }

private:
	void initDataStructures();
	void doSampling(ParticleContainer* particleContainer);
	void calculateGlobalValues();
	void resetLocalValues();
	void writeData(DomainDecompBase* domainDecomp);

private:
	uint64_t _writeFrequency;
	std::string _outputPrefix;
	double* _dMaxValuesLocal;
	double* _dMaxValuesGlobal;
	uint32_t _numQuantities;
	uint32_t _numValsPerQuantity;
	uint32_t _numValsPerComponent;
	uint32_t _numComponents;
	uint32_t _numVals;
};

#endif /*MAXWRITER_H_*/
