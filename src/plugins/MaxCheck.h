/*
 * MaxCheck.h
 *
 *  Created on: 28.05.2018
 *      Author: mheinen
 */

#ifndef MAXCHECK_H_
#define MAXCHECK_H_

#include "utils/PluginBase.h"

#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <cstdint>

enum MaxCheckMethods
{
	MCM_UNKNOWN = 0,
	MCM_LIMIT_TO_MAX_VALUE = 1,
	MCM_DELETE_PARTICLES = 2
};

struct TimestepControl
{
	uint64_t start;
	uint64_t freq;
	uint64_t stop;
};

struct MaxVals
{
	double F;
	double F2;
	double v;
	double v2;
	uint32_t method;
};
typedef std::unordered_map<uint32_t, MaxVals> maxvals_map;

class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;
class CavityEnsemble;
class MaxCheck : public PluginBase
{
public:
	// constructor and destructor
	MaxCheck();
	~MaxCheck();

	/** @brief Read in XML configuration for MaxCheck and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="MaxCheck">
			<Fmax> <DOUBLE> </Fmax>  <!-- max. allowed force -->
			<method> <INT> <method>  1:inform | 2:limit to max value | 3:delete particle
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

	/** @brief Method afterForces will be called after forcefields have been applied
	 *
	 * make pure Virtual ?
	 */
	void afterForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	void endStep(
			ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, Domain *domain,
			unsigned long simstep, std::list<ChemicalPotential> *lmu,
			std::map<unsigned, CavityEnsemble> *mcav
	) override;

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("MaxCheck");}
	static PluginBase* createInstance() {return new MaxCheck();}

private:
	double calcSquaredVectorLength(double* vec) {return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);}
	void checkMaxVals(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);

private:
	TimestepControl _control;
	maxvals_map _maxVals;
	std::vector<Molecule*> _deletions;
};

#endif /*MAXCHECK_H_*/
