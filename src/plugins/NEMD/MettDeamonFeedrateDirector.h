//Calculation of the Fluid Wall interaction by a function

#ifndef METTDEAMON_FEEDRATE_DIRECTOR_H_
#define METTDEAMON_FEEDRATE_DIRECTOR_H_

#include "plugins/PluginBase.h"
#include "utils/Random.h"
#include "plugins/NEMD/MettDeamon.h"

#include <string>
#include <map>
#include <list>
#include <cstdint>
#include <vector>

#include "utils/CommVar.h"

class ParticleContainer;
class DomainDecompBase;
class Domain;

class MettDeamonFeedrateDirector : public PluginBase
{
public:
	// constructor and destructor
	MettDeamonFeedrateDirector();
	~MettDeamonFeedrateDirector();

	/** @brief Read in XML configuration for Mirror and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="MettDeamonFeedrateDirector">
			<yPos> <float> </yPos>                     <!-- mirror position -->
			<forceConstant> <float> </forceConstant>   <!-- strength of redirection -->
			<direction> <int> </direction>             <!-- 0|1 , i.e. left |<-- or right -->| mirror -->
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

	void beforeForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

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
			unsigned long simstep) override {}

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("MettDeamonFeedrateDirector");}
	static PluginBase* createInstance() {return new MettDeamonFeedrateDirector();}

private:
	void calcFeedrate(MettDeamon* mettDeamon);
	void resetLocalValues();

private:
	uint32_t _mirror_id;
	struct UpdateControl {
		uint32_t updateFreq;
		uint32_t sampledTimestepCount;
	} _updateControl;
	double _feedrate;
	struct ParticleManipCount {
		CommVar<std::vector<uint64_t> > reflected;
		CommVar<std::vector<uint64_t> > deleted;
	} _particleManipCount;
};

#endif /*METTDEAMON_FEEDRATE_DIRECTOR_H_*/
