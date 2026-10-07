//Calculation of the Fluid Wall interaction by a function

#ifndef MIRROR_SYSTEM_H_
#define MIRROR_SYSTEM_H_

#include "PluginBase.h"

#include <string>
#include <map>
#include <list>
#include <cstdint>

enum MirrorSystemType : uint16_t {
	MST_UNKNOWN = 0,
	MST_SHIFT = 1,
	MST_ENLARGE = 2,
	MST_MIRROR = 3
};

class ParticleContainer;
class DomainDecompBase;
class Domain;

class MirrorSystem : public PluginBase
{
public:
	// constructor and destructor
	MirrorSystem();
	~MirrorSystem();

	/** @brief Read in XML configuration for MirrorSystem and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="MirrorSystem">
			<yPos> <float> </yPos>                     <!-- mirror position -->
			<forceConstant> <float> </forceConstant>   <!-- strength of redirection -->
			<direction> <int> </direction>             <!-- 0|1 , i.e. left |<-- or right -->| mirror -->
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

    /** @brief Method will be called first thing in a new timestep. */
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep, bool signalled
	) override;

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

	std::string getPluginName() override {return std::string("MirrorSystem");}
	static PluginBase* createInstance() {return new MirrorSystem();}

private:
	bool _bDone;
	MirrorSystemType _type;
	double _yPos;
	std::array<double,3> _box_old;
	std::array<double,3> _box_new;
};

#endif /*MIRROR_SYSTEM_H_*/
