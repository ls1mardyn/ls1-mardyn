//Calculation of the Fluid Wall interaction by a function

#ifndef MIRROR_H_
#define MIRROR_H_

#include "PluginBase.h"

#include <string>
#include <map>
#include <list>
#include <cstdint>

enum MirrorDirection : uint16_t {
	MD_LEFT_MIRROR = 0,
	MD_RIGHT_MIRROR = 1
};

enum MirrorType : uint16_t {
	MT_UNKNOWN = 0,
    MT_REFLECT = 1,
    MT_FORCE_CONSTANT = 2,
	MT_ZERO_GRADIENT = 3,
	MT_NORMDISTR_MB = 4
};

class ParticleContainer;
class DomainDecompBase;
class Domain;

class Mirror : public PluginBase
{
public:
	// constructor and destructor
	Mirror();
	~Mirror();

	/** @brief Read in XML configuration for Mirror and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="Mirror">
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

	std::string getPluginName() override {return std::string("Mirror");}
	static PluginBase* createInstance() {return new Mirror();}

private:
		void VelocityChange(ParticleContainer* particleContainer);
		void readNormDistr();

private:
	double _yPos;
	double _forceConstant;
	MirrorDirection _direction;
	MirrorType _type;
	struct NormMB{
		struct NormFnames{
			std::string vxz;
			std::string vy;
		} fname;
		std::list<double> vxz;
		std::list<double> vy;
	} _norm;

	/** ratio of particles to reflect */
	float _ratio;

	/** zero gradient BC */
	struct ComponentIDs {
		uint32_t original;
		uint32_t forward;
		uint32_t backward;
		uint32_t reflected;
		uint32_t permitted;
	} _cids;  // unity based

	struct ControlVolume {
		double left;
		double right;
		double left_outer;
		double right_outer;
		double width;
		double margin;
	} _cv;

	struct VelocityList {
		std::array<double, 3> initvals;
		uint32_t numvals;
		std::list<std::array<double, 3> > list;
	} _veloList;
};

#endif /*MIRROR_H_*/
