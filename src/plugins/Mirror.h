//Calculation of the Fluid Wall interaction by a function

#ifndef MIRROR_H_
#define MIRROR_H_

#include "PluginBase.h"
#include "utils/Random.h"
#include "utils/ObserverBase.h"
#include "utils/Region.h"

#include <string>
#include <map>
#include <list>
#include <cstdint>
#include <vector>
#include <memory>
#include <utility>

#include "utils/CommVar.h"

enum MirrorDirection : uint16_t {
	MD_LEFT_MIRROR = 0,
	MD_RIGHT_MIRROR = 1
};

enum MirrorType : uint16_t {
	MT_UNKNOWN = 0,
    MT_REFLECT = 1,
    MT_FORCE_CONSTANT = 2,
	MT_ZERO_GRADIENT = 3,
	MT_NORMDISTR_MB = 4,
	MT_MELAND_2004 = 5,  // Algorithm proposed by Meland et al., Phys. Fluids, Vol. 16, No. 2 (2004)
	MT_RAMPING = 6,
};

class ParticleContainer;
class DomainDecompBase;
class Domain;

class Mirror : public PluginBase, public ObserverBase, public ControlInstance
{
public:
	// constructor and destructor
	Mirror();
	~Mirror() override = default;

	/** @brief Read in XML configuration for Mirror and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="Mirror" type="5" dir="o-|">   <!-- Mirror type and direction, dir="o-|" or dir="|-o" reflecting particles to the left or right side -->
			<pluginID>INT</pluginID>   <!-- plugin id to enable communication with other plugins -->
			<cid>INT</cid>             <!-- only apply mirror to specified components; 0: all (Default); 1: component 1; etc. -->
			<position>
				<refID>INT</refID>     <!-- coordinate relative to reference point, 1:left interface | 2:right interface
				<coord>FLOAT</coord>   <!-- coordinate of Mirror position -->
			</position>
			<forceConstant>0.</forceConstant>   <!-- force added to particles in order to reflect them from Mirror plane -->
			<meland>
				<use_probability>INT</use_probability>   <!-- 0:disable | 1:enable probability factor in case of Mirror type MT_MELAND_2004 -->
				<velo_target>FLOAT</velo_target>         <!-- target hydrodynamic velocity -->
				<fixed_probability>FLOAT</fixed_probability>         <!-- (optional) fixed probability for reflection in Meland2004 mirror -->
			</meland>
			<diffuse>   <!-- particles will not sharply be reflected at mirror position x, but at x + dx, where dx is a random number between 0 and <width> -->
				<width>FLOAT</width>   <!-- width of region behind mirror in which particles will be reflected -->
			</diffuse>
			<ramping>
				<start>UNSIGNED_LONG</start>   <!-- Timestep until all particles are reflected -->
				<stop>UNSIGNED_LONG</stop>         <!-- As from this timestep all particles are treated as set in "treatment" -->
				<treatment>INT</treatment>         <!-- When not reflected, the particles are deleted (0) or transmitted (1) -->
			</ramping> 
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

	// Getters, Setters
	uint64_t getReflectedParticlesCountLocal(const uint16_t& componentid){return _particleManipCount.reflected.local.at(componentid);}
	uint64_t getDeletedParticlesCountLocal(const uint16_t& componentid) {return _particleManipCount.deleted.local.at(componentid);}
	uint32_t getPluginID() {return _pluginID;}
	void setPluginID(const uint32_t& id) {_pluginID = id;}
	double getPosition() {return _position.coord;}

	// Observer, ControlInstance
	SubjectBase* getSubject();
	void update(SubjectBase* subject) override;
	std::string getShortName() override {return "Mirr";}

private:
		void VelocityChange(ParticleContainer* particleContainer);
		void readNormDistr();

private:
	uint32_t _pluginID;
	uint32_t _targetComp;
	struct MirrorPosition {
		uint16_t axis;
		double coord;
		struct RefPoint {
			uint16_t id;
			double origin;
			double coord;
		} ref;
	} _position;
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

	std::unique_ptr<Random> _rnd;

	struct MelandParams {
		bool use_probability_factor {true};
		double velo_target {0.4};
		float fixed_probability_factor {-1};
	} _melandParams;
	
	struct RampingParams {
		unsigned long startStep {1000};
		unsigned long stopStep {2000};
		int treatment {1};
	} _rampingParams;

	struct ParticleManipCount {
		CommVar<std::vector<uint64_t> > reflected;
		CommVar<std::vector<uint64_t> > deleted;
	} _particleManipCount;
	
	struct DiffuseMirror {
		bool enabled;
		float width;
		std::map<uint64_t,double> pos_map;
	} _diffuse_mirror;
};

#endif /*MIRROR_H_*/
