//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_SPATIALPROFILE_H
#define MARDYN_TRUNK_SPATIALPROFILE_H

#include <functional>
#include <optional>
#include <vector>

#include <plugins/profiles/ProfileBase.h>
#include "PluginBase.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"


class ProfileBase;

class DensityProfile;

class Velocity3dProfile;

class VelocityAbsProfile;

class KineticProfile;

class DOFProfile;

class TemperatureProfile;

class VirialProfile;

class Virial2DProfile;


/** @brief SpatialProfile is a Plugin that is called like any other plugin derived from PluginBase. It handles all profiles in /plugins/profiles. <br>
 * New profiles must be added via the plugins/ ProfileBase to comply with this Plugin. <br>
 *
 *
 * <b>x y z</b>: Set Sampling Grid for all profiles. Output currently only readable for <b>x = 1</b>. <br>
 * <b>writefrequency</b>: Write frequency for profile output. The accumulated profile data is written out every Nth step and the profiles reset.<br>
 * <b>timesteps/init</b>: Number of timesteps to not record/output profiles at the beginning of the simulation. <br>
 * <b>timesteps/recording</b>: During the recording period between writes, select every Nth step to be recorded for the profile. <br>
 * <b>outputprefix</b>: Profile File output prefix. <br>
 * <b>profiles</b>: Enable/disable the different profiles. If one is needed for another (Density for Velocity) it is enabled automatically. If none are specified here all profiles will be written.<br>
 * \code{.xml}
* <plugin name="SpatialProfile">
 *    <mode>cylinder</mode> (cartesian)
      <x>1</x>
      <y>20</y>
      <z>20</z>
      <r>20</r>
      <h>20</h>
      <phi>1</phi>
      <writefrequency>100</writefrequency>
      <timesteps>
        <init>1</init>
        <recording>1</recording>
      </timesteps>
      <outputprefix>comparison</outputprefix>
      <profiles>
        <density>true</density>
        <temperature>false</temperature>
        <velocity>true</velocity>
        <velocity3d>true</velocity>
        <virial>true</virial>
      </profiles>
    </plugin>
* \endcode
 */
class SpatialProfile : public PluginBase {

public:

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	void readXML(XMLfileUnits& xmlconfig) override;

	void endStep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep) override;

	void finish(ParticleContainer* particleContainer,
				DomainDecompBase* domainDecomp, Domain* domain) override {};

	unsigned long getCartesianUID(ParticleIterator& thismol);

	long getCylUID(ParticleIterator& thismol);

	std::string getPluginName() override { return std::string("SpatialProfile"); }

	static PluginBase* createInstance() { return new SpatialProfile(); }

	SamplingInformation samplInfo;

	void accessAllCallbacks(const std::map<std::string, FunctionWrapper>& callbackMap) override {
		// Accesses a callback registered by FixRegion. It returns the number of molecules in the fixregion.
		std::string name{"FixRegion::getMoleculesInRegion"};
		if(callbackMap.find(name) != callbackMap.end()) {
			getNumFixRegion = callbackMap.at(name).get<unsigned long>();
		}
	}

private:

	// Profile pointers for data reuse
	DensityProfile* _densProfile; //!< Reference to DensityProfile as it is needed by most other profiles
	VelocityAbsProfile* _velAbsProfile;
	Velocity3dProfile* _vel3dProfile;
	TemperatureProfile* _tempProfile;
	DOFProfile* _dofProfile;
	KineticProfile* _kineticProfile;
	VirialProfile* _virialProfile;
	Virial2DProfile* _virial2DProfile;

	unsigned long _writeFrequency; // Write frequency for all profiles -> Length of recording frame before output
	unsigned long _initStatistics; // Timesteps to skip at start of the simulation
	unsigned long _profileRecordingTimesteps; // Record every Nth timestep during recording frame
	long _accumulatedDatasets; // Number of Datasets between output writes / profile resets
	std::string _outputPrefix; // File prefix for all profiles
	std::string _mode;
	std::string _profiledCompString;
	unsigned int _profiledComp;


	unsigned long _uIDs; //!< Total number of unique IDs with the selected Grid. This is the number of total bins in the Sampling grid.

	std::vector<ProfileBase*> _profiles; // vector holding all enabled profiles
	int _comms = 0; // total number of communications per bin needed by all profiles.

	// Needed for XML check for enabled profiles.
	bool _ALL = false;
	bool _DENSITY = false;
	bool _VELOCITY = false;
	bool _VELOCITY3D = false;
	bool _TEMPERATURE = false;
	bool _VIRIAL = false;
	bool _VIRIAL2D = false;

	void addProfile(ProfileBase* profile);

	std::optional<std::function<unsigned long(void)>> getNumFixRegion;

};


#endif //MARDYN_TRUNK_SPATIALPROFILE_H
