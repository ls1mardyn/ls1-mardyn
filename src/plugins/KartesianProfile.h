//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_KARTESIAN2DPROFILE_H
#define MARDYN_TRUNK_KARTESIAN2DPROFILE_H

#include "PluginBase.h"
#include "Domain.h"
#include "plugins/profiles/DensityProfile.h"
#include "plugins/profiles/Velocity3dProfile.h"
#include "plugins/profiles/VelocityAbsProfile.h"

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

/** @brief KartesianProfile is a Plugin that is called like any other plugin derived from PluginBase. It handles all profiles in /plugins/profiles. <br>
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
* <plugin name="KartesianProfile">
      <x>1</x>
      <y>20</y>
      <z>20</z>
      <writefrequency>100</writefrequency>
      <timesteps>
        <init>1</init>
        <recording>1</recording>
      </timesteps>
      <outputprefix>comparison</outputprefix>
      <profiles>
        <density>1</density>
        <temperature>0</temperature>
        <velocity>1</velocity>
      </profiles>
    </plugin>
* \endcode
 */
class KartesianProfile : public PluginBase {

public:
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void readXML(XMLfileUnits& xmlconfig) override;
    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep) override;

    // TODO: cleanup?
    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("KartesianProfile");}

    static PluginBase* createInstance(){return new KartesianProfile();}

    double universalInvProfileUnit[3]; // Inv. Bin Sizes
    double universalProfileUnit[3]; // Bin Sizes
    long accumulatedDatasets; // Number of Datasets between output writes / profile resets
    double globalLength[3]; // Size of Domain
    double segmentVolume; // Size of one Sampling grid bin
    Domain* dom; // Reference to global Domain
    ProfileBase* _densProfile; //!< Reference to DensityProfile as it is needed by most other profiles

private:
    unsigned long _writeFrequency; // Write frequency for all profiles -> Length of recording frame before output
    unsigned long _initStatistics; // Timesteps to skip at start of the simulation
    unsigned long _profileRecordingTimesteps; // Record every Nth timestep during recording frame
    std::string _outputPrefix; // File prefix for all profiles
    std::string _mode;

    unsigned long _uIDs; //!< Total number of unique IDs with the selected Grid. This is the number of total bins in the Sampling grid.

    vector<ProfileBase*> _profiles; // vector holding all enabled profiles
    int _comms = 0; // total number of communications per bin needed by all profiles.

    // Needed for XML check for enabled profiles.
    bool _ALL = false;
    bool _DENSITY = false;
    bool _TEMPERATURE = false;
    bool _VELOCITY = false;

};


#endif //MARDYN_TRUNK_KARTESIAN2DPROFILE_H
