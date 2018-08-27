//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_KARTESIAN2DPROFILE_H
#define MARDYN_TRUNK_KARTESIAN2DPROFILE_H

#include "PluginBase.h"
#include "Domain.h"
#include "profiles/DensityProfile.h"
#include "plugins/profiles/Velocity3dProfile.h"

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
class KartesianProfile : public PluginBase{

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

    double universalInvProfileUnit[3];
    double universalProfileUnit[3];
    long accumulatedDatasets;
    double globalLength[3];
    double segmentVolume;
    Domain* dom;
    ProfileBase* _densProfile; //!< Reference to DensityProfile as it is needed by most other profiles

private:
    unsigned long _writeFrequency;
    unsigned long _initStatistics;
    unsigned long _profileRecordingTimesteps;
    std::string _outputPrefix;
    std::string _mode;

    unsigned long _uIDs; //!< Total number of unique IDs with the selected Grid. This is the number of total cells in the grid.

    vector<ProfileBase*> _profiles;
    int _comms = 0;

    bool _ALL = false;
    bool _DENSITY = false;
    bool _TEMPERATURE = false;
    bool _VELOCITY = false;

};


#endif //MARDYN_TRUNK_KARTESIAN2DPROFILE_H
