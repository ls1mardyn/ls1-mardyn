//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_KARTESIAN2DPROFILE_H
#define MARDYN_TRUNK_KARTESIAN2DPROFILE_H

#include "PluginBase.h"
#include "Domain.h"
#include "profiles/DensityProfile.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

class KartesianProfile : public PluginBase{

public:
    KartesianProfile(){};
    ~KartesianProfile(){};

    /**
     *
     * Initialize Arrays
     *
     * @param particleContainer
     * @param domainDecomp
     * @param domain
     */
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
        global_log -> debug() << "[KartesianProfile] enabled" << std::endl;
    }

    /**
     * Read in Profile Steps
     * @param xmlconfig
     */
    void readXML(XMLfileUnits& xmlconfig) override;
    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep) override;

    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("KartesianProfile");}

    static PluginBase* createInstance(){return new KartesianProfile();}

private:
    unsigned long _writeFrequency;  // File Output / Reset frequency
    unsigned long _initStatistics;  // Timesteps to skip
    unsigned long _profileRecordingTimesteps;  // Sampling frequency
    std::string _outputPrefix;  // File name prefix

    long _universalInvProfileUnit[3];
    long _universalProfileUnit[3];

    vector<ProfileBase*> _profiles;

    DensityProfile* _density;
};


#endif //MARDYN_TRUNK_KARTESIAN2DPROFILE_H
