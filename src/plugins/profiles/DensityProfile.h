//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_DENSITYPROFILE_H
#define MARDYN_TRUNK_DENSITYPROFILE_H

#include "ProfileBase.h"
#include "../KartesianProfile.h"

class DensityProfile : public ProfileBase {
public:
	~DensityProfile(){};
    void record(ParticleIterator *mol, unsigned long uID) final  {
        _localProfile[uID] += 1;
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        domainDecomp->collCommAppendLongDouble(_localProfile[uID]);
        //global_log->info() << "[DensityProfile] localProfile " << uID << ": " << _localProfile[uID] << "\n";
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        _globalProfile[uID] = domainDecomp->collCommGetLongDouble();
        //global_log->info() << "[DensityProfile] globalProfile " << uID << ": " << _globalProfile[uID] << "\n";
    }
    void output(string prefix) override;
    void reset(unsigned long uID) final  {
        _localProfile[uID] = 0.0;
        _globalProfile[uID] = 0.0;
    }
    // set correct number of communications needed for this profile
    int comms() override {return 1;}
};


#endif //MARDYN_TRUNK_DENSITYPROFILE_H
