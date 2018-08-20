//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_DENSITYPROFILE_H
#define MARDYN_TRUNK_DENSITYPROFILE_H

#include "ProfileBase.h"
#include "../KartesianProfile.h"

class DensityProfile : public ProfileBase {
public:
    DensityProfile()
    { global_log->info() << "[DensityProfile] enabled" << std::endl; };

    void record(ParticleIterator *mol, unsigned long uID) override;
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) override;
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) override;
    void output(string prefix) override;
    void reset(unsigned long uID) override;
};


#endif //MARDYN_TRUNK_DENSITYPROFILE_H
