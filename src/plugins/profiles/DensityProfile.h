//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_DENSITYPROFILE_H
#define MARDYN_TRUNK_DENSITYPROFILE_H

#include "ProfileBase.h"

class DensityProfile : public ProfileBase {
public:
    DensityProfile(){global_log->info() << "[DensityProfile] enabled" << std::endl;};

    void record(ParticleIterator* mol, long int uID) override;
    void collectAppend() override;
    void collectRetrieve() override;
    void output() override;
    void reset() override;
};


#endif //MARDYN_TRUNK_DENSITYPROFILE_H
