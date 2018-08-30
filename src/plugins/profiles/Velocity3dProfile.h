//
// Created by Kruegener on 8/27/2018.
//

#ifndef MARDYN_TRUNK_VELOCITYPROFILE_H
#define MARDYN_TRUNK_VELOCITYPROFILE_H

#include "ProfileBase.h"
#include "../KartesianProfile.h"

/**
 * @brief Outputs the XYZ velocity components per bin specified by Sampling grid in KartesianProfile.
 */
class Velocity3dProfile : public ProfileBase {
public:
    ~Velocity3dProfile() final = default;
    void record(ParticleIterator *mol, unsigned long uID) final  {
        for(unsigned short d = 0; d < 3; d++){
            _local3dProfile[d][uID] += (*mol)->v(d);
        }
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        for(unsigned short d = 0; d < 3; d++){
            domainDecomp->collCommAppendLongDouble(_local3dProfile[d][uID]);
        }
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        for(unsigned short d = 0; d < 3; d++){
            _global3dProfile[d][uID] = domainDecomp->collCommGetLongDouble();
        }
    }
    void output(string prefix) final;
    void reset(unsigned long uID) final  {
        for(unsigned d = 0; d < 3; d++){
            _local3dProfile[d][uID] = 0.0;
            _global3dProfile[d][uID] = 0.0;
        }
    }
    // set correct number of communications needed for this profile
    int comms() final {return 3;}
};


#endif //MARDYN_TRUNK_VELOCITYPROFILE_H
