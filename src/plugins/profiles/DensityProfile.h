//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_DENSITYPROFILE_H
#define MARDYN_TRUNK_DENSITYPROFILE_H

#include "ProfileBase.h"
#include "../KartesianProfile.h"

/**
 * @brief Outputs the number density of molecules per bin specified by Sampling grid in KartesianProfile.
 */
class DensityProfile : public ProfileBase {
public:
	~DensityProfile() final {};
    void record(Molecule &mol, unsigned long uID) final  {
        _localProfile[uID] += 1;
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        domainDecomp->collCommAppendLongDouble(_localProfile[uID]);
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        _globalProfile[uID] = domainDecomp->collCommGetLongDouble();
    }
    void output(string prefix) final;
    void reset(unsigned long uID) final  {
        _localProfile[uID] = 0.0;
        _globalProfile[uID] = 0.0;
    }
    int comms() final {return 1;}
    std::map<unsigned, long double> getProfile();
    std::map<unsigned, long double>* get3dProfile();

private:
    // Local 1D Profile
    std::map<unsigned, long double> _localProfile;
    // Global 1D Profile
    std::map<unsigned, long double> _globalProfile;
};


#endif //MARDYN_TRUNK_DENSITYPROFILE_H
