//
// Created by Kruegener on 10/22/2018.
//

#ifndef MARDYN_KINETICPROFILE_H
#define MARDYN_KINETICPROFILE_H

#include "ProfileBase.h"
#include "../KartesianProfile.h"

/**
 * @brief Records (NO OUTPUT) the 2xKinetic Profile of molecules per bin specified by Sampling grid in KartesianProfile.
 */
class KineticProfile : public ProfileBase {
public:
    ~KineticProfile() final {};
    void record(Molecule &mol, unsigned long uID) final  {
        double mv2 = 0.0;
        double Iw2 = 0.0;
        mol.calculate_mv2_Iw2(mv2, Iw2);
        _localProfile[uID] += mv2 + Iw2;
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


#endif //MARDYN_KINETICPROFILE_H
