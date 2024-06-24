//
// Created by Kruegener on 10/22/2018.
//

#ifndef MARDYN_DOFPROFILE_H
#define MARDYN_DOFPROFILE_H


#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"

/**
 * @brief Records (NO OUTPUT) the DOF of molecules per bin specified by Sampling grid in KartesianProfile.
 */
class DOFProfile final : public ProfileBase {
public:
    ~DOFProfile() final = default;
    void record(Molecule &mol, unsigned long uID) final  {
        _localProfile[uID] += 3.0 + (long double) (mol.component()->getRotationalDegreesOfFreedom());
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        domainDecomp->collCommAppendInt(_localProfile[uID]);
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        _globalProfile[uID] = domainDecomp->collCommGetInt();
    }
    void output(std::string prefix, long unsigned accumulatedDatasets) final;
    void reset(unsigned long uID) final  {
        _localProfile[uID] = 0;
        _globalProfile[uID] = 0;
    }
    int comms() final {return 1;}

    int getGlobalDOF(unsigned long uid) const {
    	return _globalProfile.at(uid);
    }

private:
    // Local 1D Profile
    std::map<unsigned, int> _localProfile;
    // Global 1D Profile
    std::map<unsigned, int> _globalProfile;

    void writeDataEntry(unsigned long uID, std::ofstream &outfile) const final;
};


#endif //MARDYN_DOFPROFILE_H
