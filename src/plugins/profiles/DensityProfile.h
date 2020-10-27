//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_DENSITYPROFILE_H
#define MARDYN_TRUNK_DENSITYPROFILE_H

#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"

/**
 * @brief Outputs the number density of molecules per bin specified by Sampling grid in KartesianProfile.
 */
class DensityProfile final : public ProfileBase {
public:
	~DensityProfile() final = default;
    void record(Molecule &mol, unsigned long uID) final  {
        _localProfile[uID] += 1;
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        domainDecomp->collCommAppendInt(_localProfile[uID]);
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        _globalProfile[uID] = domainDecomp->collCommGetInt();
    }
    void output(string prefix, long unsigned accumulatedDatasets) final;
    void reset(unsigned long uID) final  {
        _localProfile[uID] = 0;
        _globalProfile[uID] = 0;
    }
    int comms() final {return 1;}

    int getGlobalNumber (unsigned long uid) const {
    	return _globalProfile.at(uid);
    }

private:
    // Local 1D Profile
    std::map<unsigned, int> _localProfile;
    // Global 1D Profile
    std::map<unsigned, int> _globalProfile;

    void writeDataEntry(unsigned long uID, ofstream &outfile) const final;
};


#endif //MARDYN_TRUNK_DENSITYPROFILE_H
