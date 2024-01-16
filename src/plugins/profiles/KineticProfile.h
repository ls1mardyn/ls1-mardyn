//
// Created by Kruegener on 10/22/2018.
//

#ifndef MARDYN_KINETICPROFILE_H
#define MARDYN_KINETICPROFILE_H

#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"

/**
 * @brief Records (NO OUTPUT) the 2xKinetic Profile of molecules per bin specified by Sampling grid in KartesianProfile.
 */
class KineticProfile final : public ProfileBase {
public:
    ~KineticProfile() final = default;
    void record(Molecule &mol, unsigned long uID) final  {
        double mv2 = 0.0;
        double Iw2 = 0.0;
        mol.calculate_mv2_Iw2(mv2, Iw2);
        _localProfile[uID] += mv2 + Iw2;
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        domainDecomp->collCommAppendDouble(_localProfile[uID]);
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        _globalProfile[uID] = domainDecomp->collCommGetDouble();
    }
    void output(std::string prefix, long unsigned accumulatedDatasets) final;
    void reset(unsigned long uID) final  {
        _localProfile[uID] = 0.0;
        _globalProfile[uID] = 0.0;
    }
    int comms() final {return 1;}

    double getGlobalKineticEnergy(unsigned long uid) const {
    	return _globalProfile.at(uid);
    }

private:
    // Local 1D Profile
    std::map<unsigned, double> _localProfile;
    // Global 1D Profile
    std::map<unsigned, double> _globalProfile;

    void writeDataEntry(unsigned long uID, std::ofstream &outfile) const final;

};


#endif //MARDYN_KINETICPROFILE_H
