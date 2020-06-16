//
// Created by Kruegener on 10/22/2018.
//

#ifndef MARDYN_TEMPERATUREPROFILE_H
#define MARDYN_TEMPERATUREPROFILE_H

#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"

class DOFProfile;
class KineticProfile;

/**
 * @brief Outputs the temperature of molecules per bin specified by Sampling grid in KartesianProfile.
 */
class TemperatureProfile final : public ProfileBase {
public:
	TemperatureProfile(DOFProfile * dofProf, KineticProfile * kinProf) :
			_dofProfile(dofProf), _kineticProfile(kinProf), _localProfile(), _globalProfile() {
	}
    ~TemperatureProfile() final = default;
    void record(Molecule &mol, unsigned long uID) final  {
        _localProfile[uID] += 1;
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        domainDecomp->collCommAppendLongDouble(_localProfile[uID]);
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        _globalProfile[uID] = domainDecomp->collCommGetLongDouble();
    }
    void output(string prefix, long unsigned accumulatedDatasets) final;
    void reset(unsigned long uID) final  {
        _localProfile[uID] = 0.0;
        _globalProfile[uID] = 0.0;
    }
    int comms() final {return 1;}

private:
    DOFProfile * _dofProfile;
    KineticProfile * _kineticProfile;

    // Local 1D Profile
    std::map<unsigned, long double> _localProfile;
    // Global 1D Profile
    std::map<unsigned, long double> _globalProfile;

    void writeDataEntry(unsigned long uID, ofstream &outfile) const final;
};


#endif //MARDYN_TEMPERATUREPROFILE_H
