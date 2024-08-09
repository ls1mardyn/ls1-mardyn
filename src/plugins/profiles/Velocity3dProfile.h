//
// Created by Kruegener on 8/27/2018.
//

#ifndef MARDYN_TRUNK_VELOCITYPROFILE_H
#define MARDYN_TRUNK_VELOCITYPROFILE_H

#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"

#include <array>

class DensityProfile;

/**
 * @brief Outputs the XYZ velocity components per bin specified by Sampling grid in KartesianProfile.
 */
class Velocity3dProfile final : public ProfileBase {
public:
	Velocity3dProfile(DensityProfile * densProf) :
			_densityProfile(densProf), _local3dProfile(), _global3dProfile() {
	}
    ~Velocity3dProfile() final = default;
    void record(Molecule &mol, unsigned long uID) final  {
        for(unsigned short d = 0; d < 3; d++){
            _local3dProfile[uID][d] += mol.v(d);
        }
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        for(unsigned short d = 0; d < 3; d++){
            domainDecomp->collCommAppendDouble(_local3dProfile[uID][d]);
        }
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        for(unsigned short d = 0; d < 3; d++){
            _global3dProfile[uID][d] = domainDecomp->collCommGetDouble();
        }
    }
    void output(std::string prefix, long unsigned accumulatedDatasets) final;
    void reset(unsigned long uID) final  {
        for(unsigned d = 0; d < 3; d++){
            _local3dProfile[uID][d] = 0.0;
            _global3dProfile[uID][d] = 0.0;
        }
    }
    // set correct number of communications needed for this profile
    int comms() final {return 3;}

private:
    DensityProfile * _densityProfile;

    // Local 3D Profile
    std::map<unsigned, std::array<double,3>> _local3dProfile;
    // Global 3D Profile
    std::map<unsigned, std::array<double,3>> _global3dProfile;

    void writeDataEntry(unsigned long uID, std::ofstream &outfile) const final;
};


#endif //MARDYN_TRUNK_VELOCITYPROFILE_H
