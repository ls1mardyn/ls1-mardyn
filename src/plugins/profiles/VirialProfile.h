//
// Created by Kruegener on 1/15/2019.
//

#ifndef MARDYN_VIRIAL_H
#define MARDYN_VIRIAL_H

#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"


class VirialProfile : public ProfileBase {
public:
	VirialProfile (DensityProfile* densProf) :
			_densityProfile{densProf}, _local3dProfile(), _global3dProfile() {};

	~VirialProfile () = default;

	void record (Molecule& mol, unsigned long uID) final {
		for(unsigned short d = 0; d < 3; d++){
			_local3dProfile[uID][d] += mol.Vi(d);
		}
	}

	void collectAppend (DomainDecompBase* domainDecomp, unsigned long uID) final {
		for(unsigned short d = 0; d < 3; d++){
			domainDecomp->collCommAppendDouble(_local3dProfile[uID][d]);
		}
	}

	void collectRetrieve (DomainDecompBase* domainDecomp, unsigned long uID) final {
		for(unsigned short d = 0; d < 3; d++){
			_global3dProfile[uID][d] = domainDecomp->collCommGetDouble();
		}
	}

	void output (string prefix, long unsigned accumulatedDatasets) final;

	void reset (unsigned long uID) final {
		for(unsigned d = 0; d < 3; d++){
			_local3dProfile[uID][d] = 0.0;
			_global3dProfile[uID][d] = 0.0;
		}
	}

	int comms () final { return 3; }

private:
	DensityProfile* _densityProfile;

	// Local 3D Profile
	std::map<unsigned, std::array<double, 3>> _local3dProfile;
	// Global 3D Profile
	std::map<unsigned, std::array<double, 3>> _global3dProfile;

	void writeDataEntry (unsigned long uID, ofstream& outfile) const final;
};


#endif //MARDYN_VIRIAL_H
