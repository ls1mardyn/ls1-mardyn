//
// Created by Heier on 07/29/2020.
//

#ifndef MARDYN_VIRIAL2D_H
#define MARDYN_VIRIAL2D_H

#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"


class DensityProfile;
class DOFProfile;
class KineticProfile;


class Virial2DProfile final : public ProfileBase {
public:
	Virial2DProfile(DensityProfile* densProf, DOFProfile * dofProf, KineticProfile * kinProf) :
			_densityProfile(densProf), _dofProfile(dofProf), _kineticProfile(kinProf), _local3dProfile(), _global3dProfile() {
			}

	~Virial2DProfile() final = default;


	void record(Molecule& mol, unsigned long uID) final {
		for (unsigned short d = 0; d < 3; d++) {
			_local3dProfile[uID][d] += mol.Vi(d);
		}
	}

	void collectAppend(DomainDecompBase* domainDecomp, unsigned long uID) final {
		for (unsigned short d = 0; d < 3; d++) {
			domainDecomp->collCommAppendDouble(_local3dProfile[uID][d]);
		}
	}

	void collectRetrieve(DomainDecompBase* domainDecomp, unsigned long uID) final {
		for (unsigned short d = 0; d < 3; d++) {
			_global3dProfile[uID][d] = domainDecomp->collCommGetDouble();
		}
	}


	void output(std::string prefix, long unsigned accumulatedDatasets) final;

	void reset(unsigned long uID) final {
		for (unsigned d = 0; d < 3; d++) {
			_local3dProfile[uID][d] = 0.0;
			_global3dProfile[uID][d] = 0.0;
		}
	}

	int comms() final { return 3; }

private:
	DensityProfile* _densityProfile;
	DOFProfile* _dofProfile;
	KineticProfile* _kineticProfile;

	// Local 3D Profile
	std::map<unsigned, std::array<double, 3>> _local3dProfile;
	// Global 3D Profile
	std::map<unsigned, std::array<double, 3>> _global3dProfile;

	// Only needed because its abstract, all output handled by output()
	void writeDataEntry(unsigned long uID, std::ofstream& outfile) const final;

};


#endif //MARDYN_VIRIAL_H
