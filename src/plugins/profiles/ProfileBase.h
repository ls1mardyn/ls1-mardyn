//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_PROFILEBASE_H
#define MARDYN_TRUNK_PROFILEBASE_H

#include "../../Domain.h"
#include "../../parallel/DomainDecompBase.h"
class KartesianProfile;

class ProfileBase {

public:
	virtual ~ProfileBase(){};
    void init(KartesianProfile* kartProf) {_kartProf = kartProf;};
    virtual void record(ParticleIterator *mol, unsigned long uID) = 0;
    virtual void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) = 0;
    virtual void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) = 0;
    virtual void output(string prefix) = 0;
    virtual void reset(unsigned long uID) = 0;
    // number of needed communications
    virtual int comms() = 0;

protected:
    // TODO: Add necessary maps for profiles
    // Local 1D Profile
    std::map<unsigned, long double> _localProfile;
    // Global 1D Profile
    std::map<unsigned, long double> _globalProfile;
    // Local 3D Profile
    std::map<unsigned, long double> _local3dProfile;
    // Global 3D Profile
    std::map<unsigned, long double> _global3dProfile;

    // output file prefix
    string _profilePrefix;
    // KartesianProfile managing class pointer for meta info
    KartesianProfile* _kartProf;
};


#endif //MARDYN_TRUNK_PROFILEBASE_H
