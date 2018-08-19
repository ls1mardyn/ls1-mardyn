//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_PROFILEBASE_H
#define MARDYN_TRUNK_PROFILEBASE_H

#include "../../Domain.h"

class ProfileBase {

public:
    virtual void record(ParticleIterator* mol, long int uID) = 0;
    virtual void collectAppend() = 0;
    virtual void collectRetrieve() = 0;
    virtual void output() = 0;
    virtual void reset() = 0;
private:
    // TODO: Add necessary maps for profiles

};


#endif //MARDYN_TRUNK_PROFILEBASE_H
