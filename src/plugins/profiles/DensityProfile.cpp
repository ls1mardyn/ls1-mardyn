//
// Created by Kruegener on 8/19/2018.
//

#include "DensityProfile.h"

inline void DensityProfile::record(ParticleIterator *mol, long int uID) {
    //global_log->info() << "[DensityProfile] record" << std::endl;
}

void DensityProfile::collectAppend() {
    global_log->info() << "[DensityProfile] collectAppend" << std::endl;
}

void DensityProfile::collectRetrieve() {
    global_log->info() << "[DensityProfile] collectRetrieve" << std::endl;
}

void DensityProfile::output() {
    global_log->info() << "[DensityProfile] output" << std::endl;
}

void DensityProfile::reset() {
    global_log->info() << "[DensityProfile] reset" << std::endl;
}
