//
// Created by Kruegener on 10/22/2018.
//

#include "KineticProfile.h"

std::map<unsigned, long double> KineticProfile::getProfile() {
    return _globalProfile;
}

std::map<unsigned, long double>* KineticProfile::get3dProfile() {
    global_log->error() << "[KineticProfile] Does not support/use 3D Map. Use getProfile for 1D map instead." << std::endl;
    Simulation::exit(-1);
}

void KineticProfile::output(string prefix) {
    global_log->info() << "[KineticProfile] has no output.\n";
}
