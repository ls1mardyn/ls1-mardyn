//
// Created by Kruegener on 10/22/2018.
//

#include "DOFProfile.h"

std::map<unsigned, long double> DOFProfile::getProfile() {
    return _globalProfile;
}

std::map<unsigned, long double>* DOFProfile::get3dProfile() {
    global_log->error() << "[DOFProfile] Does not support/use 3D Map. Use getProfile for 1D map instead." << std::endl;
    Simulation::exit(-1);
}

void DOFProfile::output(string prefix) {
    global_log->info() << "[DOFProfile] has no output.\n";
}
