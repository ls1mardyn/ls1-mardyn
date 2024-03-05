//
// Created by Kruegener on 10/22/2018.
//

#include "KineticProfile.h"

void KineticProfile::output(std::string prefix, long unsigned accumulatedDatasets) {
    Log::global_log->info() << "[KineticProfile] has no output.\n";
}

void KineticProfile::writeDataEntry(unsigned long uID, std::ofstream &outfile) const {
    Log::global_log->debug() << "[KineticProfile] has no writeDataEntry.\n";
}
