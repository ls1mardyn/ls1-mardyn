//
// Created by Kruegener on 10/22/2018.
//

#include "KineticProfile.h"

void KineticProfile::output(string prefix) {
    global_log->info() << "[KineticProfile] has no output.\n";
}

void KineticProfile::writeDataEntry(unsigned long uID, ofstream &outfile) const {
    global_log->debug() << "[KineticProfile] has no writeDataEntry.\n";
}
