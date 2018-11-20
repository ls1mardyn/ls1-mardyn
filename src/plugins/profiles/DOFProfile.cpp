//
// Created by Kruegener on 10/22/2018.
//

#include "DOFProfile.h"

void DOFProfile::output(string prefix, long unsigned accumulatedDatasets) {
    global_log->debug() << "[DOFProfile] has no output.\n";
}

void DOFProfile::writeDataEntry(unsigned long uID, ofstream &outfile) const {
    global_log->debug() << "[DOFProfile] has no writeDataEntry.\n";
}
