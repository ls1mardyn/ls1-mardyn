//
// Created by Kruegener on 10/22/2018.
//

#include "DOFProfile.h"

void DOFProfile::output(std::string prefix, long unsigned accumulatedDatasets) {
    Log::global_log->debug() << "[DOFProfile] has no output.\n";
}

void DOFProfile::writeDataEntry(unsigned long uID, std::ofstream &outfile) const {
    Log::global_log->debug() << "[DOFProfile] has no writeDataEntry.\n";
}
