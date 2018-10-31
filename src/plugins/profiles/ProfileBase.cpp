//
// Created by Kruegener on 10/25/2018.
//

#include "ProfileBase.h"
#include "../KartesianProfile.h"

void ProfileBase::writeMatrix(ofstream &outfile) {
    // Write Data
    // Assuming x = 1 -> projection along x axis
    // Z - axis label
    //global_log->info() << "[DensityProfile] output -> units: " << _kartProf->universalProfileUnit[0] << " " << _kartProf->universalProfileUnit[1] << " "<< _kartProf->universalProfileUnit[2] << "\n";
    for(unsigned z = 0; z < _samplInfo.universalProfileUnit[2]; z++){
        outfile << (z+0.5) / _samplInfo.universalInvProfileUnit[2] <<"  \t"; // Eintragen der z Koordinaten in Header
    }
    outfile << "\n";
    // Y - axis label
    for(unsigned y = 0; y < _samplInfo.universalProfileUnit[1]; y++){
        double hval = (y + 0.5) / _samplInfo.universalInvProfileUnit[1];
        outfile << hval << "  \t";
        // number density values
        for(unsigned z = 0; z < _samplInfo.universalProfileUnit[2]; z++){
            for(unsigned x = 0; x < _samplInfo.universalProfileUnit[0]; x++){
                auto unID = (unsigned long) (x * _samplInfo.universalProfileUnit[0] * _samplInfo.universalProfileUnit[2] + y * _samplInfo.universalProfileUnit[1] + z);
                this->writeDataEntry(unID, outfile);
            }
        }
        outfile << "\n";
    }

}

void ProfileBase::writeSimpleMatrix(ofstream &outfile) {
    global_log->error() << "IMPLEMENT ME BEFORE USING\n";
}

void ProfileBase::writeCylMatrix(ofstream &outfile) {
    global_log->error() << "IMPLEMENT ME BEFORE USING\n";
    // TODO: varying order of loops in Domain.cpp:
    // TODO: line 1415ff / 1502ff: Phi -> H -> R^2
    // TODO: Line 1452ff: H -> R^2 -> Phi
    // Actually wrong in both cases. Only works if number of bins in Phi direction == 1 !!!
    // Otherwise, radial header gets rewritten all the time.
}
