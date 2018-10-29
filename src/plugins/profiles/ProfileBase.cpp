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
    for(unsigned z = 0; z < _kartProf->universalProfileUnit[2]; z++){
        outfile << (z+0.5) / _kartProf->universalInvProfileUnit[2] <<"  \t"; // Eintragen der z Koordinaten in Header
    }
    outfile << "\n";
    // Y - axis label
    for(unsigned y = 0; y < _kartProf->universalProfileUnit[1]; y++){
        double hval = (y + 0.5) / _kartProf->universalInvProfileUnit[1];
        outfile << hval << "  \t";
        // number density values
        for(unsigned z = 0; z < _kartProf->universalProfileUnit[2]; z++){
            for(unsigned x = 0; x < _kartProf->universalProfileUnit[0]; x++){
                auto unID = (unsigned long) (x * _kartProf->universalProfileUnit[0] * _kartProf->universalProfileUnit[2] + y * _kartProf->universalProfileUnit[1] + z);
                this->writeDataEntry(unID, outfile);
            }
        }
        outfile << "\n";
    }

}

void ProfileBase::writeSimpleMatrix(ofstream &outfile) {
    global_log->info() << "IMPLEMENT ME BEFORE USING\n";
}