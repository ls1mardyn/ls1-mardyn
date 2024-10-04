//
// Created by Kruegener on 10/25/2018.
//

#include <sstream>

#include "ProfileBase.h"
#include "plugins/SpatialProfile.h"
#include "utils/mardyn_assert.h"

void ProfileBase::writeMatrix (std::ofstream& outfile) {
	if (_samplInfo.cylinder) {
		writeCylMatrix(outfile);
	} else {
		writeKartMatrix(outfile);
	}
}

void ProfileBase::writeKartMatrix (std::ofstream& outfile) {
	// Write Data
	// Assuming x = 1 -> projection along x axis
	// Z - axis label
	//global_log->info() << "[DensityProfile] output -> units: " << _kartProf->universalProfileUnit[0] << " " << _kartProf->universalProfileUnit[1] << " "<< _kartProf->universalProfileUnit[2] << "\n";
	for (unsigned z = 0; z < _samplInfo.universalProfileUnit[2]; z++) {
		outfile << (z + 0.5) / _samplInfo.universalInvProfileUnit[2] << "  \t"; // Eintragen der z Koordinaten in Header
	}
	outfile << "\n";
	// Y - axis label
	for (unsigned y = 0; y < _samplInfo.universalProfileUnit[1]; y++) {
		double hval = (y + 0.5) / _samplInfo.universalInvProfileUnit[1];
		outfile << hval << "  \t";
		// number density values
		for (unsigned z = 0; z < _samplInfo.universalProfileUnit[2]; z++) {
			for (unsigned x = 0; x < _samplInfo.universalProfileUnit[0]; x++) {
				// CRUCIAL:
				// Do not change unID calculation. Has to be the same as in SpatialProfile.cpp
				auto unID = (unsigned long) (
						x * _samplInfo.universalProfileUnit[1] * _samplInfo.universalProfileUnit[2] +
						y * _samplInfo.universalProfileUnit[2] + z);
				this->writeDataEntry(unID, outfile);
			}
		}
		outfile << "\n";
	}

}

void ProfileBase::writeSimpleMatrix (std::ofstream& outfile) {
	std::ostringstream error_message;
	error_message << "SIMPLE MATRIX OUTPUT NOT IMPLEMENTED!\n";
	MARDYN_EXIT(error_message);
}

void ProfileBase::writeCylMatrix (std::ofstream& outfile) {
	// Write Data
	// Assuming phi = 1 -> projection along phi axis
	// R2 - axis label
	//outfile << "R2 labels\n";
	for (unsigned r = 0; r < _samplInfo.universalProfileUnit[0]; r++) {
		outfile << 0.5*(sqrt(r+1)+sqrt(r)) / sqrt(_samplInfo.universalInvProfileUnit[0]) << " \t"; // Eintragen der R-Koordinate in Header korrigiert
	}
	//outfile << "END R2";
	outfile << "\n";
	// Y - axis label
	for (unsigned h = 0; h < _samplInfo.universalProfileUnit[1]; h++) {
		double hval = (h + 0.5) / _samplInfo.universalInvProfileUnit[1];
		outfile << hval << "  \t";
		for (unsigned phi = 0; phi < _samplInfo.universalProfileUnit[2]; phi++) {
			for (unsigned r = 0; r < _samplInfo.universalProfileUnit[0]; r++) {
				// CRUCIAL:
				// Do not change unID calculation. Has to be the same as in SpatialProfile.cpp
				auto unID = (long) (h * _samplInfo.universalProfileUnit[0] * _samplInfo.universalProfileUnit[2]
									+ r * _samplInfo.universalProfileUnit[2] + phi);
				this->writeDataEntry(unID, outfile);
			}
		}
		outfile << "\n";
	}
}
