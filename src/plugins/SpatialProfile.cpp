//
// Created by Kruegener on 8/19/2018.
//

#include "SpatialProfile.h"
#include "plugins/profiles/ProfileBase.h"
#include "plugins/profiles/DensityProfile.h"
#include "plugins/profiles/Velocity3dProfile.h"
#include "plugins/profiles/VelocityAbsProfile.h"
#include "plugins/profiles/TemperatureProfile.h"
#include "plugins/profiles/KineticProfile.h"
#include "plugins/profiles/DOFProfile.h"
#include "plugins/profiles/VirialProfile.h"
#include "plugins/profiles/Virial2DProfile.h"

/**
* @brief Read in Information about write/record frequencies, Sampling Grid and which profiles are enabled.
 * Also create needed profiles and initialize them. New Profiles need to be handled via the XML here as well.
* @param xmlconfig
*/

void SpatialProfile::readXML(XMLfileUnits& xmlconfig) {
	global_log->debug() << "[SpatialProfile] enabled" << std::endl;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "[SpatialProfile] Write frequency: " << _writeFrequency << std::endl;
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "[SpatialProfile] Output prefix: " << _outputPrefix << std::endl;
	xmlconfig.getNodeValue("mode", _mode);
	global_log->info() << "[SpatialProfile] Mode " << _mode << std::endl;
	xmlconfig.getNodeValue("profiledComponent", _profiledCompString);
	global_log->info() << "[SpatialProfile] Profiled Component:" << _profiledCompString << std::endl;
	
	if (_profiledCompString != "all") {
		_profiledComp = std::stoi(_profiledCompString);
	}

	if (_mode == "cylinder") {
		xmlconfig.getNodeValue("r", samplInfo.universalProfileUnit[0]);
		xmlconfig.getNodeValue("h", samplInfo.universalProfileUnit[1]);
		xmlconfig.getNodeValue("phi", samplInfo.universalProfileUnit[2]);
		samplInfo.cylinder = true;
	} else if (_mode == "cartesian") {
		xmlconfig.getNodeValue("x", samplInfo.universalProfileUnit[0]);
		xmlconfig.getNodeValue("y", samplInfo.universalProfileUnit[1]);
		xmlconfig.getNodeValue("z", samplInfo.universalProfileUnit[2]);
		samplInfo.cylinder = false;
	} else {
		global_log->error() << "[SpatialProfile] Invalid mode. cylinder/cartesian" << std::endl;
		Simulation::exit(-1);
	}

	global_log->info() << "[SpatialProfile] Binning units: " << samplInfo.universalProfileUnit[0] << " "
					   << samplInfo.universalProfileUnit[1] << " " << samplInfo.universalProfileUnit[2] << "\n";

	// CHECKING FOR ENABLED PROFILES
	int numProfiles = 0;
	xmlconfig.getNodeValue("profiles/density", _DENSITY);
	if (_DENSITY) {
		global_log->info() << "[SpatialProfile] DENSITY PROFILE ENABLED\n";
		numProfiles++;
	}
	xmlconfig.getNodeValue("profiles/velocity", _VELOCITY);
	if (_VELOCITY) {
		global_log->info() << "[SpatialProfile] VELOCITY PROFILE ENABLED\n";
		numProfiles++;
	}
	xmlconfig.getNodeValue("profiles/velocity3d", _VELOCITY3D);
	if (_VELOCITY) {
		global_log->info() << "[SpatialProfile] VELOCITY3D PROFILE ENABLED\n";
		numProfiles++;
	}
	xmlconfig.getNodeValue("profiles/temperature", _TEMPERATURE);
	if (_TEMPERATURE) {
		global_log->info() << "[SpatialProfile] TEMPERATURE PROFILE ENABLED\n";
		numProfiles++;
	}
	xmlconfig.getNodeValue("profiles/virial", _VIRIAL);
	if (_VIRIAL) {
		global_log->info() << "[SpatialProfile] VIRIAL PROFILE ENABLED\n";
		numProfiles++;
	}
	xmlconfig.getNodeValue("profiles/virial2D", _VIRIAL2D);
	if (_VIRIAL2D) {
		global_log->info() << "[SpatialProfile] 2D VIRIAL PROFILE ENABLED\n";
		numProfiles++;
	}
	global_log->info() << "[SpatialProfile] Number of profiles: " << numProfiles << "\n";
	if (numProfiles < 1) {
		global_log->warning() << "[SpatialProfile] NO PROFILES SPECIFIED -> Outputting all\n";
		_ALL = true;
	}

	// ADDING PROFILES
	// Need DensityProfile for Velocity*Profile / Temperature / Virial
	if (_DENSITY || _VELOCITY || _VELOCITY3D || _VIRIAL || _VIRIAL2D || _ALL) {
		_densProfile = new DensityProfile();
		addProfile(_densProfile);
	}
	// Need Velocity for Temperature
	if (_VELOCITY || _ALL) {
		_velAbsProfile = new VelocityAbsProfile(_densProfile);
		addProfile(_velAbsProfile);
	}
	if (_VELOCITY3D || _ALL) {
		_vel3dProfile = new Velocity3dProfile(_densProfile);
		addProfile(_vel3dProfile);
	}
	if (_VIRIAL || _ALL) {
		_virialProfile = new VirialProfile(_densProfile);
		addProfile(_virialProfile);
	}
	if (_VIRIAL2D || _ALL) {
		_dofProfile = new DOFProfile();
		addProfile(_dofProfile);
		_kineticProfile = new KineticProfile();
		addProfile(_kineticProfile);
		_virial2DProfile = new Virial2DProfile(_densProfile, _dofProfile, _kineticProfile);
		addProfile(_virial2DProfile);
	}
	if (_TEMPERATURE || _ALL) {
		_dofProfile = new DOFProfile();
		addProfile(_dofProfile);
		_kineticProfile = new KineticProfile();
		addProfile(_kineticProfile);
		_tempProfile = new TemperatureProfile(_dofProfile, _kineticProfile);
		addProfile(_tempProfile);
	}

	xmlconfig.getNodeValue("timesteps/init", _initStatistics);
	global_log->info() << "[SpatialProfile] init statistics: " << _initStatistics << std::endl;
	xmlconfig.getNodeValue("timesteps/recording", _profileRecordingTimesteps);
	global_log->info() << "[SpatialProfile] profile recording timesteps: " << _profileRecordingTimesteps << std::endl;

}

/**
 *
 * @brief Initialize Arrays needed for calculating the profiles. Also get reference to domain for specific quantities.
 * All profiles will be reset here to 0 before starting the recording frame. The uIDs are can be calculated from the
 * globalLength and the number of divisions specified for the Sampling grid.
 *
 * @param particleContainer
 * @param domainDecomp
 * @param domain
 */
void SpatialProfile::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	// Get Global length
	for (unsigned d = 0; d < 3; d++) {
		samplInfo.globalLength[d] = domain->getGlobalLength(d);
		global_log->info() << "[SpatialProfile] globalLength " << samplInfo.globalLength[d] << "\n";
	}
	// Get global number of molecules
	samplInfo.globalNumMolecules = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);
	// Get number of molecules in FixRegion
	if (getNumFixRegion) {
		samplInfo.numMolFixRegion = (*getNumFixRegion)();
	} else {
		samplInfo.numMolFixRegion = 0;
	}
	global_log->info() << "getDomain was called with Molecules in Fix Region: " << samplInfo.numMolFixRegion << std::endl;

	// Calculate sampling units
	if (samplInfo.cylinder) {
		double minXZ = this->samplInfo.globalLength[0];
		if (this->samplInfo.globalLength[2] < minXZ) {
			minXZ = this->samplInfo.globalLength[2];
		}
		// R < .5minXZ -> R2max < .25minXZminXZ
		//double Rmax = .5*minXZ;
		double R2max = .24 * minXZ * minXZ;

		samplInfo.universalInvProfileUnit[0] =
				this->samplInfo.universalProfileUnit[0] / (R2max);                   // delta_R2
		samplInfo.universalInvProfileUnit[1] =
				this->samplInfo.universalProfileUnit[1] / (samplInfo.globalLength[1]);  // delta_H
		samplInfo.universalInvProfileUnit[2] = this->samplInfo.universalProfileUnit[2] / (2 * M_PI); // delta_Phi
		global_log->info() << "[CylinderProfile] dR: " << samplInfo.universalInvProfileUnit[0]
						   << " dH: " << samplInfo.universalInvProfileUnit[1]
						   << " dPhi: " << samplInfo.universalInvProfileUnit[2];
	} else {
		for (unsigned i = 0; i < 3; i++) {
			samplInfo.universalInvProfileUnit[i] = samplInfo.universalProfileUnit[i] / samplInfo.globalLength[i];
			global_log->info() << "[SpatialProfile] universalInvProfileUnit " << samplInfo.universalInvProfileUnit[i]
							   << "\n";
		}
	}

	// Calculate total number of bins
	_uIDs = (unsigned long) (this->samplInfo.universalProfileUnit[0] * this->samplInfo.universalProfileUnit[1]
							 * this->samplInfo.universalProfileUnit[2]);


	global_log->info() << "[SpatialProfile] number uID " << _uIDs << "\n";

	// Calculate bin Volume
	if (samplInfo.cylinder) {
		// When linearly sampling in R^2 domain -> Slice size stays constant with V = PI * (R_max^2 / nR) * (H / nH) * (1 / nPhi)
		samplInfo.segmentVolume = M_PI / (this->samplInfo.universalInvProfileUnit[0] *
										  this->samplInfo.universalInvProfileUnit[1] *
										  this->samplInfo.universalProfileUnit[2]);
	} else {
		samplInfo.segmentVolume =
				this->samplInfo.globalLength[0] * this->samplInfo.globalLength[1] * this->samplInfo.globalLength[2]
				/ (this->samplInfo.universalProfileUnit[0] * this->samplInfo.universalProfileUnit[1] *
				   this->samplInfo.universalProfileUnit[2]);
	}
	global_log->info() << "[SpatialProfile] segmentVolume " << samplInfo.segmentVolume << "\n";

	// Calculate Centre for cylinder coords
	samplInfo.universalCentre[0] = 0.5 * samplInfo.globalLength[0];
	samplInfo.universalCentre[1] = 0;
	samplInfo.universalCentre[2] = 0.5 * samplInfo.globalLength[2];

	global_log->info() << "[SpatialProfile] profile init" << std::endl;
	// Init profiles with sampling information and reset maps
	for (unsigned long uID = 0; uID < _uIDs; uID++) {
		for (unsigned i = 0; i < _profiles.size(); i++) {
			_profiles[i]->init(samplInfo);
			_profiles[i]->reset(uID);
		}
	}
}

/**
 * @brief Iterates over all molecules and passes them together with their Bin ID to the profiles for further processing.
 * If the current timestep hits the writefrequency the profile writes/resets are triggered here.
 * All of this only occurs after the initStatistics are passed.
 * @param particleContainer
 * @param domainDecomp
 * @param domain
 * @param simstep
 */
void SpatialProfile::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
							 unsigned long simstep) {
	int mpi_rank = domainDecomp->getRank();

	unsigned xun, yun, zun;
	if ((simstep >= _initStatistics) && (simstep % _profileRecordingTimesteps == 0)) {
		long uID;

		// Loop over all particles and bin them with uIDs
		for (auto thismol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); thismol.isValid(); ++thismol) {
			
			if ((_profiledCompString != "all") && (thismol->componentid() == _profiledComp-1)) {
				
				// Get uID
				if (samplInfo.cylinder) {
					uID = getCylUID(thismol);
					if (uID == -1) {
						// Invalid uID -> Molecule not in cylinder -> continue
						continue;
					}
				} else {
					uID = getCartesianUID(thismol);
				}
				// pass mol + uID to all profiles
				for (unsigned i = 0; i < _profiles.size(); i++) {
					_profiles[i]->record(*thismol, (unsigned) uID);
				}
			}
			
			if (_profiledCompString == "all"){
							
				// Get uID
				if (samplInfo.cylinder) {
					uID = getCylUID(thismol);
					if (uID == -1) {
						// Invalid uID -> Molecule not in cylinder -> continue
						continue;
					}
				} else {
					uID = getCartesianUID(thismol);
				}
				// pass mol + uID to all profiles
				for (unsigned i = 0; i < _profiles.size(); i++) {
					_profiles[i]->record(*thismol, (unsigned) uID);
				}
			}
		}

		// Record number of Timesteps recorded since last output write
		_accumulatedDatasets++;
	}
	if ((simstep >= _initStatistics) && (simstep % _writeFrequency == 0)) {

		// COLLECTIVE COMMUNICATION
		global_log->info() << "[SpatialProfile] uIDs: " << _uIDs << " acc. Data: " << _accumulatedDatasets << "\n";

		// Initialize Communication with number of bins * number of total comms needed per bin by all profiles.
		domainDecomp->collCommInit(_comms * _uIDs);

		// Append Communications
		for (unsigned long uID = 0; uID < _uIDs; uID++) {
			for (unsigned i = 0; i < _profiles.size(); i++) {
				_profiles[i]->collectAppend(domainDecomp, uID);
			}
		}

		// Reduction Communication to get global values
		domainDecomp->collCommAllreduceSum();

		// Write global values in all bins in all profiles
		for (unsigned long uID = 0; uID < _uIDs; uID++) {
			for (unsigned i = 0; i < _profiles.size(); i++) {
				_profiles[i]->collectRetrieve(domainDecomp, uID);
			}
		}

		// Finalize Communication
		domainDecomp->collCommFinalize();

		// Initialize Output from rank 0 process
		if (mpi_rank == 0) {
			global_log->info() << "[SpatialProfile] Writing profile output" << std::endl;
			for (unsigned i = 0; i < _profiles.size(); i++) {
				_profiles[i]->output(_outputPrefix + "_" + std::to_string(simstep), _accumulatedDatasets);
			}
		}

		// Reset profile arrays for next recording frame.
		for (unsigned long uID = 0; uID < _uIDs; uID++) {
			for (unsigned i = 0; i < _profiles.size(); i++) {
				_profiles[i]->reset(uID);
			}
		}
		_accumulatedDatasets = 0;
	}
}

/**
 * @brief getCartesianUID samples the domain cartesian coordinate bins.
 *
 * The calculation of uID has to be the same as used in the actual sampling profiles.
 * E.g. ProfileBase / Child classes of Profile Base.
 *
 * @param thismol
 * @return
 */
unsigned long SpatialProfile::getCartesianUID(ParticleIterator& thismol) {
	auto xun = (unsigned) floor(thismol->r(0) * samplInfo.universalInvProfileUnit[0]);
	auto yun = (unsigned) floor(thismol->r(1) * samplInfo.universalInvProfileUnit[1]);
	auto zun = (unsigned) floor(thismol->r(2) * samplInfo.universalInvProfileUnit[2]);
	auto uID = (unsigned long) (xun * samplInfo.universalProfileUnit[1] * samplInfo.universalProfileUnit[2]
								+ yun * samplInfo.universalProfileUnit[2] + zun);
	return uID;
}

/**
 *
 * @brief getCylUID samples the domain in cylinder coordinate bins.
 *
 * The sampling is done in linear increments both in the h- and phi-axis.
 * In the r-axis, the sampling is done linearly in R^2 space, which means the distance between sampling bins
 * R_i ^2 - R_i+1 ^2 is always the same. This results in smaller r distance between bins in outer regions.
 * This has the benefit, that the volume of all sampling bins is equal, meaning the volume of each ring slice, no matter
 * if sampled in the inner region of the cylinder or in the outer is the same.
 * This results from V_segment = PI*h*(r_i+1 ^2 - r_i ^2)
 * As said earlier, the R^2 difference terms are always the same if sampling is done linearly in R^2 - space.
 * The h term is a linear sampling size in h-axis, so also a constant as is PI.
 *
 * @param thismol
 * @return
 */
long SpatialProfile::getCylUID(ParticleIterator& thismol) {

	int phiUn, rUn, hUn;// (phiUn,rUn,yun): bin number in a special direction, e.g. rUn==5 corresponds to the 5th bin in the radial direction,
	long unID;    // as usual
	double xc, yc, zc; // distance of a particle with respect to the origin of the cylindrical coordinate system

	unID = -1; // initialization, causes an error message, if unID is not calculated in this method but used in record profile

	xc = thismol->r(0) - samplInfo.universalCentre[0];
	yc = thismol->r(1) - samplInfo.universalCentre[1];
	zc = thismol->r(2) - samplInfo.universalCentre[2];

	// transformation in polar coordinates
	double R2 = xc * xc + zc * zc;

	double phi = atan2(zc, xc); // asin(zc/sqrt(R2)) + ((xc>=0.0) ? 0:M_PI);
	if (phi < 0.0) { phi = phi + 2.0 * M_PI; }

	rUn = (int) floor(R2 * samplInfo.universalInvProfileUnit[0]);   // bin no. in R-direction
	hUn = (int) floor(yc * samplInfo.universalInvProfileUnit[1]);   // bin no. in H-direction
	phiUn = (int) floor(phi * samplInfo.universalInvProfileUnit[2]);   // bin no. in Phi-direction

	// Check if R is inside cylinder
	if (rUn >= (int) samplInfo.universalProfileUnit[0]) {
		return -1;
	}

	if ((rUn >= 0) && (hUn >= 0) && (phiUn >= 0) &&
		(rUn < (int) samplInfo.universalProfileUnit[0]) && (hUn < (int) samplInfo.universalProfileUnit[1]) &&
		(phiUn < (int) samplInfo.universalProfileUnit[2])) {
		unID = (long) (hUn * samplInfo.universalProfileUnit[0] * samplInfo.universalProfileUnit[2]
					   + rUn * samplInfo.universalProfileUnit[2] + phiUn);
	} else {
		global_log->error() << "INV PROFILE UNITS " << samplInfo.universalInvProfileUnit[0] << " "
							<< samplInfo.universalInvProfileUnit[1] << " " << samplInfo.universalInvProfileUnit[2]
							<< "\n";
		global_log->error() << "PROFILE UNITS " << samplInfo.universalProfileUnit[0] << " "
							<< samplInfo.universalProfileUnit[1] << " " << samplInfo.universalProfileUnit[2] << "\n";
		global_log->error() << "Severe error!! Invalid profile ID (" << rUn << " / " << hUn << " / " << phiUn
							<< ").\n\n";
		global_log->error() << "Severe error!! Invalid profile unit (" << R2 << " / " << yc << " / " << phi << ").\n\n";
		global_log->error() << "Coordinates off center (" << xc << " / " << yc << " / " << zc << ").\n";
		global_log->error() << "unID = " << unID << "\n";
		Simulation::exit(707);
	}
	return unID;
}

void SpatialProfile::addProfile(ProfileBase* profile) {
	global_log->info() << "[SpatialProfile] Profile added: \n";
	_profiles.push_back(profile);
	_comms += profile->comms();
}

