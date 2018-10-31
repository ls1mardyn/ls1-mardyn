//
// Created by Kruegener on 8/19/2018.
//

#include "KartesianProfile.h"
#include "plugins/profiles/ProfileBase.h"
#include "plugins/profiles/DensityProfile.h"
#include "plugins/profiles/Velocity3dProfile.h"
#include "plugins/profiles/VelocityAbsProfile.h"
#include "plugins/profiles/TemperatureProfile.h"
#include "plugins/profiles/KineticProfile.h"
#include "plugins/profiles/DOFProfile.h"

// DEBUGGING ONLY

bool CYLINDER_DEBUG = true;

/**
* @brief Read in Information about write/record frequencies, Sampling Grid and which profiles are enabled.
 * Also create needed profiles and initialize them. New Profiles need to be handled via the XML here as well.
* @param xmlconfig
*/
void KartesianProfile::readXML(XMLfileUnits &xmlconfig) {
    global_log -> debug() << "[KartesianProfile] enabled" << std::endl;
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    global_log->info() << "[KartesianProfile] Write frequency: " << _writeFrequency << endl;
    xmlconfig.getNodeValue("outputprefix", _outputPrefix);
    global_log->info() << "[KartesianProfile] Output prefix: " << _outputPrefix << endl;

    // TODO: options for cylinder / CLEAR WHAT UNIT DOES WHAT / DO WE WANT LINEAR RADIAL SAMPLING?
    xmlconfig.getNodeValue("x", samplInfo.universalProfileUnit[0]); // Doubles as Phi
    xmlconfig.getNodeValue("y", samplInfo.universalProfileUnit[1]); // Doubles as R
    xmlconfig.getNodeValue("z", samplInfo.universalProfileUnit[2]); // Doubles as H
    global_log->info() << "[KartesianProfile] Binning units: " << samplInfo.universalProfileUnit[0] << " " << samplInfo.universalProfileUnit[1] << " " << samplInfo.universalProfileUnit[2] << "\n";

    // CHECKING FOR ENABLED PROFILES
    int numProfiles = 0;
    xmlconfig.getNodeValue("profiles/density", _DENSITY);
    if(_DENSITY){
        global_log->info() << "[KartesianProfile] DENSITY PROFILE ENABLED\n";
        numProfiles++;
    }
    xmlconfig.getNodeValue("profiles/velocity", _VELOCITY);
    if(_VELOCITY){
        global_log->info() << "[KartesianProfile] VELOCITY PROFILE ENABLED\n";
        numProfiles++;
    }
    xmlconfig.getNodeValue("profiles/velocity3d", _VELOCITY3D);
    if(_VELOCITY){
        global_log->info() << "[KartesianProfile] VELOCITY3D PROFILE ENABLED\n";
        numProfiles++;
    }
    xmlconfig.getNodeValue("profiles/temperature", _TEMPERATURE);
    if(_TEMPERATURE){
        global_log->info() << "[KartesianProfile] TEMPERATURE PROFILE ENABLED\n";
        numProfiles++;
    }
    global_log->info() << "[KartesianProfile] Number of profiles: " << numProfiles << "\n";
    if(numProfiles < 1){
        global_log->warning() << "[KartesianProfile] NO PROFILES SPECIFIED -> Outputting all\n";
        _ALL = true;
    }

    // ADDING PROFILES
    // Need DensityProfile for Velocity*Profile / Temperature
    if(_DENSITY || _VELOCITY || _VELOCITY3D || _ALL){
        _densProfile = new DensityProfile();
        addProfile(_densProfile);
    }
    // Need Velocity for Temperature
    if(_VELOCITY || _ALL){
        _velAbsProfile = new VelocityAbsProfile(_densProfile);
        addProfile(_velAbsProfile);
    }
    if(_VELOCITY3D || _ALL){
        _vel3dProfile = new Velocity3dProfile(_densProfile);
        addProfile(_vel3dProfile);
    }
    if(_TEMPERATURE || _ALL){
        _dofProfile = new DOFProfile();
        addProfile(_dofProfile);
        _kineticProfile = new KineticProfile();
        addProfile(_kineticProfile);
        _tempProfile = new TemperatureProfile(_dofProfile, _kineticProfile);
        addProfile(_tempProfile);
    }

    xmlconfig.getNodeValue("timesteps/init", _initStatistics);
    global_log->info() << "[KartesianProfile] init statistics: " << _initStatistics << endl;
    xmlconfig.getNodeValue("timesteps/recording", _profileRecordingTimesteps);
    global_log->info() << "[KartesianProfile] profile recording timesteps: " << _profileRecordingTimesteps << endl;

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
void KartesianProfile::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
    // Get Global length
    for(unsigned d = 0; d < 3; d ++){
        samplInfo.globalLength[d] = domain->getGlobalLength(d);
        global_log->info() << "[KartesianProfile] globalLength " << samplInfo.globalLength[d] << "\n";
    }
    // Calculate sampling units
    // TODO: Adapt for option of cylinder uID, none of this actually makes too much sense for phi-bins > 1
    if(CYLINDER_DEBUG){
        double minXZ = this->samplInfo.globalLength[0];
        if(this->samplInfo.globalLength[2]<minXZ){
            minXZ = this->samplInfo.globalLength[2];
        }
        // R < .5minXZ -> R2max < .25minXZminXZ
        double Rmax = .5*minXZ;
        // TODO: WHY ARE THE DIMENSIONS SO JUMBLED???
        samplInfo.universalInvProfileUnit[0] = this->samplInfo.universalProfileUnit[0]/(2*M_PI);                   // delta_phi
        samplInfo.universalInvProfileUnit[1] = this->samplInfo.universalProfileUnit[1]/(Rmax);  // delta_R^2 -> delta R
        samplInfo.universalInvProfileUnit[2] = this->samplInfo.universalProfileUnit[2]/(samplInfo.globalLength[1]); // delta_H
    }
    else{
        for(unsigned i = 0; i < 3; i++){
            samplInfo.universalInvProfileUnit[i] = samplInfo.universalProfileUnit[i] / samplInfo.globalLength[i];
            global_log->info() << "[KartesianProfile] universalInvProfileUnit " << samplInfo.universalInvProfileUnit[i] << "\n";
        }
    }

    // Calculate total number of bins
    _uIDs = (unsigned long) (this->samplInfo.universalProfileUnit[0] * this->samplInfo.universalProfileUnit[1]
                                                                           * this->samplInfo.universalProfileUnit[2]);


    global_log->info() << "[KartesianProfile] number uID " << _uIDs << "\n";

    // Calculate bin Volume
    // TODO: again, cylinder, whats happening here?
    if(CYLINDER_DEBUG){
        samplInfo.segmentVolume = M_PI / (this->samplInfo.universalInvProfileUnit[1] * this->samplInfo.universalInvProfileUnit[2] * this->samplInfo.universalProfileUnit[0]);
    }
    else{
        samplInfo.segmentVolume = this->samplInfo.globalLength[0] * this->samplInfo.globalLength[1] * this->samplInfo.globalLength[2]
                                  / (this->samplInfo.universalProfileUnit[0]*this->samplInfo.universalProfileUnit[1]*this->samplInfo.universalProfileUnit[2]);
    }
    global_log->info() << "[KartesianProfile] segmentVolume " << samplInfo.segmentVolume << "\n";

    // Calculate Centre for cylinder coords
    samplInfo.universalCentre[0] = 0.5*samplInfo.globalLength[0];
    samplInfo.universalCentre[1] = 0;
    samplInfo.universalCentre[2] = 0.5*samplInfo.globalLength[2];

    global_log->info() << "[KartesianProfile] profile init" << std::endl;
    // Init profiles with sampling information and reset maps
    for(unsigned long uID = 0; uID < _uIDs; uID++) {
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
void KartesianProfile::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                               unsigned long simstep) {
    int mpi_rank = domainDecomp->getRank();

    unsigned xun, yun, zun;
    if ((simstep >= _initStatistics) && (simstep % _profileRecordingTimesteps == 0)) {
        unsigned long uID;

        // Loop over all particles and bin them with uIDs
        for(auto thismol = particleContainer->iterator(); thismol.isValid(); ++thismol){
            // Get uID
            // TODO: SELECT OPTION
            if(CYLINDER_DEBUG){
                uID = getCylUID(thismol);
            }
            else{
                uID = getUID(thismol);
            }
            // pass mol + uID to all profiles
            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->record(*thismol, uID);
            }
        }

        // Record number of Timesteps recorded since last output write
        _accumulatedDatasets++;
    }
    if ((simstep >= _initStatistics) && (simstep % _writeFrequency == 0)) {

        // COLLECTIVE COMMUNICATION
        global_log->info() << "[KartesianProfile] uIDs: " << _uIDs << " acc. Data: " << _accumulatedDatasets << "\n";

        // Initialize Communication with number of bins * number of total comms needed per bin by all profiles.
		domainDecomp->collCommInit(_comms * _uIDs);

		// Append Communications
        for(unsigned long uID = 0; uID < _uIDs; uID++){
            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->collectAppend(domainDecomp, uID);
            }
        }

        // Reduction Communication to get global values
        domainDecomp->collCommAllreduceSum();

        // Write global values in all bins in all profiles
        for(unsigned long uID = 0; uID < _uIDs; uID++) {
            for (unsigned i = 0; i < _profiles.size(); i++) {
                _profiles[i]->collectRetrieve(domainDecomp, uID);
            }
        }

        // Finalize Communication
        domainDecomp->collCommFinalize();

        // Initialize Output from rank 0 process
        if (mpi_rank == 0) {
            global_log->info() << "[KartesianProfile] Writing profile output" << std::endl;
            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->output(_outputPrefix + "_" + std::to_string(simstep), _accumulatedDatasets);
            }
        }

        // Reset profile arrays for next recording frame.
        for(unsigned long uID = 0; uID < _uIDs; uID++) {
            for (unsigned i = 0; i < _profiles.size(); i++) {
                _profiles[i]->reset(uID);
            }
        }
        _accumulatedDatasets = 0;
    }
}

unsigned long KartesianProfile::getUID(ParticleIterator& thismol) {
    auto xun = (unsigned) floor(thismol->r(0) * this->samplInfo.universalInvProfileUnit[0]);
    auto yun = (unsigned) floor(thismol->r(1) * this->samplInfo.universalInvProfileUnit[1]);
    auto zun = (unsigned) floor(thismol->r(2) * this->samplInfo.universalInvProfileUnit[2]);
    auto uID = (unsigned long) (xun * this->samplInfo.universalProfileUnit[1] * this->samplInfo.universalProfileUnit[2]
                           + yun * this->samplInfo.universalProfileUnit[2] + zun);
    return uID;
}

unsigned long KartesianProfile::getCylUID(ParticleIterator& thismol) {

    int xun,yun,zun;// (phiUn,r2Un,yun): bin number in a special direction, e.g. r2Un==5 corresponds to the 5th bin in the radial direction,
    long unID;	// as usual
    double xc,yc,zc; // distance of a particle with respect to the origin of the cylindrical coordinate system

    unID = -1; // initialization, causes an error message, if unID is not calculated in this method but used in record profile

    xc = thismol->r(0) - samplInfo.universalCentre[0];
    yc = thismol->r(1) - samplInfo.universalCentre[1];
    zc = thismol->r(2) - samplInfo.universalCentre[2];

    // transformation in polar coordinates
    double R2 = xc*xc + zc*zc;
    // TODO: CHANGED TO R
    double R = sqrt(R);
    double phi = asin(zc/sqrt(R2)) + ((xc>=0.0) ? 0:M_PI);
    if(phi<0.0) {phi = phi + 2.0*M_PI;}

    xun = (int)floor(phi * samplInfo.universalInvProfileUnit[0]);   // bin no. in phi-direction
    yun = (int)floor(R *  samplInfo.universalInvProfileUnit[1]);   // bin no. in R-direction
    zun = (int)floor(yc *  samplInfo.universalInvProfileUnit[2]);   // bin no. in H-direction

    if((xun >= 0) && (yun >= 0) && (zun >= 0) &&
       (xun < (int)samplInfo.universalProfileUnit[0]) && (yun < (int)samplInfo.universalProfileUnit[2]) && (zun < (int)samplInfo.universalProfileUnit[1]))
    {
        unID = (long) (xun * samplInfo.universalProfileUnit[1] * samplInfo.universalProfileUnit[2]
               + yun * samplInfo.universalProfileUnit[2] + zun);
    }
    else
    {
        global_log->error() << "INV PROFILE UNITS " << samplInfo.universalInvProfileUnit[0] << " " << samplInfo.universalInvProfileUnit[1] << " " << samplInfo.universalInvProfileUnit[2] << "\n";
        global_log->error() << "PROFILE UNITS " << samplInfo.universalProfileUnit[0] << " " << samplInfo.universalProfileUnit[1] << " " << samplInfo.universalProfileUnit[2] << "\n";
        global_log->error() << "Severe error!! Invalid profile unit (" << xun << " / " << yun << " / " << zun << ").\n\n";
        global_log->error() << "Coordinates off center (" << xc << " / " << yc << " / " << zc << ").\n";
        global_log->error() << "unID = " << unID << "\n";
        Simulation::exit(707);
    }
    return unID;
}

void KartesianProfile::addProfile(ProfileBase *profile) {
    global_log->info() << "PROFILE ADD\n";
    _profiles.push_back(profile);
    _comms += profile->comms();
}

