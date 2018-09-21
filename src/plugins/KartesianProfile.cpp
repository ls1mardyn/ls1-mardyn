//
// Created by Kruegener on 8/19/2018.
//

#include "KartesianProfile.h"

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

    xmlconfig.getNodeValue("x", universalProfileUnit[0]);
    xmlconfig.getNodeValue("y", universalProfileUnit[1]);
    xmlconfig.getNodeValue("z", universalProfileUnit[2]);
    global_log->info() << "[KartesianProfile] Binning units: " << universalProfileUnit[0] << " " << universalProfileUnit[1] << " " << universalProfileUnit[2] << "\n";

    // CHECKING FOR ENABLED PROFILES
    int numProfiles = 0;
    xmlconfig.getNodeValue("profiles/density", _DENSITY);
    if(_DENSITY){
        global_log->info() << "[KartesianProfile] DENSITY PROFILE ENABLED\n";
        numProfiles++;
    }
    xmlconfig.getNodeValue("profiles/temperature", _TEMPERATURE);
    if(_TEMPERATURE){
        global_log->info() << "[KartesianProfile] TEMPERATURE PROFILE ENABLED\n";
        numProfiles++;
    }
    xmlconfig.getNodeValue("profiles/velocity", _VELOCITY);
    if(_VELOCITY){
        global_log->info() << "[KartesianProfile] VELOCITY PROFILE ENABLED\n";
        numProfiles++;
    }
    global_log->info() << "[KartesianProfile] Number of profiles: " << numProfiles << "\n";
    if(numProfiles < 1){
        global_log->warning() << "[KartesianProfile] NO PROFILES SPECIFIED -> Outputting all\n";
        _ALL = true;
    }
    // ADDING PROFILES
    // Need DensityProfile for Velocity3dProfile
    if(_DENSITY || _VELOCITY || _ALL){
        _densProfile = new DensityProfile();
        _densProfile->init(this);
        _profiles.push_back(_densProfile);
        _comms += _densProfile->comms();
    }
    if(_TEMPERATURE || _ALL){
        // TODO
    }
    // TODO: Different XML for different velocity profiles
    if(_VELOCITY || _ALL){
        ProfileBase* profile = new Velocity3dProfile();
        profile->init(this);
        _profiles.push_back(profile);
        _comms += profile->comms();
        ProfileBase* profile_abs = new VelocityAbsProfile();
        profile_abs->init(this);
        _profiles.push_back(profile_abs);
        _comms += profile_abs->comms();
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
    for(unsigned d = 0; d < 3; d ++){
        globalLength[d] = domain->getGlobalLength(d);
        global_log->info() << "[KartesianProfile] globalLength " << globalLength[d] << "\n";
    }
    for(unsigned i = 0; i < 3; i++){
        universalInvProfileUnit[i] = universalProfileUnit[i] / globalLength[i];
        global_log->info() << "[KartesianProfile] universalInvProfileUnit " << universalInvProfileUnit[i] << "\n";
    }
    _uIDs = (unsigned long) (this->universalProfileUnit[0] * this->universalProfileUnit[1]
                             * this->universalProfileUnit[2]);
    global_log->info() << "[KartesianProfile] number uID " << _uIDs << "\n";

    segmentVolume = this->globalLength[0] * this->globalLength[1] * this->globalLength[2]
                    / (this->universalProfileUnit[0]*this->universalProfileUnit[1]*this->universalProfileUnit[2]);
    global_log->info() << "[KartesianProfile] segmentVolume " << segmentVolume << "\n";
    dom = domain;

    global_log->info() << "[KartesianProfile] profile init" << std::endl;
    for(unsigned long uID = 0; uID < _uIDs; uID++) {
        for (unsigned i = 0; i < _profiles.size(); i++) {
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
        long int uID;
        // Loop over all particles and bin them with uIDs
        for(ParticleIterator thismol = particleContainer->iterator(); thismol.isValid(); ++thismol){
            // Calculate uID
            xun = (unsigned) floor(thismol->r(0) * this->universalInvProfileUnit[0]);
            yun = (unsigned) floor(thismol->r(1) * this->universalInvProfileUnit[1]);
            zun = (unsigned) floor(thismol->r(2) * this->universalInvProfileUnit[2]);
            uID = (unsigned long) (xun * this->universalProfileUnit[1] * this->universalProfileUnit[2]
                                   + yun * this->universalProfileUnit[2] + zun);
            // pass mol + uID to all profiles
            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->record(&thismol, uID);
            }
        }
        // Record number of Timesteps recorded since last output write
        accumulatedDatasets++;
    }
    if ((simstep >= _initStatistics) && (simstep % _writeFrequency == 0)) {
        // COLLECTIVE COMMUNICATION
        global_log->info() << "[KartesianProfile] uIDs: " << _uIDs << " acc. Data: " << accumulatedDatasets << "\n";
        // Initialize Communication with number of bins * number of total comms needed per bin by all profiles.
        domainDecomp->collCommInit(_comms*_uIDs);
        //global_log->info() << "[KartesianProfile] profile collectAppend" << std::endl;
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
                _profiles[i]->output(_outputPrefix + "_" + std::to_string(simstep));
            }
        }
        // Reset profile arrays for next recording frame.
        for(unsigned long uID = 0; uID < _uIDs; uID++) {
            for (unsigned i = 0; i < _profiles.size(); i++) {
                _profiles[i]->reset(uID);
            }
        }
        accumulatedDatasets = 0;
    }
}
