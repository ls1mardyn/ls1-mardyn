//
// Created by Kruegener on 8/19/2018.
//

#include "KartesianProfile.h"

bool DENSITY = true;

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
    // TODO: add options to enable different profiles
    //xmlconfig.getNodeValue("options/option@[keyword='profileVirial']", doRecordVirialProfile);
    if(DENSITY){
        ProfileBase* profile = new DensityProfile();
        profile->init(this);
        _profiles.push_back(profile);
    }

    xmlconfig.getNodeValue("timesteps/init", _initStatistics);
    global_log->info() << "[KartesianProfile] init statistics: " << _initStatistics << endl;
    xmlconfig.getNodeValue("timesteps/recording", _profileRecordingTimesteps);
    global_log->info() << "[KartesianProfile] profile recording timesteps: " << _profileRecordingTimesteps << endl;


}

void KartesianProfile::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
    for(unsigned d = 0; d < 3; d ++){
        globalLength[d] = domain->getGlobalLength(d);
        global_log->info() << "[KartesianProfile] globalLength " << globalLength[d] << "\n";
    }
    for(unsigned i = 0; i < 3; i++){
        universalInvProfileUnit[i] = universalProfileUnit[i] / globalLength[i];
        global_log->info() << "[KartesianProfile] universalInvProfileUnit " << universalInvProfileUnit[i] << "\n";
    }
    _uIDs = this->universalProfileUnit[0] * this->universalProfileUnit[1]
            * this->universalProfileUnit[2];
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

void KartesianProfile::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                               unsigned long simstep) {
    int mpi_rank = domainDecomp->getRank();

    unsigned xun, yun, zun;
    if ((simstep >= _initStatistics) && (simstep % _profileRecordingTimesteps == 0)) {
        // TODO: RECORD PROFILES
        long int uID;
        // Loop over all particles and bin them
        for(ParticleIterator thismol = particleContainer->iterator(); thismol.hasNext(); thismol.next()){
            // Calculate uID
            xun = (unsigned) floor(thismol->r(0) * this->universalInvProfileUnit[0]);
            yun = (unsigned) floor(thismol->r(1) * this->universalInvProfileUnit[1]);
            zun = (unsigned) floor(thismol->r(2) * this->universalInvProfileUnit[2]);
            uID = xun * this->universalProfileUnit[1] * this->universalProfileUnit[2]
                  + yun * this->universalProfileUnit[2] + zun;

            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->record(&thismol, uID);
            }
        }
        accumulatedDatasets++;
    }
    if ((simstep >= _initStatistics) && (simstep % _writeFrequency == 0)) {
        // COLLECTIVE COMMUNICATION
        global_log->info() << "[KartesianProfile] profile collectAppend" << std::endl;
        for(unsigned long uID = 0; uID < _uIDs; uID++){
            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->collectAppend(domainDecomp, uID);
            }
        }
        global_log->info() << "[KartesianProfile] Allreduce" << std::endl;
        domainDecomp->collCommAllreduceSum();
        global_log->info() << "[KartesianProfile] profile collectRetrieve" << std::endl;
        for(unsigned long uID = 0; uID < _uIDs; uID++) {
            for (unsigned i = 0; i < _profiles.size(); i++) {
                _profiles[i]->collectRetrieve(domainDecomp, uID);
            }
        }
        if (mpi_rank == 0) {
            global_log->info() << "[KartesianProfile] profile output" << std::endl;
            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->output(_outputPrefix + "_" + std::to_string(simstep));
            }
        }
        global_log->info() << "[KartesianProfile] profile reset" << std::endl;
        for(unsigned long uID = 0; uID < _uIDs; uID++) {
            for (unsigned i = 0; i < _profiles.size(); i++) {
                _profiles[i]->reset(uID);
            }
        }
        accumulatedDatasets = 0;
    }
}
