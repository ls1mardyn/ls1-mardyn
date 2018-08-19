//
// Created by Kruegener on 8/19/2018.
//

#include "KartesianProfile.h"

bool DENSITY = true;

void KartesianProfile::readXML(XMLfileUnits &xmlconfig) {
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    global_log->info() << "[KartesianProfile] Write frequency: " << _writeFrequency << endl;
    xmlconfig.getNodeValue("outputprefix", _outputPrefix);
    global_log->info() << "[KartesianProfile] Output prefix: " << _outputPrefix << endl;

    xmlconfig.getNodeValue("x", _universalProfileUnit[0]);
    xmlconfig.getNodeValue("y", _universalProfileUnit[1]);
    xmlconfig.getNodeValue("z", _universalProfileUnit[2]);
    for(unsigned i = 0; i < 3; i++){
        _universalInvProfileUnit[i] = 1/_universalProfileUnit[i];
    }

    // TODO: add options to enable different profiles
    //xmlconfig.getNodeValue("options/option@[keyword='profileVirial']", doRecordVirialProfile);
    if(DENSITY){
        ProfileBase* profile = new DensityProfile();
        _profiles.push_back(profile);
    }
    //_density = new DensityProfile();

    xmlconfig.getNodeValue("timesteps/init", _initStatistics);
    global_log->info() << "[KartesianProfile] init statistics: " << _initStatistics << endl;
    xmlconfig.getNodeValue("timesteps/recording", _profileRecordingTimesteps);
    global_log->info() << "[KartesianProfile] profile recording timesteps: " << _profileRecordingTimesteps << endl;

}

void KartesianProfile::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                                 unsigned long simstep) {
    int mpi_rank = domainDecomp->getRank();
    long int uID;
    unsigned xun, yun, zun;
    if ((simstep >= _initStatistics) && (simstep % _profileRecordingTimesteps == 0)) {
        // TODO: RECORD PROFILES

        // Loop over all particles and bin them
        for(ParticleIterator thismol = particleContainer->iterator(); thismol.hasNext(); thismol.next()){
            // Calculate uID
            xun = (unsigned) floor(thismol->r(0) * this->_universalInvProfileUnit[0]);
            yun = (unsigned) floor(thismol->r(1) * this->_universalInvProfileUnit[1]);
            zun = (unsigned) floor(thismol->r(2) * this->_universalInvProfileUnit[2]);
            uID = xun * this->_universalProfileUnit[1] * this->_universalProfileUnit[2]
                   + yun * this->_universalProfileUnit[2] + zun;

            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->record(&thismol, uID);
            }
        }
    }
    if ((simstep >= _initStatistics) && (simstep % _writeFrequency == 0)) {
        for(unsigned i = 0; i < _profiles.size(); i++){
            _profiles[i]->collectAppend();
        }
        domainDecomp->collCommAllreduceSum();
        for(unsigned i = 0; i < _profiles.size(); i++){
            _profiles[i]->collectRetrieve();
        }
        if (mpi_rank == 0) {
            for(unsigned i = 0; i < _profiles.size(); i++){
                _profiles[i]->output();
            }
        }
        for(unsigned i = 0; i < _profiles.size(); i++){
            _profiles[i]->reset();
        }
    }
}
