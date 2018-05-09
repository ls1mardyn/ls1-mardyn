//
// Created by Moritz Krügener
//

#include "COMaligner.h"

//! @brief will be called to read configuration
//!
//! All values have defaults and are not mandatory to be supplied<br>
//!
//! Defaults are:<br>
//!     bool X, Y, Z: all true (alignment on all Axis)<br>
//!     int interval: 25 (alignment every 25th simstep)<br>
//!     float correctionFactor: 1 (full alignment)<br>
//!
//! \param xmlconfig config.xml
void COMaligner::readXML(XMLfileUnits& xmlconfig){

    xmlconfig.getNodeValue("x", _alignX);
    xmlconfig.getNodeValue("y", _alignY);
    xmlconfig.getNodeValue("z", _alignZ);
    xmlconfig.getNodeValue("interval", _interval);
    xmlconfig.getNodeValue("correctionFactor", _alignmentCorrection);

    global_log -> info() << "[COMaligner] settings:" << std::endl;
    global_log -> info() << "                  x: " << _alignX << std::endl;
    global_log -> info() << "                  y: " << _alignY << std::endl;
    global_log -> info() << "                  z: " << _alignZ << std::endl;
    global_log -> info() << "                  interval: " << _interval << std::endl;
    global_log -> info() << "                  correctionFactor: " << _alignmentCorrection << std::endl;

}

//! @brief called before Forces are applied
//! calculates realignment motion that is applied after the forces have been applied
//!
//! \param particleContainer
//! \param domainDecomp
//! \param simstep
void COMaligner::beforeForces(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp,
                              unsigned long simstep) {

    global_log -> debug() << "[COMaligner] before forces called" << std::endl;

    if((simstep - 1) % _interval != 0){
        return;
    }

    // RESET
    for(unsigned d = 0; d < 3; d++){
        _balance[d] = 0.0;
        _motion[d] = 0.0;
    }
    _mass = 0;

    // ALL DIMENSIONS
    for(ParticleIterator tm = particleContainer->iterator(); tm.hasNext(); tm.next()){
        double partMass = tm->mass();
        _mass += partMass;
        for(unsigned d = 0; d < 3; d++){
            _balance[d] += tm->r(d) * partMass;
        }
    }

    // calculate Motion
    for(unsigned d = 0; d < 3; d++){
        _motion[d] = -_alignmentCorrection*((_balance[d] / _mass)-.5*_boxLength[d]);
    }
    global_log -> info() << "[COMaligner] motion is x: " << _motion[0] << " y: " << _motion[1] << " z: " << _motion[2] << std::endl;

    // TODO: Check for OpenMP implementation of above for-loop
    // TODO: Check if ComponentID check is important
    // TODO: CHECK IF MPI IMPLEMENTATION NECESSARY
}

//! @brief called after Forces are applied
//! applies the motion calculated earlier
//!
//! \param particleContainer
//! \param domainDecomp
//! \param domain
//! \param simstep
//! \param lmu
//! \param mcav
void COMaligner::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                         unsigned long simstep, std::list<ChemicalPotential> *lmu,
                         std::map<unsigned, CavityEnsemble> *mcav) {
    for(ParticleIterator tm = particleContainer->iterator(); tm.hasNext(); tm.next()){
        for(unsigned d = 0; d < 3; d++){
            tm->move(d, _motion[d]);
        }
    }
}