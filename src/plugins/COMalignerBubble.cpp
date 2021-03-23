/*
 * COMaligner.cpp
 *
 *  Created on: 7 May 2018
 *      Author: kruegener
 */

#include "COMalignerBubble.h"

//! @brief will be called to read configuration
//!
//! All values have defaults and are not mandatory to be supplied<br>
//!
//! Defaults are: <br>
//!     bool X, Y, Z: all true (alignment on all Axis) <br>
//!     int interval: 25 (alignment every 25th simstep) <br>
//!     float correctionFactor: 1 (full alignment) <br>
//!
//! \param xmlconfig  read from config.xml
void COMalignerBubble::readXML(XMLfileUnits& xmlconfig){

    xmlconfig.getNodeValue("x", _alignX);
    xmlconfig.getNodeValue("y", _alignY);
    xmlconfig.getNodeValue("z", _alignZ);
    xmlconfig.getNodeValue("interval", _interval);
    xmlconfig.getNodeValue("correctionFactor", _alignmentCorrection);

    // SANITY CHECK
    if(_interval < 1 || _alignmentCorrection < 0 || _alignmentCorrection > 1){
        global_log -> error() << "[COMalignerB] INVALID CONFIGURATION!!! DISABLED!" << std::endl;
        global_log -> error() << "[COMalignerB] HALTING SIMULATION" << std::endl;
        _enabled = false;
        // HALT SIM
        Simulation::exit(1);
        return;
    }

    global_log -> info() << "[COMalignerB] settings:" << std::endl;
    global_log -> info() << "                  x: " << _alignX << std::endl;
    global_log -> info() << "                  y: " << _alignY << std::endl;
    global_log -> info() << "                  z: " << _alignZ << std::endl;
    global_log -> info() << "                  interval: " << _interval << std::endl;
    global_log -> info() << "                  correctionFactor: " << _alignmentCorrection << std::endl;

    // Setting up different cases here to save on if statements in the simulation phase
    _dim_step = 1;
    if(_alignX){
        _dim_start = 0;
        if(_alignY){
            if(_alignZ){
                // X Y Z
                _dim_end = 3;
            }
            else{
                // X Y
                _dim_end = 2;
            }
        }
        else if(_alignZ){
            // X Z
            _dim_step = 2;
            _dim_end = 3;
        }
        else{
            // X
            _dim_end = 1;
        }
    }
    else if(_alignY){
        _dim_start = 1;
        if(_alignZ){
            // Y Z
            _dim_end = 3;
        }
        else{
            // Y
            _dim_end = 2;
        }
    }
    else if(_alignZ){
        // Z
        _dim_start = 2;
        _dim_end = 3;
    }
    else{
        _enabled = false;
    }

    global_log -> debug() << "[COMalignerB] dim settings are: " << _dim_start << " " << _dim_end << " " << _dim_step << std::endl;

}

//! @brief called before Forces are applied
//! calculates realignment motion that is applied after the forces have been applied<br>
//! only calculates motion on specified dimensions
//!
//! \param particleContainer
//! \param domainDecomp
//! \param simstep
void COMalignerBubble::beforeForces(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp,
                              unsigned long simstep) {

    if(_enabled) {

        global_log->debug() << "[COMalignerB] before forces called" << std::endl;

        if ((simstep - 1) % _interval != 0) {
            return;
        }

        // RESET
        for (unsigned d = 0; d < 3; d++) {
            _balance[d] = 0.0;
            _motion[d] = 0.0;
        }
        _mass = 0;

        // ITERATE OVER PARTICLES
        for (auto tm = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tm.isValid(); ++tm) {
            double partMass = tm->mass();
            _mass += partMass;
            for (int d = _dim_start; d < _dim_end; d += _dim_step) {
                _balance[d] += (tm->r(d)-_boxLength[d]*(std::floor(2*tm->r(d)/_boxLength[d]) ))* partMass;
            }
        }

        // COMMUNICATION
        domainDecomp->collCommInit(4);
        for (int d = 0; d < 3; d++) {
            domainDecomp->collCommAppendDouble(_balance[d]);
        }
        domainDecomp->collCommAppendDouble(_mass);
        domainDecomp->collCommAllreduceSum();
        for (int d = 0; d < 3; d++) {
            _balance[d] = domainDecomp->collCommGetDouble();
        }
        _mass = domainDecomp->collCommGetDouble();
        domainDecomp->collCommFinalize();

        // CALCULATE MOTION
        for (int d = _dim_start; d < _dim_end; d += _dim_step) {
            _motion[d] = -_alignmentCorrection * (_balance[d] / _mass);
        }
        global_log->info() << "[COMalignerB] motion is x: " << _motion[0] << " y: " << _motion[1] << " z: " << _motion[2]
                           << std::endl;

        // AVOID MOVES LARGER THAN ONE CUTOFF RADIUS
        double totalMotion = sqrt(_motion[0]*_motion[0]+_motion[1]*_motion[1]+_motion[2]*_motion[2]);
        if(totalMotion > _cutoff){
            double factor = _cutoff/totalMotion;
            for(int d = 0; d < 3; d++){
                _motion[d] *= factor;
            }
            global_log->info() << "[COMalignerB] Motion larger than Cutoff Radius. Reducing Motion" << _motion[2]
                               << std::endl;
            global_log->info() << "[COMalignerB] New motion is x: " << _motion[0] << " y: " << _motion[1] << " z: " << _motion[2]
                               << std::endl;
        }

        // MOVE
        for(auto tm = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tm.isValid(); ++tm){
            for (int d = _dim_start; d < _dim_end; d += _dim_step){
                tm->move(d, _motion[d]);
            }
        }

    }
    else{
        global_log->info() << "[COMalignerB] DISABLED, all dims set to 0" << std::endl;
    }

    // TODO: Check for OpenMP implementation of above for-loop
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
void COMalignerBubble::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                         unsigned long simstep) {

    // Moved to before Forces
    /*if(_enabled){
        for(auto tm = particleContainer->iterator(); tm.isValid(); ++tm){
            for (int d = _dim_start; d < _dim_end; d += _dim_step){
                tm->move(d, _motion[d]);
            }
        }
    }*/
}
