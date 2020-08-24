/*
 * Dropveloc.h
 *
 *  Created on: 11 July 2020
 *      Author: Koch
 */

#include "Dropveloc.h"

//! @brief will be called to read configuration
//!
//!
//!
//! \param xmlconfig  read from config.xml
void Dropveloc::readXML(XMLfileUnits& xmlconfig){

    xmlconfig.getNodeValue("velocityradius", _velocityradius);
    xmlconfig.getNodeValue("wallheight", _wallheight);
    xmlconfig.getNodeValue("component", _component);
    xmlconfig.getNodeValue("accuracy", _accuracy);
    xmlconfig.getNodeValue("measureinterval", _measureinterval);
    xmlconfig.getNodeValue("maxiteration", _maxiteration);


	// SANITY CHECK
    if(_measureinterval < 1 ){
        global_log -> error() << "[Dropveloc] INVALID CONFIGURATION!!! DISABLED!" << std::endl;
        global_log -> error() << "[Dropveloc] HALTING SIMULATION" << std::endl;
        _enabled = false;
        // HALT SIM
        Simulation::exit(1);
        return;
    }

	global_log -> info() << "[Dropveloc] settings:" << std::endl;
    global_log->info() << "                  velocityradius: " << _velocityradius << std::endl;
    global_log->info() << "                  wallheight: " << _wallheight << std::endl;
    global_log->info() << "                  component: " << _component << std::endl;
    global_log->info() << "                  accuracy: " << _accuracy << std::endl;
    global_log -> info() << "                measureinterval: " << _measureinterval << std::endl;
    global_log->info() << "                maxiteration: " << _maxiteration << std::endl;

}



void Dropveloc::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    if (_enabled) {

        global_log->debug() << "[Dropveloc] before forces called" << std::endl;

        if ((simstep - 1) % _measureinterval != 0) {
            return;
        }

        double y1 = _boxLength[1];
        double y2 = _wallheight;
        double molposition;
        double DROPLETCENTER;
        double Dropletcenter = 0;
        double a1;
        double a2;
        double approxstep = _velocityradius / 10;
        bool insidedroplet = false;
        double j = 0;
        double l = 0;
        int w =0;

        do
        {
            // RESET

            _balance = 0.0;
            _mass = 0.0;

            // ITERATE OVER PARTICLES
            for (auto tm = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tm.isValid(); ++tm) {

                unsigned COMPONENT = tm->getComponentLookUpID();
                molposition = tm->r(1);

                if (COMPONENT == (_component-1) && molposition<y1 && molposition>y2) {
                    double partMass = tm->mass();
                    _mass += partMass;
                    _balance += molposition*partMass;

                }

            }


            // COMMUNICATION
            domainDecomp->collCommInit(2);
            domainDecomp->collCommAppendDouble(_balance);
            domainDecomp->collCommAppendDouble(_mass);
            domainDecomp->collCommAllreduceSum();
            _balance = domainDecomp->collCommGetDouble();
            _mass = domainDecomp->collCommGetDouble();
            domainDecomp->collCommFinalize();

            // CALCULATE POSITION OF DROPLETCENTER
            DROPLETCENTER = _balance / _mass;

            a1 = y1 - DROPLETCENTER;
            a2 = DROPLETCENTER - y2;
            

            if (a1 > _velocityradius || a2 > _velocityradius) {
                if (a1 > a2) {
                    y1 = y1 - approxstep;
                    insidedroplet = false;
                    w += 1;
                }
                else if (a1 < a2) {
                    y2 = y2 + approxstep;
                    insidedroplet = false;
                    w += 1;
                }
                else {
                    y1 = y1 - approxstep;
                    y2 = y2 + approxstep;
                    insidedroplet = false;
                    w += 1;
                }
            }
            else
            {
                y1 = DROPLETCENTER + _velocityradius;           // damit ist a1=velocityradius
                y2 = DROPLETCENTER - _velocityradius;           // damit ist a2=velocityradius
                                                                //      --> immer in diese else Schleife
                if (Dropletcenter == 0) {

                    Dropletcenter = DROPLETCENTER;
                    insidedroplet = false;
                    approxstep = approxstep/10;
                    global_log->info() << "Wiederholung: "  << w <<"\t Dropletcenter (first)= "<< DROPLETCENTER << endl;
                    w += 1;
                }

                else if ((DROPLETCENTER<(Dropletcenter + _accuracy)) && (DROPLETCENTER>(Dropletcenter - _accuracy))) {
                    Dropletcenter= DROPLETCENTER;
                    insidedroplet = true;
                    global_log->info() << "Wiederholung: " << w << "\t Dropletcenter (end)= " << DROPLETCENTER << endl;
                    w += 1;
                }

                else {
                     Dropletcenter = DROPLETCENTER;
                    insidedroplet = false; 
                    approxstep = approxstep/2;
                    global_log->info() << "Wiederholung: " << w << "\t Dropletcenter ("<<w<<")= " << DROPLETCENTER << endl;
                    w += 1;
                }
               
            }
             l += 1;
        } while ((insidedroplet == false)&&(l<_maxiteration));

        double dropvelocity=0;

        if (simstep/_measureinterval >1) {

            global_log->info() << "Measureintervall of Velocity = " << simstep << endl;

         //CALCULATE VELOCITY OUT OF THE DROPLETCENTER-POSITIONS

            dropvelocity = ((DROPLETCENTER - _centerarray[(simstep - _measureinterval) / _measureinterval]) / (simstep - _simsteparray[(simstep-_measureinterval) / _measureinterval]))/0.0005;

         //SAFE AND WRITE VELOITY (AND DROPLETCENTER-Position) IN A NEW FILE
            global_log->info() << "vdrop = "<< dropvelocity;
            global_log->info() << "\t y-Dropposition = " <<DROPLETCENTER<<endl;


        }

        //SAFE ACTUAL DROPLETCENTER-POSITION, ACTUAL SIMSTEP AND VELOCITYm

        _simsteparray[simstep/_measureinterval] = simstep;
        _centerarray[simstep / _measureinterval] = DROPLETCENTER;
        _velocityarray[simstep / _measureinterval] = dropvelocity;

        
	}
}		

void Dropveloc::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                         unsigned long simstep) {

    _velocStream.open("Dropvelocity.txt", std::ios::out);
    _velocStream << "Dropletvelocity and movement of y-Dropletcenter: " << endl;
    _velocStream << endl;
    _velocStream << "Timestep \t y-Dropcenterposition \t Velocity in y-Direction" << endl;
    if (simstep / _measureinterval > 1) {

        for (int forstep = 1; forstep < (simstep / _measureinterval); ++forstep) {
            _velocStream << _simsteparray[forstep] << " \t\t\t" << _centerarray[forstep] << " \t\t\t " << _velocityarray[forstep] << endl;
        }
    }
    _velocStream.close();
    
}

        
