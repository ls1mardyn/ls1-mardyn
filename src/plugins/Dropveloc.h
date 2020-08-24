/*
 * Dropveloc.h
 *
 *  Created on: 11 July 2020
 *      Author: Koch
 */

#ifndef MARDYN_TRUNK_DROPVELOC_H
#define MARDYN_TRUNK_DROPVELOC_H

//class DropvelocTest;
#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Calculates Center of mass of a spherical droplet and gradually approaches the drop diameter in the y direction<br>
* Then it takes the average speeds of all the particles in the drop and measures the speed of the moving droplet<br>
* The velocity of the surrounding gas phase is also determined, but can be neglected<br>
* Calculation happens every desired interval<br>
* Accuracy describes how exactly the approximate Center of mass reflects the real center of mass.<br>
* The unit of "accuracy" refers to the length units. A Value <1 should be given.<br>
* \code{.xml}
* <plugin>
*			<velocityradius>r-5</velocityradius>
*           <wallheight>4.75+?</wallheight>
*	        <component>1</component>
*           <accuracy>0.1</accuracy>
*   	    <measureinterval>100</measureinterval>
*           <maxiteration>100</maxiteration>
*			
* </plugin>
* \endcode
*/

class Dropveloc : public PluginBase{

private:
    //friend DropvelocTest;
	
	bool _enabled = true;

	int _measureinterval = 1;

    double _accuracy;
    double _wallheight;
    double _velocityradius;
    unsigned _component;
    double _alignmentCorrection = 1.0;
	double _motion[3];
    double _balance;
    double _mass = 0.0;
    double _boxLength[3];
	double _radius;
	double _xPos;
	double _yPos;
	double _zPos;
    double _maxiteration;
    std::map<unsigned long,double> _simsteparray;
    std::map<unsigned long,double> _centerarray;
    std::map<unsigned long, double> _velocityarray;
    std::ofstream _velocStream;

public:
 /*  Dropveloc(){
	 // SETUP
             _balance = 0.0;
            _motion[d] = 0.0;
    };
    ~Dropveloc(){};*/

    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
        global_log->debug() << "DropletVelocity enabled" << std::endl;
        for (unsigned d = 0; d < 3; d++) {
            _boxLength[d] = domain->getGlobalLength(d);
        }

           
        };

	void readXML (XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep) override;


    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("Dropveloc");}

    static PluginBase* createInstance(){return new Dropveloc();}

};


#endif //MARDYN_TRUNK_DROPVELOC_H
