/*
 * COMaligner.h
 *
 *  Created on: 7 May 2018
 *      Author: kruegener
 */

#ifndef MARDYN_TRUNK_COMALIGNER_H
#define MARDYN_TRUNK_COMALIGNER_H

class COMalignerTest;
#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin that moves all particles so the center of mass aligns with the center of the simulation box.
*
* Individual dimensions X,Y,Z can be toggled on/off for the alignment.<br>
* Alignment happens once every interval-simsteps.<br>
* A correction factor can be set from 0-1 to control the alignment with 1 being full alignment and 0 no alignment at all.<br>
* <b>HALO must not be present</b> for the alignment. Halo would lead to incorrect alignment.<br>
* This is guaranteed by calling the alignment in the beforeForces step of the simulation.
* \endcode
*/
class COMaligner : public PluginBase{

private:
    friend COMalignerTest;

    // DEFAULT: ALIGN IN ALL DIMENSIONS
    bool _alignX = true;
    bool _alignY = true;
    bool _alignZ = true;

    bool _enabled = true;

    int _dim_start = 0;
    int _dim_end = 3;
    int _dim_step = 1;

    // DEFAULT: EVERY FRAME FULL ALIGNMENT
    int _interval = 1;
    float _alignmentCorrection = 1.0f;

    double _motion[3];
    double _balance[3];
    double _mass = 0.0;
    double _boxLength[3];
    double _cutoff;

public:
    COMaligner(){
        // SETUP
        for (unsigned d = 0; d < 3; d++) {
            _balance[d] = 0.0;
            _motion[d] = 0.0;
        }
    };
    ~COMaligner(){};

    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
        Log::global_log -> debug() << "COM Realignment enabled" << std::endl;

        for(unsigned d = 0; d < 3; d++){
            _boxLength[d] = domain->getGlobalLength(d);
        }

        _cutoff = .9*particleContainer->getCutoff();

    }

	/** @brief Read in XML configuration for COMAligner
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
        <plugin name="COMaligner">
                <x>BOOL</x>  <!-- align along x-axis (default): true | do no align along x-axis: false -->
                <y>BOOL</y>  <!-- align along y-axis (default): true | do no align along y-axis: false -->
                <z>BOOL</z>  <!-- align along z-axis (default): true | do no align along z-axis: false -->
                <interval>INT</interval>  <!-- frequency of algignment in number of time steps between alignments -->
                <correctionFactor>FLOAT</correctionFactor>  <!-- correction factor [0-1] to apply with 1 meaning full alignment and 0 no alignment at all -->
        </plugin>
	   \endcode
	 */
    void readXML (XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep) override;


    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("COMaligner");}

    static PluginBase* createInstance(){return new COMaligner();}

};


#endif //MARDYN_TRUNK_COMALIGNER_H
