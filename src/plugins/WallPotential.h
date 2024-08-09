//
// Created by kruegener on 6/8/2018.
// moved from /molecules/Wall.h by mheinen
//

#ifndef MARDYN_TRUNK_WALLPOTENTIAL_H
#define MARDYN_TRUNK_WALLPOTENTIAL_H

#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

/**
 * @brief WallPotential exerts the force of a Lennard-Jones potential on the Lennard-Jones centers of particles
 *
 * The user can choose between a 9-3 or 10-4 potential via the config.xml <br>
 * All significant values of these potentials can also be set there.<br>
 * Components to be considered can be set via <component> <br>
 * Their individual xi and eta can also be set in xml-format<br>
 * The Potential is applied once per simulation step in the Plugin::afterForces call<br>
 *
 * Usage:<br>
 * <br>
 * \code{.xml}
 * <plugin name="WallPotential">
    	<potential>93</potential>
        <density>1.0</density>
        <sigma>3.499976678</sigma>
        <epsilon>0.8</epsilon>
        <yoff>105.23725</yoff>
        <ycut>17.49988339</ycut>
        <width>34.99976678</width>
        <delta>1</delta>
        <component id="1">
          <xi>1.0</xi>
          <eta>1.0</eta>
        </component>
    </plugin>
 * \endcode
 */
class WallPotential  : public PluginBase{

private:
    double _rhoW, _yc, _yOff, _delta;
    double* _eps_wi;
    double* _sig3_wi;
    double* _sig2_wi;
    double* _sig_wi;
    double* _uShift_9_3;
    double* _uPot_9_3;
    double* _uShift_10_4_3;
    double* _uPot_10_4_3;
    unsigned _nc;
    double _dWidth;
    double _dWidthHalf;
    std::vector<bool> _bConsiderComponent;
    Domain* _domain;

    enum Potential {
        LJ9_3,
        LJ10_4,
    };
    int _potential;

public:
    WallPotential(): _rhoW(0.),
                     _yc(0.),
                     _yOff(0.),
                     _delta(0.),
                     _eps_wi(nullptr),
                     _sig3_wi(nullptr),
                     _sig2_wi(nullptr),
                     _sig_wi(nullptr),
                     _uShift_9_3(nullptr),
                     _uPot_9_3(nullptr),
                     _uShift_10_4_3(nullptr),
                     _uPot_10_4_3(nullptr),
                     _nc(0),
                     _dWidth(0.),
                     _dWidthHalf(0.){

    };
    ~WallPotential(){
        // free memory
        delete [] _eps_wi;
        delete [] _sig3_wi;
        delete [] _sig2_wi;
        delete [] _sig_wi;
        delete [] _uShift_9_3;
        delete [] _uPot_9_3;
        delete [] _uShift_10_4_3;
        delete [] _uPot_10_4_3;
    };

    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
        Log::global_log -> debug() << "[WallPotential] Wall enabled" << std::endl;
        _domain = domain;
    }

    void readXML (XMLfileUnits& xmlconfig) override;

    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep) override {};

    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("WallPotential");}

    static PluginBase* createInstance(){return new WallPotential();}

    void initializeLJ93(const std::vector<Component> *components, double in_rhoWall, double in_sigWall, double in_epsWall,
                        std::vector<double> in_xi, std::vector<double> in_eta, double in_yOffWall, double in_yWallCut);

    void initializeLJ1043(const std::vector<Component> *components, double in_rhoWall, double in_sigWall, double in_epsWall,
                          std::vector<double> in_xi, std::vector<double> in_eta, double in_yOffWall, double in_yWallCut);

    void calcTSLJ_9_3(ParticleContainer *partContainer);

    void calcTSLJ_10_4(ParticleContainer *partContainer);

    void siteWiseForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
                   unsigned long simstep) override;
};


#endif //MARDYN_TRUNK_WALLPOTENTIAL_H
