//
// Created by kruegener on 6/8/2018.
//

#ifndef MARDYN_TRUNK_WALLPOTENTIAL_H
#define MARDYN_TRUNK_WALLPOTENTIAL_H

#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

class WallPotential  : public PluginBase{

private:
    double _rhoW, _yc, _yOff, Delta;
    double* _eps_wi;
    double* _sig3_wi;
    double* _sig2_wi;
    double* _sig_wi;
    double* _uShift_9_3;
    double* _uPot_9_3;
    double* _uShift_10_4;
    double* _uPot_10_4;
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
                     Delta(0.),
                     _eps_wi(nullptr),
                     _sig3_wi(nullptr),
                     _sig2_wi(nullptr),
                     _sig_wi(nullptr),
                     _uShift_9_3(nullptr),
                     _uPot_9_3(nullptr),
                     _uShift_10_4(nullptr),
                     _uPot_10_4(nullptr),
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
        delete [] _uShift_10_4;
        delete [] _uPot_10_4;
    };

    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
        global_log -> debug() << "Wall enabled" << std::endl;
        _domain = domain;
    }

    void readXML (XMLfileUnits& xmlconfig);

    void beforeEventNewTimestep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ){}

    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {};

    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep){};

    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("WallPotential");}

    static PluginBase* createInstance(){return new WallPotential();}

    void initializeLJ93(const vector<Component> *components, double in_rhoWall, double in_sigWall, double in_epsWall,
                        vector<double> in_xi, vector<double> in_eta, double in_yOffWall, double in_yWallCut);

    void initializeLJ104(const vector<Component> *components, double in_rhoWall, double in_sigWall, double in_epsWall,
                         vector<double> in_xi, vector<double> in_eta, double in_yOffWall, double in_yWallCut);

    void calcTSLJ_9_3(ParticleContainer *partContainer);

    void calcTSLJ_10_4(ParticleContainer *partContainer);
};


#endif //MARDYN_TRUNK_WALLPOTENTIAL_H
