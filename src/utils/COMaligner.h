//
// Created by Moritz Krügener
//

#ifndef MARDYN_TRUNK_COMALIGNER_H
#define MARDYN_TRUNK_COMALIGNER_H

#include "utils/PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"

class COMaligner : public PluginBase{

public:
    COMaligner(){};
    ~COMaligner(){};

    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
        global_log -> debug() << "COM Realignment enabled" << std::endl;

        for(unsigned d = 0; d < 3; d++){
            _boxLength[d] = domain->getGlobalLength(d);
        }

    }


    void readXML (XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep,
            std::list<ChemicalPotential> *lmu,
            std::map<unsigned, CavityEnsemble> *mcav
    );


    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("COMaligner");}

    static PluginBase* createInstance(){return new COMaligner();}

private:

    bool _alignX = true;
    bool _alignY = true;
    bool _alignZ = true;

    int _interval = 25;
    float _alignmentCorrection = 1.0f;

    double _motion[3];
    double _balance[3];
    double _mass;
    double _boxLength[3];
};


#endif //MARDYN_TRUNK_COMALIGNER_H
