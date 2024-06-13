/*
 * Created on Tue Jun 11 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"
#include "Region.h"
#include "InteractionCellProcessor.h"
#include "Resolution.h"

//TODO: who owns the InteractionCellProcessor, AdResS or Simulation?

class InteractionSite;
class PMF:public PluginBase{
    using tracker = std::pair<InteractionSite,ResolutionType>;
    private:

    std::vector<double> nodes;
    std::vector<FPRegion> regions;//should create AT region with
    InteractionCellProcessor* adres_cell_processor;
    std::map<unsigned long, tracker> sites;
    ResolutionHandler resolution_handler;

    public:
    PMF();
    ~PMF(){}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override{}
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{};

    std::string getPluginName(){ return "PMF";}
    static PluginBase* createInstance() {return new PMF(); }

    public:

    std::vector<FPRegion>& GetRegions();
    ResolutionType GetMoleculeResolution(unsigned long idx);

};

class InteractionSite:public Site{

};