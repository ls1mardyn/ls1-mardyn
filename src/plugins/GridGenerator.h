/*
 * Created on Wed Feb 14 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#include <functional>
#include <optional>
#include <vector>

#include <plugins/profiles/ProfileBase.h>
#include "PluginBase.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

struct ElementInfo{
    int index;
    //std::vector<double> 
    double x_width;
    double y_width;
    double z_width;
    double volume;

};

class GridGenerator: public PluginBase{

    private:
    ElementInfo element_information;
    std::vector<int> elements_per_dimension{1,1,1};
    std::vector<double> element_width_per_dimension{0,0,0};
    std::vector<double> lower_corner;
    std::vector<double> upper_corner;
    int total_elements;


    public:
        void readXML(XMLfileUnits& xmlconfig) override;
        void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
        void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
        void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
        void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
        void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override{}
        void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{}
        std::string getPluginName() override;



        static PluginBase* createInstance() { return new GridGenerator(); }

    private:

    void MeshEntireDomain();
    void SetElementInfo();
    void SetTotalElements();

};