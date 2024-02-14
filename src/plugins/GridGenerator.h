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

struct Element{
    int index;
};

class GridGenerator: public PluginBase{

    private:

    std::vector<Element> elements;

    public:
        void readXML(XMLfileUnits& xmlconfig) override;
        void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
        std::string getPluginName() override;

};