#pragma once

#include "Region.h"
#include "molecules/Component.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"

#include <vector>
class InteractionSite;

//Based on Handler by Alex

class ResolutionHandler{

    public:

    void CheckResolution(ParticleContainer* pc, std::map<unsigned long, std::pair<InteractionSite,ResolutionType>>& res_map, std::vector<FPRegion>& regions);
    void CheckAndModifyMoleculeResolution(std::pair<InteractionSite,ResolutionType>& val, ResolutionType target_resolution);

};
