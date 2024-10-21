#pragma once

#include "Region.h"
#include "molecules/Component.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PotentialOfMeanForce/common.h"
#include <vector>
class InteractionSite;

//Based on Handler by Alex

class ResolutionHandler{
    private:

    Component cg;
    Component hy;

    public:

    void AddComponents2Ensemble();

    void AddCGComponent();

    void CheckResolution(ParticleContainer* pc, std::map<unsigned long, std::pair<InteractionSite,ResolutionType>>& res_map, std::vector<FPRegion>& regions);
    void CheckAndModifyMoleculeResolution(std::pair<InteractionSite,ResolutionType>& val, ResolutionType target_resolution);

};

class ResolutionComponentHandler{
    private:

    Component cg;
    Component hy;


    public:

    void init();

    void AddComponents2Ensemble();

    void AddCGComponent();
    ResolutionType GetMoleculeResolution(const Molecule& m);
    void CheckResolution(ParticleContainer* pc, std::vector<FPRegion>& regions);
    void CheckAndModifyMoleculeResolution(Molecule& mol, ResolutionType target_resolution);

    private:

    void Coarsen(Molecule& mol);
};
