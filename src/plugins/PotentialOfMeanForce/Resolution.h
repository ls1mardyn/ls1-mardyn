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
    public:

    using tracker = std::pair<InteractionSite,ResolutionType>;
    using cg_map = std::map<unsigned long, tracker>;
    private:

    /**
     * Stores a COM site for every molecule (needs improvement)
     */
    std::map<unsigned long, tracker> sites;

    public:
    void SetMoleculeTrackerPosition(unsigned long idx, std::array<double, 3>& pos);
    cg_map& GetCGMap();
    std::array<double, 3>& GetMoleculeTrackerPosition(unsigned long idx);
    void CheckResolution(ParticleContainer* pc, std::vector<FPRegion>& regions);
    void CheckAndModifyMoleculeResolution(std::pair<InteractionSite,ResolutionType>& val, ResolutionType target_resolution);

    ResolutionType GetMoleculeResolution(unsigned long idx);
    InteractionSite GetMoleculeCOMSite(unsigned long idx);
    ResolutionType GetCOMResolution(std::array<double,3>& com, std::vector<FPRegion>& regions);

};

class ResolutionComponentHandler{
    private:

    Component cg;
    Component hy;

    public:

    void init();
    ResolutionType GetMoleculeResolution(const Molecule& m);
    void CheckResolution(ParticleContainer* pc, std::vector<FPRegion>& regions);
    void CheckAndModifyMoleculeResolution(Molecule& mol, ResolutionType target_resolution);

    private:
    /**
     * Called in init to create the cg site
     */
    void AddCGComponent();
    /**
     * Sets COM as r
     */
    void Coarsen(Molecule& mol);
};


/**
 * Stores velocity, position, force, potential
 */
class InteractionSite:public Site{
    private:
    double u_com;
    std::array<double,3> f_com;
    std::array<double,3> v_com;//not used
    std::array<double,3> r_com;

    public: 
    void SubForce(std::array<double, 3> f);
    void AddForce(std::array<double,3> f);
    void AddPotential(double pot);
    void SetPosition(std::array<double,3>& pos);
    void SetVelocity(std::array<double,3>& vel);

    std::array<double,3>& GetPosition();
    std::array<double,3>& GetVelocity();


};