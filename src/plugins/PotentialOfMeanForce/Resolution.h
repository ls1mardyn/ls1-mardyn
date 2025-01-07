#pragma once

#include "Region.h"
#include "molecules/Component.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PotentialOfMeanForce/common.h"

class ResolutionHandlerBase{
    public:
    
    virtual void init()=0;
    // virtual void CheckMoleculeResolution(Molecule& mol, ResolutionType target) =0;
    virtual void CheckContainerResolution(ParticleContainer* pc)=0;
    virtual ResolutionType GetMoleculeResolution(Molecule& mol)=0;

    std::vector<FPRegion>& GetRegions(){
        return regions;
    }

    protected:
    std::vector<FPRegion> regions;

};

class ResolutionHandler:public ResolutionHandlerBase{

    public:
    virtual void init()override;
    virtual void CheckContainerResolution(ParticleContainer* pc)override{
        
    }
    virtual ResolutionType GetMoleculeResolution(Molecule& mol) override;

    private:

    ResolutionType GetCOMResolution(std::array<double,3>& com, std::vector<FPRegion>& regions);

};

class ResolutionComponentHandler:public ResolutionHandlerBase{
    public:

    virtual void init() override;
    virtual ResolutionType GetMoleculeResolution(Molecule& m) override;
    virtual void CheckContainerResolution(ParticleContainer* pc) override;


    private:

    void CheckMoleculeResolution(Molecule& mol, ResolutionType target);
    void ModifyMoleculeResolution(Molecule& m, ResolutionType target);
    void AddCGComponent();
    void AddHYComponent();
    void Coarsen(Molecule& mol);
    void Refine(Molecule& mol);
    void Hybridize(Molecule& mol);
};