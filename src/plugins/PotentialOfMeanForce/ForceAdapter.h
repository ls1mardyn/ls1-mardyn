/*
 * Created on Wed Jun 12 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#pragma once

#include "molecules/Molecule.h"
#include "molecules/Comp2Param.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

#include "Region.h"
#include "Resolution.h"
#include "WeightFunction.h"

class PMF;

class InteractionForceAdapter:public ParticlePairsHandler{
    private:

    struct DomainValues{

        double virial =0;
        double U_lj=0;
        double U_poles=0;
        double RF=0;

        void SetInDomain(Domain* domain){
            domain->setLocalUpot(U_lj/6.0 + U_poles + RF + domain->getLocalUpot());
            domain->setLocalVirial(virial + 3.0*RF+domain->getLocalVirial());
        }

        void ClearAll(){
            virial =0;
            U_lj=0;
            U_poles=0;
            RF=0;
        }

    };

    enum class InteractionType {onlyfp=0, mixed=1, onlycg=2, unspecified=3};
    public: 

    InteractionForceAdapter(ResolutionHandlerBase* res,PMF* pmf);
    virtual ~InteractionForceAdapter() override {}
    //Interface methods
    void init() override;
    void finish() override;
    /**
     * The distance vector wont be used for the vector in the hybrid case...or think about it
     */
    double processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ) override;
    ResolutionHandlerBase* GetResolutionHandler(){
        return resolution_handler;
    }

    private:

    //Has the case statement
    double processPairBackend(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool calcLJ, InteractionType interaction);
    
    //Checks interaction type and calls subroutines
    inline void PotForceType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, InteractionType interaction);
    inline void FluidPotType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, InteractionType interaction);


    //At some point there needs to be a mapping
    inline void PotForceOnlyCG(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ);
    inline void PotForceHybrid(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ);
    inline void PotForcePureHybridBackend(Molecule& HyM, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region);
    inline void PotForceHybridFPBackend(Molecule& HyM, Molecule& FPm2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region);
    inline void PotForceHybridCGBackend(Molecule& HyM, Molecule& FPm2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region);
    inline void FluidPotType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, bool hybrid);

    double PotentialOfMeanForce(double r);
    void ForceOfPotentialOfMeanForce(std::array<double,3>& f_com, double r);

    void AddAndMapForceToFP(std::array<double,3>& force, Molecule& molecule);
    void SubtractAndMapForceToFP(std::array<double,3>& force, Molecule& molecule);
    private:

    ResolutionHandlerBase* resolution_handler;
    PMF* adres;
    std::vector<ParticlePairs2PotForceAdapter::PP2PFAThreadData*> thread_data;
    DomainValues values;
};