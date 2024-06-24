/*
 * Created on Wed Jun 12 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

//based on the AdResSForceAdapter by Alex 

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
    enum class InteractionType {onlyfp=0, mixed=1, onlycg=2, unspecified=3};
    public: 

    InteractionForceAdapter(ResolutionHandler& res, PMF* pmf);
    virtual ~InteractionForceAdapter() override {}
    //Interface methods
    void init() override;
    void finish() override;
    /**
     * The distance vector wont be used for the vector in the hybrid case...or think about it
     */
    double processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ) override;



    private:

    //Has the case statement
    double processPairBackend(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool calcLJ, InteractionType interaction);
    
    //Checks interaction type and calls subroutines
    inline void PotForceType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, InteractionType interaction);

    //At some point there needs to be a mapping
    inline void PotForceOnlyCG(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ);
    inline void PotForceHybrid(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region);

    inline void HybridFluidPot(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, bool hybrid);

    /**
     * Uses formula: U_pmf = -k_b*T*log(g(r))
     * with k_b=1, T the ensemble temperature, g(r) the at-reference pair dist. function (rdf)
     */
    double PotentialOfMeanForce(double r);
    double SqrdDistanceBetweenCOMs(std::array<double,3> c1,std::array<double,3> c2);
    void ForceOfPotentialOfMeanForce(std::array<double,3>& f_com, double r);

    private:

    ResolutionHandler& resolution_handler;
    PMF* adres;
    std::vector<ParticlePairs2PotForceAdapter::PP2PFAThreadData*> thread_data;
};