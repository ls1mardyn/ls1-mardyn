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
#include "IBI_Math.h"

class IBIPairsHandler : public ParticlePairsHandler {
public:
    IBIPairsHandler();
    virtual ~IBIPairsHandler() override;
    //Interface methods
    void init() override;
    void finish() override;
    /**
     * The distance vector wont be used for the vector in the hybrid case...or think about it
     */
    double processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ) override;

    FunctionPL& getForceFunction() { return force_function; }
    FunctionPL& getPotentialFunction() { return potential_function; }

private:
    inline void PotForceOnlyCG(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ);
    std::array<double, 3> getActingForce(const std::array<double, 3>& dr, double r);

    /// is either PMF or updated version
    FunctionPL force_function;
    /// is either PMF or updated version
    FunctionPL potential_function;
    /// MT potential data
    std::vector<ParticlePairs2PotForceAdapter::PP2PFAThreadData*> thread_data;
};