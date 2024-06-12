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

class InteractionForceAdapter:public ParticlePairsHandler{
    public: 

    InteractionForceAdapter(ResolutionHandler& res);
    ~InteractionForceAdapter();
    void init() override;
    void finish() override;
    double processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ) override;


    private:

    ResolutionHandler& resolution_handler;
};