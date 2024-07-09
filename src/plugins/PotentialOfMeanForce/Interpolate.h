#pragma once
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

/*
 * Created on Thu Jun 20 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

class Interpolate{
    
    private:

    std::vector<double> r_nodes;//always increases
    std::vector<double> g_nodes;

    public:

    std::vector<double>& GetRValues();
    std::vector<double>& GetGValues();
    double GetRDFAt(double r);
    double LinearInterpolation(int a, int b, double x);
    void ReadInRDF();
    /**
     * Essentially figure out between which values r is, and do central differences
     */
    double CentralFiniteDifference(double r);

    private:

    int GetUpperLimit(double r);
};