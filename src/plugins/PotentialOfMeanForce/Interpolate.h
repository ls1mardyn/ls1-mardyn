/*
 * Created on Thu Jun 20 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"


class Interpolate{
    
    private:

    std::vector<double> x_values;//function value
    std::vector<double> y_values;
    double default_value;

    public:
    Interpolate()=default;
    Interpolate(double def);     
    void SetXValues(std::vector<double>& v);
    void SetYValues(std::vector<double>& v);

    std::vector<double>& GetXValues();
    std::vector<double>& GetYValues();
    double InterpolateAt(double r);
    double LinearInterpolation(int a, int b, double x);
    /**
     * Essentially figure out between which values r is, and do central differences
     */
    double CentralFiniteDifference(double r);

    private:

    int GetUpperLimit(double r);
};