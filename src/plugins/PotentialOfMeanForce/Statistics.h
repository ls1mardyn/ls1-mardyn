/*
 * Created on Mon Jul 22 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"
#include "plugins/PotentialOfMeanForce/Region.h"

class StatisticsAdResS{

    private:

    struct ResRegion{
        std::array<double,3> low;
        std::array<double,3> high;
    };
    ResRegion fp, cg1, cg2, hy1, hy2;

    std::vector<double> fp_temp, cg1_temp, cg2_temp, hy1_temp, hy2_temp;
    std::string fp_component = "LJ";
    std::string file_name = "AdResSStatistics.output";
    int output_stride = 10;
    public:
    /**
     * Set the regions
     */
    void init(FPRegion& region);
    void Output2File(long step);
    void MeasureStatistics(ParticleContainer* pc){
        MeasureFPTemperature(pc);
        MeasureCGTemperature(pc);
        MeasureHYTemperature(pc);
    }

    private:
    void MeasureHYTemperature(ParticleContainer* pc);
    void MeasureFPTemperature(ParticleContainer* particleContainer);
    void MeasureCGTemperature(ParticleContainer* particleContainer);
};