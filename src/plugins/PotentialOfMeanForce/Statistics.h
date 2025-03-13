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
#include "plugins/PotentialOfMeanForce/common.h"
#include "TheoreticalRDF.h"

struct ResRegion{
    std::array<double,3> low;
    std::array<double,3> high;
    std::string component_name = "not-set";
};

class StatisticsAdResS{

    private:

    ResRegion fp, cg1, cg2, hy1, hy2;
    double fp_temp, cg1_temp, cg2_temp, hy1_temp, hy2_temp;//per region
    double cg_temp, hy_temp;//average component
    int N_fp, N_cg1, N_cg2, N_hy1, N_hy2;
    std::string file_name = "AdResSStatistics.output";
    std::ofstream statistics;
    int output_stride = 1;

    TheoreticalRDF fp_theoreticalRdf;
    TheoreticalRDF cg_theoreticalRdf;

    public:
    /**
     * Set the regions
     */
    void init(FPRegion& region);
    void Output2File(long step);
    void PrintRDFs2File(unsigned long step,ParticleContainer* pc){
        fp_theoreticalRdf.writeFile("fp_theo_rdf_", step);
        cg_theoreticalRdf.writeFile("cg_theo_rdf_", step);
    }
    void MeasureStatistics(ParticleContainer* pc){
        ClearAll();
        MeasureFPTemperature(pc);
        MeasureCGTemperature(pc);
        MeasureHYTemperature(pc);
    }

    void MeasureLocalRDFs(ParticleContainer* pc, unsigned long simstep){
        fp_theoreticalRdf.measure(pc);
        cg_theoreticalRdf.measure(pc);
    }
    void ClearAll();
    private:
    void MeasureHYTemperature(ParticleContainer* pc);
    void MeasureFPTemperature(ParticleContainer* particleContainer);
    void MeasureCGTemperature(ParticleContainer* particleContainer);

    
};