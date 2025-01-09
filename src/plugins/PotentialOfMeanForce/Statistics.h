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

struct ResRegion{
    std::array<double,3> low;
    std::array<double,3> high;
};

class StatisticsAdResS{

    private:
    ResRegion fp, cg1, cg2, hy1, hy2;

    double fp_temp, cg1_temp, cg2_temp, hy1_temp, hy2_temp;
    double cg_temp, hy_temp;
    int N_fp, N_cg1, N_cg2, N_hy1, N_hy2;
    std::string file_name = "AdResSStatistics.output";
    std::ofstream statistics;
    int output_stride = 1;

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

    void ClearAll();
};

class RegionRDFProfiler{

    private:
    ResRegion& region;
    std::vector<double> centers;
    std::vector<double> values;
    int total_bins;
    double bin_width;
    public:

    RegionRDFProfiler(ResRegion& r);
    void init(int tb);
    void MeasureRDF(ParticleContainer* pc);

    private:
    void InitCenters(){
        for(int i=0;i<centers.size();++i){
            double center;
            center=(i+0.5)*bin_width;
            centers[i] = center;
        }
    }

};

class RegionCellProcessor:public CellProcessor{

    public:
    using Data = std::vector<double>;
    std::vector<Data> thread_data;
    std::vector<double> global_buffer;
    FPRegion& region;

    double bin_width;

    public:
    RegionCellProcessor& operator=(const RegionCellProcessor&)=delete;
    RegionCellProcessor(const double cutoff, int bins, double width, FPRegion& region);
    ~RegionCellProcessor(){};

    void initTraversal() override {
        #pragma omp parallel
        {
            int own_id = mardyn_get_thread_num();
            std::fill(thread_data[own_id].begin(),thread_data[own_id].end(),0.0);
        }
    }

    void preprocessCell(ParticleCell& cell) override {}
    void processCell(ParticleCell& cell) override;
    void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll=false) override;
    double processSingleMolecule(Molecule* m1, ParticleCell& cell) override {return 0.0;}
    void postprocessCell(ParticleCell& cell) override {}
    void endTraversal() override;

    std::vector<double>& GetPairCountBuffer(){
        return global_buffer;
    }


    private:

    double DistanceBetweenCOMs(std::array<double,3>& com1, std::array<double,3>& com2);
    void ProcessPairData(Molecule& m1, Molecule& m2);
    void ProcessSingleData(Molecule& m2);


};