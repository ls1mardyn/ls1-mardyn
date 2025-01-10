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

struct ResRegion{
    std::array<double,3> low;
    std::array<double,3> high;
};

class RegionCellProcessor;

class RegionRDFProfiler{

    private:
    ResRegion& region;
    std::vector<double> centers;
    int total_bins=50;
    double bin_width;
    RegionCellProcessor* processor;
    int output_frequency=1;
    int sample_frequency=1;
    int measured_steps=0;//increases on every call to ProfileData
    std::string file_prefix="not-set";
    int N;//total molecules for the respective component
    int component_index=0;//component index in vector of ensemble
    std::string component_name="FP";

    public:

    RegionRDFProfiler(ResRegion& r);
    void init(int rc);
    void MeasureRDF(ParticleContainer* pc, unsigned long step);
    void SetFilePrefix(const char* name){
        file_prefix = name;
    }
    void SetTotalMolecules(int n){
        this->N=n;
    }
    void PrintOutput2Files(unsigned long simstep, ParticleContainer* pc);
    private:
    double RegionVolume(){
        std::array<double,3> low{0,0,0};
        std::array<double,3> high{0,0,0};
        low[1] = region.low[1]-_simulation.getcutoffRadius();
        low[2] = region.low[2]-_simulation.getcutoffRadius();
        high[1] = region.high[1]+_simulation.getcutoffRadius();
        high[2] = region.high[2]+_simulation.getcutoffRadius();

        double lx =region.high[0]-region.low[0];
        double ly =high[1]-low[1];
        double lz =high[2]-low[2];

        return lx*ly*lz;
    }
    void InitCenters(){
        for(int i=0;i<centers.size();++i){
            double center;
            center=(i+0.5)*bin_width;
            centers[i] = center;
        }
    }

    int CountMoleculesInRegion(ParticleContainer* pc);




};

class RegionCellProcessor:public CellProcessor{

    public:
    using Data = std::vector<double>;
    std::vector<Data> thread_data;
    std::vector<double> global_buffer;
    ResRegion& region;
    std::string component_name="FP";

    double bin_width;
    public:
    RegionCellProcessor& operator=(const RegionCellProcessor&)=delete;
    RegionCellProcessor(const double cutoff, int bins, double width, ResRegion& region);
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
    bool CellInRegion(ParticleCell& cell);
    void ProcessPairData(Molecule& m1, Molecule& m2);


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

    RegionRDFProfiler fp_profiler{fp};

    public:
    /**
     * Set the regions
     */
    void init(FPRegion& region);
    void Output2File(long step);
    void PrintRDFs2File(unsigned long step,ParticleContainer* pc){
        fp_profiler.PrintOutput2Files(step, pc);
    }
    void MeasureStatistics(ParticleContainer* pc){
        ClearAll();
        MeasureFPTemperature(pc);
        MeasureCGTemperature(pc);
        MeasureHYTemperature(pc);
    }

    void MeasureLocalRDFs(ParticleContainer* pc, unsigned long simstep){
        fp_profiler.MeasureRDF(pc,simstep);
    }
    void ClearAll();
    private:
    void MeasureHYTemperature(ParticleContainer* pc);
    void MeasureFPTemperature(ParticleContainer* particleContainer);
    void MeasureCGTemperature(ParticleContainer* particleContainer);

    
};