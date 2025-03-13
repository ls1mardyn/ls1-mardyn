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

class RegionCellProcessor;

class RegionRDFProfiler{

    private:
    ResRegion& region;
    std::vector<double> centers;
    int total_bins=100;
    double bin_width;
    RegionCellProcessor* processor;
    int output_frequency=1000;
    int sample_frequency=1;
    int measured_steps=0;//increases on every call to ProfileData
    std::string file_prefix="not-set";

    bool measure_local_rdfs=true;

    public:

    RegionRDFProfiler(ResRegion& r);
    void init(int rc);
    void MeasureRDF(ParticleContainer* pc, unsigned long step);
    void SetFilePrefix(const char* name){
        file_prefix = name;
    }
    void PrintOutput2Files(unsigned long simstep, ParticleContainer* pc);
    private:
    //TODO: do on every step?
    double RegionVolume() {
        // TODO Note: we do not need to correct for halo volume, at least regular RDF does not do so
        double lx = region.high[0] - region.low[0];
        double ly = region.high[1] - region.low[1];
        double lz = region.high[2] - region.low[2];

        return lx*ly*lz;
    }

    void InitCenters() { for(int i = 0; i < centers.size(); ++i) centers[i] = (i+0.5)*bin_width; }

    /**
     * Extends region to include halo and gets N molecules, is there better option?
     */
    int CountMoleculesInRegion(ParticleContainer* pc);




};

//TODO: verification of components seems pretty primitive, better alternative?
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

    void SetComponentName(const char* name){
        component_name=name;
    }


    private:
    /*TODO: can this be better, the volume and the visited cells do not match.
            i.e., incomplete cells also belong in the region, but the region 
            dpes not match with the cells. What to do?
    */
    
    bool CellInRegion(ParticleCell& cell);
    void ProcessPairData(Molecule& m1, Molecule& m2);


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

    RegionRDFProfiler fp_profiler{fp};
    RegionRDFProfiler cg1_profiler{cg1};
    TheoreticalRDF fp_theoreticalRdf;
    TheoreticalRDF cg_theoreticalRdf;

    public:
    /**
     * Set the regions
     */
    void init(FPRegion& region);
    void Output2File(long step);
    void PrintRDFs2File(unsigned long step,ParticleContainer* pc){
        fp_profiler.PrintOutput2Files(step, pc);
        cg1_profiler.PrintOutput2Files(step,pc);
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
        fp_profiler.MeasureRDF(pc,simstep);
        cg1_profiler.MeasureRDF(pc,simstep);
        fp_theoreticalRdf.measure(pc);
        cg_theoreticalRdf.measure(pc);
    }
    void ClearAll();
    private:
    void MeasureHYTemperature(ParticleContainer* pc);
    void MeasureFPTemperature(ParticleContainer* particleContainer);
    void MeasureCGTemperature(ParticleContainer* particleContainer);

    
};