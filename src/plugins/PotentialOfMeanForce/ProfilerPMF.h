/*
 * Created on Mon Aug 26 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"
#include "Region.h"
#include "InteractionCellProcessor.h"
#include "Resolution.h"
#include "WeightFunction.h"
#include "Interpolate.h"
#include "plugins/RDFAtCOM.h"
#include "ForceAdapter.h"
#include "common.h"

class InternalCellProcessor;
class InternalProfiler{

    private:

    struct ProfileBins{
        std::string name="not-set";
        int total_bins=0;
        double bin_width=0;
        std::vector<double> bin_centers;

        ProfileBins(const char* n):name{n}{

        }

        void InitBins(double dimension){
            bin_width = dimension/(double)total_bins;
            bin_centers.resize(total_bins,0.0);
            SetBinCenters();
            Log::global_log->info()<<"[PMF] Profiling "<<name<<std::endl;
            Log::global_log->info()<<"[PMF] Profer "<<name<<" has total bins"<< total_bins<<std::endl;

        }

        private:
        void SetBinCenters(){
            for(int i=0;i<total_bins;++i){
                double rmid;
                rmid = (i+0.5)*bin_width;
                bin_centers[i] = rmid;
            }
        }
    };
    bool measure_density=false;
    ProfileBins density{"density"};
    std::string rdf_file_name_="avg_rdf";
    InternalCellProcessor* cell_processor;
    int number_bins;
    double bin_width;
    std::vector<double> r_nodes;
    int sample_frequency;
    int measured_steps;//increases on every call to ProfileData
    double measured_distance_squared;
    int output_frequency;
    public:

    InternalProfiler();
    ~InternalProfiler(){

    }
    void ReadXML(XMLfileUnits& xmlconfig);
    void init(double rc);
    
    /**
     * Carries out the traversal and profiling
     */
    void ProfileData(ParticleContainer* pc, unsigned long simstep);
    void SetMeasureDensity(bool d){
        measure_density=d;
    }
    /**
     * Evaluates RDF and U(r) into the buffers, data is then ruined and clear is required
     */
    std::vector<double> GetRDF();    
    std::vector<double> GetDensity();    
    std::vector<double>& GetBinCounts();
    std::vector<double>& GetDensityCounts();
    std::vector<double>& GetBinCenters();
    std::vector<double>& GetDensityBinCenters(){
        return density.bin_centers;
    }
    double GetMeasuredSteps();
    void PrintOutput2Files(unsigned long simstep);

    private:

    void InitRNodes();
    void SetBinContainer(double rc);
};

class InternalCellProcessor: public CellProcessor{

    private: 
    using Data = std::vector<double>;
    std::vector<Data> thread_data;
    std::vector<double> global_buffer;
    std::vector<Data> thread_density_data;
    std::vector<double> global_density_buffer;
    double bin_width;
    double density_bin_width;//For indexing
    bool profile_density=true;
    public:
    InternalCellProcessor& operator=(const InternalCellProcessor&)=delete;
    InternalCellProcessor(const double cutoff, int bins, double width, double dwidth=0);
    ~InternalCellProcessor(){};

    void initTraversal() override {
        // std::fill(global_buffer.begin(),global_buffer.end(),0.0);
        // std::fill(global_density_buffer.begin(),global_density_buffer.end(),0.0);
        #pragma omp parallel
        {
            int own_id = mardyn_get_thread_num();
            std::fill(thread_data[own_id].begin(),thread_data[own_id].end(),0.0);
            std::fill(thread_density_data[own_id].begin(),thread_density_data[own_id].end(),0.0);
        }
    }
    void preprocessCell(ParticleCell& cell) override {}
    void processCell(ParticleCell& cell) override;
    void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll=false) override;
    double processSingleMolecule(Molecule* m1, ParticleCell& cell) override {return 0.0;}
    void postprocessCell(ParticleCell& cell) override {}
    void endTraversal() override;

    void SetMeasureDensity(bool b){
        profile_density = b;
    }

    void SetBinWidth(double w){
        density_bin_width = w;
    }

    std::vector<double>& GetPairCountBuffer(){
        return global_buffer;
    }

    std::vector<double>& GetDensityBuffer(){
        return global_density_buffer;
    }

    private:

    double DistanceBetweenCOMs(std::array<double,3>& com1, std::array<double,3>& com2);
    void ProcessPairData(Molecule& m1, Molecule& m2);
    void ProcessSingleData(Molecule& m2);

};