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

    InternalCellProcessor* cell_processor;
    int number_bins;
    double bin_width;
    std::vector<double> rdf_buffer;
    std::vector<double> r_nodes;
    int sample_frequency;
    int measured_steps;//increases on every call to ProfileData
    double measured_distance_squared;

    public:

    InternalProfiler();
    ~InternalProfiler(){

    }

    void init(ParticleContainer* pc, int bins, int freq);
    
    /**
     * Carries out the traversal and profiling
     */
    void ProfileData(ParticleContainer* pc, unsigned long simstep);
    
    /**
     * Evaluates RDF and U(r) into the buffers, data is then ruined and clear is required
     */
    std::vector<double> GetInstantaneousData(Domain* domain);    
    std::vector<double>& GetBinCounts();
    std::vector<double>& GetBinCenters();
    double GetMeasuredSteps();


    private:
    /**
     * Fills rdf buffer with 0s
     */
    void ResetBuffer();
    void InitRNodes();
    void SetBinContainer(ParticleContainer* pc);
};

class InternalCellProcessor: public CellProcessor{

    private: 
    using Data = std::vector<double>;
    std::vector<Data> thread_data;
    std::vector<double> global_buffer;
    double bin_width;

    public:
    InternalCellProcessor& operator=(const InternalCellProcessor&)=delete;
    InternalCellProcessor(const double cutoff, int bins, double width);
    ~InternalCellProcessor(){};

    void initTraversal() override {
        std::fill(global_buffer.begin(),global_buffer.end(),0.0);
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

    std::vector<double>& GetBuffer(){
        return global_buffer;
    }

    private:

    double DistanceBetweenCOMs(std::array<double,3>& com1, std::array<double,3>& com2);
    void ProcessPairData(Molecule& m1, Molecule& m2);

};