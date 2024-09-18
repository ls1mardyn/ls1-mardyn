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
#include "Filter.h"

class InternalCellProcessor;
class InternalProfiler{

    private:

    Filter filter;
    CellProcessor* cell_processor;
    int number_bins;
    double bin_width;
    /**
     * Is time average
     */
    std::vector<double> rdf_buffer;
    /**
     * Will be instantaneous
     */

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
     * Computes center of mass
     */
    std::array<double,3> GetCOM(Molecule* m);
    
    /**
     * Bins distances
     */
    void ProcessDistance(double distance, double pot);
    
    /**
     * Evaluates RDF and U(r) into the buffers, data is then ruined and clear is required
     */
    void GenerateInstantaneousData(ParticleContainer* particleContainer, Domain* domain);
    
    /**
     * Fills rdf buffer with 0s
     */
    void ResetBuffers();
    
    std::vector<double>& GetRDFValues();
    std::vector<double>& GetRNodes();
    double GetMeasuredSteps();


    private:
    void InitRNodes();
    void SetBinContainer(ParticleContainer* pc);
};

class InternalCellProcessor: public CellProcessor{

    private: 

    InternalProfiler* const my_profiler;

    public:
    InternalCellProcessor& operator=(const InternalCellProcessor&)=delete;
    InternalCellProcessor(const double cutoff, InternalProfiler* r);
    ~InternalCellProcessor(){};

    void initTraversal() override {}
    void preprocessCell(ParticleCell& cell) override {}
    void processCell(ParticleCell& cell) override;
    void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll=false) override;
    double processSingleMolecule(Molecule* m1, ParticleCell& cell) override {return 0.0;}
    void postprocessCell(ParticleCell& cell) override {}
    void endTraversal() override {}

    private:

    double DistanceBetweenCOMs(std::array<double,3>& com1, std::array<double,3>& com2);

    double PotentialCallBack(double eps, double sigma, double r);


};