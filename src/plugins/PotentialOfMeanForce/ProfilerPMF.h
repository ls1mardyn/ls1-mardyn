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

class InternalCellProcessor;
class InternalProfiler{

    private:

    CellProcessor* cell_processor;
    int number_bins;
    double bin_width;
    /**
     * Is time average
     */
    std::vector<int> bin_counts;
    /**
     * Will be instantaneous
     */
    std::vector<double> pot_vals;
    std::vector<double> pot_counts;
    int sample_frequency;
    int measured_steps;
    double measured_distance_squared;

    public:

    InternalProfiler();
    ~InternalProfiler(){

    }

    void init(ParticleContainer* pc, int bins, int freq);
    void SetBinContainer(ParticleContainer* pc);
    std::array<double,3> GetCOM(Molecule* m);
    void ProcessDistance(double distance, double pot);

    void GenerateInstantaneousData(ParticleContainer* particleContainer, Domain* domain);

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