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
#include "plugins/RDFAtCOM.h"
#include "IBIPairsHandler.h"

class RDFProfiler {
public:
    RDFProfiler();
    ~RDFProfiler() = default;
    void init(int bins);
    /// Carries out the traversal and profiling
    void ProfileData(ParticleContainer* pc);
    /// Fills rdf buffer with 0s and set measured steps to 0
    void ResetBuffers() { std::fill(rdf_buffer.begin(), rdf_buffer.end(), 0.0); measured_steps = 0; }
    /// Returns not normalized RDF bin counts across all timesteps
    std::vector<double>& GetRDFTotalCounts() { return this->rdf_buffer; }
    /// Returns not normalized RDF bin counts of last measurement
    std::vector<double>& GetRDFPrevStepCounts() { return this->cell_processor->get_sample(); }
    /// Returns normalized RDF bin counts across all timesteps
    void GetRDFTotal(std::vector<double>& buffer);
    /// Returns normalized RDF bin counts of last measurement
    void GetRDFPrevStep(std::vector<double>& buffer);
    /// Returns RDF r nodes
    std::vector<double>& GetRNodes() { return this->r_nodes; }
    /// Returns number of measured steps
    double GetMeasuredSteps() const { return measured_steps; }

private:
    /// normalizes the provided rdf buffer
    void normalize(std::vector<double>& buffer, int steps);
    /// number of r nodes
    int number_bins;
    /// distance between two r nodes
    double bin_width;
    /// RDF count buffer across all timesteps
    std::vector<double> rdf_buffer;
    /// RDF r nodes
    std::vector<double> r_nodes;
    /// increases on every call to ProfileData
    int measured_steps;
    /// maximum measure distance
    double max_r;
    /// performs actual rdf measurement
    class InternalCellProcessor : public CellProcessor {
    public:
        InternalCellProcessor& operator=(const InternalCellProcessor&)=delete;
        InternalCellProcessor(double cutoff, double bin_width, int num_bins, double max_r);
        ~InternalCellProcessor() override = default;

        void initTraversal() override;
        void preprocessCell(ParticleCell& cell) override {}
        void processCell(ParticleCell& cell) override;
        void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll=false) override;
        double processSingleMolecule(Molecule* m1, ParticleCell& cell) override {return 0.0;}
        void postprocessCell(ParticleCell& cell) override {}
        void endTraversal() override;

        std::vector<double>& get_sample() { return single_sample_counts; }
    private:
        /// Checks distance and adds to RDF buffer if in range
        void handleMoleculePair(Molecule& m1, Molecule& m2, int thread);

        /// RDF count buffer across a single timestep
        std::vector<double> single_sample_counts;
        /// thread local data of single_sample_counts
        std::vector<std::vector<double>> thread_single_sample_counts;
        /// distance between two r nodes
        double bin_width;
        /// maximum measure distance
        double max_r;
    } *cell_processor;
};