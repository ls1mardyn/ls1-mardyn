#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

class PotentialProfiler : public PluginBase {
public:
    PotentialProfiler() = default;
    ~PotentialProfiler() override { delete cell_processor; }
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void afterForces(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override;

    //==========================
    // Unused methods
    //==========================
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {};
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override {};
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override {};

    //==========================
    // Util methods
    //==========================
    std::string getPluginName()  {return "PotentialProfiler";}
    static PluginBase* createInstance() { return new PotentialProfiler(); }

private:
    void WriteDataToFile(unsigned long simstep);

    std::vector<double> current_potential_average {};
    int sample_frequency = 0;
    int avg_frequency = 0;
    int measured_steps = 0;
    int number_bins = 0;
    double bin_width = 0;

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

        std::vector<double>& get_sample_counts() { return single_sample_counts; }
        std::vector<double>& get_sample_pots() { return single_sample_pots; }
    private:
        /// Checks distance and adds to RDF buffer if in range
        void handleMoleculePair(Molecule& m1, Molecule& m2, int thread);
        /// Calculates potential with given parameters
        double calcPot(double sigma, double epsilon, double r2);

        /// count buffer across a single timestep
        std::vector<double> single_sample_counts;
        /// pot buffer across a single timestep
        std::vector<double> single_sample_pots;
        /// thread local data of single_sample_counts
        std::vector<std::vector<double>> thread_single_sample_counts;
        /// thread local data of single_sample_pots
        std::vector<std::vector<double>> thread_single_sample_pots;
        /// distance between two r nodes
        double bin_width;
        /// maximum measure distance
        double max_r;
    } *cell_processor = nullptr;
};