/*
 * Created on Tue Jun 11 2024
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
#include "RDFProfiler.h"
#include "../ExtraConsole.h"

/**
 * Performs IBI on the current phasespace.
 * */
class IBI : public PluginBase {
public:
    //========================================
    // INTERFACE METHODS
    //========================================
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

    std::string getPluginName() { return "IBI"; }
    static PluginBase* createInstance() {return dynamic_cast<PluginBase*>( new IBI()); }

    //========================================
    // NOT USED
    //========================================

    /// Print avg rdf value to file, checks convergence. Then sets up interpolation buffers and accumulates to those.
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override { };
    void siteWiseForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step) override { };
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override { };
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override { };

    //========================================
    // FUNCTIONS NOT FROM INTERFACE
    //========================================
    FunctionPL& GetRDFInterpolation() { return reference_rdf; };

private:
    /// Implements U(r)_0 = -T^*ln(g(r)_0)
    void InitializePotentialValues();
    /// Implements U(r)_{i+1} =U(r)_i + alpha*T^*ln(g(r)_i/g(r)_*)
    void AddPotentialCorrection();
    /// Implements F(r) = - d/dr U(r)
    void DerivativeOfPotential();
    /// Writes the total average RDF
    void WriteRDF();
    /// Creates a new single site component and replaces the comp pointer of all molecules
    void CreateSwapComponent();
    /// Generates file paths
    [[nodiscard]] std::string createFilepath(const std::string& prefix) const;
    /// Creates the optimizer based on the current configuration
    void CreateOptimizer();

    /// g*(r)
    FunctionPL reference_rdf {0.0, 1.0};
    /// V0(r)
    FunctionPL reference_potential {1e+6, 0.0};
    /// measuring tool
    RDFProfiler profiler;
    /// pair handler to comp forces
    IBIPairsHandler* pairs_handler = nullptr;
    /// disables most functionality and only outputs rdf at the end of the simulation
    bool mode_initial_rdf = false;
    /// Target simulation temperature
    double T = 0.0;
    /// path to rdf0 file
    std::string rdf_path = "";
    /// path to active potential file (may not be set, then U0 is used)
    std::string pot_path = "";
    /// IBI iteration counter
    int ibi_iteration = 0;
    /// disables potential update and rdf sampling
    bool mode_equilibrate = true;
    /// number of equilibration steps
    int steps_equilibration = 1e+6;
    /// number of rdf measurement steps
    int steps_measurement = 1e+5;
    /// different execution path -> do all steps in one simulation
    bool mode_single_run = false;
    /// counter for elapsed simulation steps of current IBI phase
    int current_steps = 0;
    /// single run: current ibi phase
    enum { FIRST_INIT, EQUILIBRATE, MEASURE } ibi_phase = FIRST_INIT;
    /// configuration parameters for optimization
    struct OptConfig {
        enum { DEFAULT, ADAM } type;
        double alpha;
        double beta1;
        double beta2;
        double eps;
    } optConfig;

    /**
     * Provides methods for convergence checking
     * */
    struct Convergence {
        void init(double threshold, const std::string& mode_str = "l2", const std::string& stop_str = "worse", int window = 10, int ignore = 0);

        /**
         * Computes integrals over r_i and r_0... see https://arxiv.org/pdf/1410.1853
         * */
        std::pair<bool, double> integral(const FunctionPL& ref, RDFProfiler& profiler);

        /**
         * Computes simple l2 norm over functions
         * */
        std::pair<bool, double> l2(const FunctionPL& ref, RDFProfiler& profiler);

        std::pair<bool, double> operator()(const FunctionPL &ref, RDFProfiler& profiler);

        void LogValues(std::ostream& ostream);

        bool ShouldStop();

    private:
        /// Convergence values of each iteration
        std::vector<double> conv_values;
        /// mode selection
        enum {INTEGRAL, L2} mode = L2;
        /// conv threshold
        double threshold;
        /// early stopping method
        enum {ON_WORSE, WINDOW} stopping_mode = ON_WORSE;
        /// window size for window stopping method
        int window_size;
        /// offset counter, will not return true for should stop as long as this is > 0
        int ignore_counter = 0;

    } ConvergenceCheck;

    /// handles the update step
    std::unique_ptr<IBIOptimizer> optimizer;

    /// external logging
    std::unique_ptr<Terminal::PlotLogger<double>> extLog_conv;
    std::unique_ptr<Terminal::PlotLogger<double>> extLog_force;
    std::unique_ptr<Terminal::PlotLogger<double>> extLog_pot;
    std::unique_ptr<Terminal::PlotLogger<double>> extLog_rdf;
};