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

/**
 * Performs IBI on the current phasespace.
 * */
class PMF : public PluginBase {
public:
    //========================================
    // INTERFACE METHODS
    //========================================
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

    std::string getPluginName() { return "PMF"; }
    static PluginBase* createInstance() {return dynamic_cast<PluginBase*>( new PMF()); }

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
    double ConvergenceCheck();

private:
    /// Implements U(r)_0 = -T^*ln(g(r)_0)
    void InitializePotentialValues();
    /// Implements U(r)_{i+1} =U(r)_i - alpha*T^*ln(g(r)_i/g(r)_*)
    void AddPotentialCorrection();
    /// Implements F(r) = - d/dr U(r)
    void DerivativeOfPotential();
    /// Writes the total average RDF
    void WriteRDF();
    /// Generates file paths
    [[nodiscard]] std::string createFilepath(const std::string& prefix) const;

    /// g*(r)
    FunctionPL reference_rdf {0.0, 1.0};
    /// V0(r)
    FunctionPL reference_potential {1e+6, 0.0};
    /// measuring tool
    RDFProfiler profiler;
    /// pair handler to comp forces
    IBIPairsHandler* pairs_handler = nullptr;
    /// alpha for step size
    double alpha = 0.0;
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
    /// collection of convergence values
    std::vector<double> convergence_steps;
};