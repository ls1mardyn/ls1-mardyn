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
#include "Region.h"
#include "Resolution.h"
#include "WeightFunction.h"
#include "Interpolate.h"
#include "plugins/RDFAtCOM.h"
#include "ForceAdapter.h"
#include "ProfilerPMF.h"
#include "IterativeBoltzmannInversion/Convergence.h"
#include "IterativeBoltzmannInversion/ConvergenceCheck.h"
#include "plugins/PotentialOfMeanForce/common.h"
#include "Statistics.h"
#include "common.h"

class InteractionSite;

enum class Mode{Equilibration, EffectivePotential, Production};

class PMF:public PluginBase{

    private:
    Mode mode = Mode::EffectivePotential;
    Interpolate reference_rdf_interpolation;//should not touch
    Interpolate potential_interpolation;//updated on every step or stride
    Interpolate avg_rdf_interpolation;//current avg value so far 
    Interpolate derivative_interpolation;//stores derivatives potential
    InternalProfiler profiler;//measuring tool
    Convergence convergence;
    ResolutionHandlerBase* resolution_handler;
    /**
     * L2 norm, used for hybrid
     */
    WeightFunction weight_function;
    InteractionForceAdapter* pairs_handler;

    StatisticsAdResS adres_statistics;

    int update_stride=1000;//how often to IBI
    int output_freq=1;
    double multiplier;//alpha for step size
    bool output;
    bool do_nothing;//kill switch for all functionalities

    std::string ref_rdf_path;
    std::string effective_potential_path;

    public:
    PMF();
    ~PMF(){}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override;
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{
    };
    void siteWiseForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step) override;
    std::string getPluginName(){ return "PMF";}
    static PluginBase* createInstance() {return new PMF(); }

    /***********************
     * 
     * FUNCTIONS NOT FROM INTERFACE
     * 
     ***********************/
    public:
    std::vector<FPRegion>& GetRegions();
    // ResolutionType GetMoleculeResolution(unsigned long idx);
    // InteractionSite GetMoleculeCOMSite(unsigned long idx);
    double WeightValue(const std::array<double,3>& pos, FPRegion& region);
    Interpolate& GetRDFInterpolation();
    Interpolate& GetAVGRDFInterpolation();
    Interpolate& GetPotentialInterpolation();

    void MapToAtomistic(std::array<double,3> f, Molecule& m1, Molecule& m2);


    double GetMultiplier(){
        return multiplier;
    }

    private:
    /**
     * Read in reference RDF
     */
    void ReadRDF();
    void ReadEffectivePotential();
    /**
     * Implements U(r)_0 = -T^*ln(g(r)_0)
     */
    void SetPotentialInitialGuess();
    /**
     * Implements U(r)_{i+1} =U(r)_i - alpha*T^*ln(g(r)_i/g(r)_*)
     */
    void AddPotentialCorrection(unsigned long step);
    void UpdateRDFInterpolation();



};
