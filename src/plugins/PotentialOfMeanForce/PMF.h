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
#include "InteractionCellProcessor.h"
#include "Resolution.h"
#include "WeightFunction.h"
#include "Interpolate.h"
#include "plugins/RDFAtCOM.h"
#include "ForceAdapter.h"
#include "ProfilerPMF.h"

//TODO: who owns the InteractionCellProcessor, AdResS or Simulation?

class InteractionSite;
class PMF:public PluginBase{
    using tracker = std::pair<InteractionSite,ResolutionType>;
    private:



    //std::vector<double> r_nodes;//stores distance values
    //std::vector<double> v_nodes;//stores g(r) values
    Interpolate reference_rdf_interpolation;//should not touch
    Interpolate potential_interpolation;//updated on every step
    Interpolate current_rdf_interpolation;//not needed
    Interpolate avg_rdf_interpolation;//current avg value so far 

    InternalProfiler profiler;//measuring tool

    std::vector<FPRegion> regions;//should create AT region with
    InteractionCellProcessor* adres_cell_processor;
    std::map<unsigned long, tracker> sites;
    ResolutionHandler resolution_handler;
    WeightFunction weight_function;
    RadialDFCOM rdf;//used to generate COM mapping and refresh COMs
    InteractionForceAdapter* pairs_handler;

    double internal_bins;//used for profiler
    int measure_frequency;//used for profiler

    double multiplier;//alpha for step size
    /**
     * Intermediate buffer, is a sum of all rdf values per step
     * Does NOT give an avg
     * Must not be used as data
     * Only through GetAverageRDF method if required
     */
    std::vector<double> accumulate_rdf_buffer;
    

    public:
    PMF();
    ~PMF(){}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    /**
     * Updates resolution and tracker of each molecule
     */
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}


    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

    /**
     * 
     * Print avg rdf value to file
     * Convergence check
     * Set up interpolation buffers
     * Accumulate to internal buffer
     */
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override;
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{};
    void siteWiseForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step) override;

    std::string getPluginName(){ return "PMF";}
    static PluginBase* createInstance() {return new PMF(); }

    public:

    std::vector<FPRegion>& GetRegions();
    ResolutionType GetMoleculeResolution(unsigned long idx);
    InteractionSite GetMoleculeCOMSite(unsigned long idx);
    double WeightValue(const std::array<double,3>& pos, FPRegion& region);
    void ReadRDF();
    Interpolate& GetRDFInterpolation();
    Interpolate& GetAVGRDFInterpolation();
    Interpolate& GetPotentialInterpolation();
    Interpolate& GetCurrentRDFInterpolation();

    /**
     * Implements U(r)_0 = -T^*ln(g(r)_0)
     */
    void InitializePotentialValues();
    void AddPotentialCorrection();
    void AccumulateRDF(std::vector<double>& current_rdf);
    /**
     * 
     * Returns accumulated RDF divided by measured steps from profiler
     * Returns a full copy vector
     * 
     */
    std::vector<double> GetAverageRDF();
    double ConvergenceCheck();
    double GetMultiplier(){
        return multiplier;
    }

    public: 
    /**
     * Maps the com forces to FP, formula:
     * F_{i\alpha\beta} = \num{m_{i\alpha}}\den{\sum_{i\alpha}{m_{i\alpha}}}F^{cm}_{\alpha\beta}
     * 
     * Receives f already weighted
     * Same callbacks as in potforce
     * Assume both m1 and m2 have the same component class...
     * Comply with newton 3rd law
     */
    void MapToAtomistic(std::array<double,3> f, Molecule& m1, Molecule& m2);


};

/**
 * Stores velocity, position, force, potential
 */
class InteractionSite:public Site{
    private:
    double u_com;
    std::array<double,3> f_com;
    std::array<double,3> v_com;//not used

    public: 
    void SubForce(std::array<double, 3> f);
    void AddForce(std::array<double,3> f);
    void AddPotential(double pot);
    void SetPosition(std::array<double,3> pos);
    void SetVelocity(std::array<double,3> pos);

    std::array<double,3>& GetPosition();
    std::array<double,3>& GetVelocity();


};