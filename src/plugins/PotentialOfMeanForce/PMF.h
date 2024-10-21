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
#include "Filter.h"
#include "Convergence.h"

class InteractionSite;

class PMF:public PluginBase{

    // using tracker = std::pair<InteractionSite,ResolutionType>;
    private:

    Interpolate reference_rdf_interpolation;//should not touch
    Interpolate potential_interpolation;//updated on every step or stride
    Interpolate acc_rdf_interpolation;//current avg value so far 
    Interpolate derivative_interpolation;//stores derivatives potential

    InternalProfiler profiler;//measuring tool
    Filter filter;
    Convergence convergence;
    std::vector<FPRegion> regions;
    /**
     * Stores a COM site for every molecule (needs improvement)
     */
    // std::map<unsigned long, tracker> sites;
    ResolutionHandler resolution_handler;
    ResolutionComponentHandler res_comp_handler;
    /**
     * L2 norm, used for hybrid
     */
    WeightFunction weight_function;
    InteractionForceAdapter* pairs_handler;

    double internal_bins;//used for profiler
    int measure_frequency;//used for profiler
    int update_stride=300;//how often to IBI
    double multiplier;//alpha for step size
    bool output;


    // DataPostProcess post_processing;

    struct ConvergenceTest{
        bool ConvergenceCheck(std::vector<double>& ref, std::vector<double>& crrnt){
            double conv=0.0;
            for(int i=0;i<ref.size();++i){
                conv += std::abs(ref[i]-crrnt[i]);
            }
            convergence_per_step.emplace_back(conv);
            ++ibi_iteration;
            if(conv<tolerance)
            return true;

            return false;
        }
        double tolerance=0.05;
        std::vector<double> convergence_per_step;
        int ibi_iteration=0;
    };

    ConvergenceTest convergence_check;

    public:
    PMF();
    ~PMF(){}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    /**
     * Updates resolution and tracker of each molecule
     */
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
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
    ResolutionComponentHandler& GetResolutionHandler(){
        return res_comp_handler;
    }

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

    /**
     * 
     * Returns accumulated RDF divided by measured steps from profiler
     * Returns a full copy vector
     * 
     */
    std::vector<double> GetAverageRDF();
    double GetMultiplier(){
        return multiplier;
    }

    private:
    /**
     * Read in reference RDF
     */
    void ReadRDF();
    /**
     * Implements U(r)_0 = -T^*ln(g(r)_0)
     */
    void InitializePotentialValues();
    /**
     * Implements U(r)_{i+1} =U(r)_i - alpha*T^*ln(g(r)_i/g(r)_*)
     */
    void AddPotentialCorrection(unsigned long step);
    void AccumulateRDF(ParticleContainer* pc, Domain* domain);



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