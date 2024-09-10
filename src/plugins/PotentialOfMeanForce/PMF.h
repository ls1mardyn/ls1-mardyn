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
    Interpolate reference_rdf_interpolation;
    Interpolate potential_interpolation;
    Interpolate current_rdf_interpolation;
    std::vector<FPRegion> regions;//should create AT region with
    InteractionCellProcessor* adres_cell_processor;
    std::map<unsigned long, tracker> sites;
    ResolutionHandler resolution_handler;
    WeightFunction weight_function;
    RadialDFCOM rdf;
    InteractionForceAdapter* pairs_handler;
    InternalProfiler profiler;
    double multiplier;
    

    public:
    PMF();
    ~PMF(){}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

    /**
     * Updates resolution and tracker of each molecule
     */
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{
        //pass buffers to interpolations
        current_rdf_interpolation.SetYValues(profiler.GetRDFValues());
        potential_interpolation.SetYValues(profiler.GetPotentialValues());
        //transfer buffers to interpolation
        profiler.ResetBuffers();
    }


    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{
        this->profiler.ProfileData(particleContainer,simstep);
    }
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override{
        this->profiler.GenerateInstantaneousData(particleContainer,domain);
        //convergence check??
        if(global_simulation->getSimulationStep()>0){
            Log::global_log->info()<<"[PMF] Convergence check: "<<ConvergenceCheck()<<std::endl;
            std::string filename="rdf_"+std::to_string(simstep);
            std::ofstream rdf_file(filename);
            for(int i=0;i<current_rdf_interpolation.GetGValues().size();++i){
                rdf_file<<std::setw(8)<<std::left<<current_rdf_interpolation.GetRValues()[i]<<"\t"<<std::setw(8)<<std::left<<current_rdf_interpolation.GetGValues()[i]<<std::endl;
            }
        }
    }
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
    Interpolate& GetPotentialInterpolation();
    Interpolate& GetCurrentRDFInterpolation();
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