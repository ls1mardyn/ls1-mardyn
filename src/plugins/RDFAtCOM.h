/*
 * Created on Thu Jun 06 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

class RadialDFCOM:public PluginBase{

    private:

    CellProcessor* cell_processor;
    int number_bins;
    double bin_width;
    std::vector<double> pairs_per_bin;//total number of COMs in bin[i]
    std::vector<std::vector<double>> average_pairs_per_bin;//pairs per measured time step
    std::vector<double> accumulated_com_rdf;//accumulated over all measured time steps
    
    int sample_frequency;
    int measured_steps;
    int average_frequency;//controls when the stored vectors are gonna be used
    double measured_distance_squared;

    public:
    RadialDFCOM();
    ~RadialDFCOM(){delete cell_processor;}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void afterForces(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep) override{   
        if(simstep%sample_frequency ==0 && simstep > global_simulation->getInitStatistics()){
        measured_steps++;
        pc->traverseCells(*cell_processor);    
    }}
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override{
        //Generate && Print instantaneous data
        if(simstep%sample_frequency ==0 && 
        simstep > global_simulation->getInitStatistics()
        ){

            average_pairs_per_bin.emplace_back(pairs_per_bin);
            WriteRDFToFile(particleContainer, domain, simstep);

            //reset instantaneous buffer after printout and end of step
            std::fill(pairs_per_bin.begin(),pairs_per_bin.end(),0.0);
        }   
    }
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{
        WriteRDFToFile(particleContainer, domain, global_simulation->getSimulationStep());
    }
    std::string getPluginName()  {return "RadialDFCOM";}
    static PluginBase* createInstance() {return new RadialDFCOM(); }

    public:

    std::array<double,3> GetCOM(Molecule* m);
    void ProcessDistance(double distance);

    private:
    void SetBinContainer(ParticleContainer* pc);
    void WriteRDFToFile(ParticleContainer* particleContainer, Domain* domain, unsigned long simstep);
    

};

/**
 * Based on RDF cell processor
 * 
 */
class COMDistanceCellProcessor: public CellProcessor{

    private: 

    RadialDFCOM* const rdf;

    public:
    COMDistanceCellProcessor& operator=(const COMDistanceCellProcessor&)=delete;
    COMDistanceCellProcessor(const double cutoff, RadialDFCOM* r);
    ~COMDistanceCellProcessor(){};

    void initTraversal() override {}
    void preprocessCell(ParticleCell& cell) override {}
    void processCell(ParticleCell& cell) override;
    void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll=false) override;
    double processSingleMolecule(Molecule* m1, ParticleCell& cell) override {return 0.0;}
    void postprocessCell(ParticleCell& cell) override {}
    void endTraversal() override {}

    private:

    double DistanceBetweenCOMs(std::array<double,3>& com1, std::array<double,3>& com2);

};