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

/**
 * Calculate radial distribution function on multi-site molecule center of mass.
 * Uses R = 
 */

class RadialDFCOM:public PluginBase{

    private:

    CellProcessor* cell_processor;
    int number_bins;
    double bin_width;
    std::vector<int> bin_counts;//total number of COMs in bin[i]
    int sample_frequency;
    int measured_steps;
    public:
    RadialDFCOM();
    ~RadialDFCOM(){}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override;
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{
        std::cout<<bin_counts[0]<<"\n";
        WriteRDFToFile(particleContainer, domain);
    }
    std::string getPluginName()  {return "RadialDFCOM";}
    static PluginBase* createInstance() {return new RadialDFCOM(); }

    public:

    std::array<double,3> GetCOM(Molecule* m);
    void ProcessDistance(double distance);

    private:
    void SetBinContainer(ParticleContainer* pc);
    void WriteRDFToFile(ParticleContainer* particleContainer, Domain* domain);
    

};

/**
 * Based on RDF cell processor
 * 
 */
class COMDistanceCellProcessor: public CellProcessor{

    private: 

    RadialDFCOM* rdf;

    public:

    COMDistanceCellProcessor(RadialDFCOM* r);
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