#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

class PotentialProfiler:public PluginBase{
    private:

    CellProcessor* cell_processor;
    int number_bins;
    double bin_width;
    std::vector<int> bin_counts;//total number of COMs in bin[i]
    std::vector<std::vector<double>> pot_values_per_timestep;
    std::vector<double> pot_vals;//U values within the bin[i]
    std::vector<double> pot_counts;
    int sample_frequency;
    int measured_steps;
    double measured_distance_squared;

    public:
    
    PotentialProfiler();
    ~PotentialProfiler(){delete cell_processor;}
    void readXML(XMLfileUnits& xmlconfig) override;
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
    void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{
        std::cout<<"End of step\n";
        for(int i=0;i<pot_vals.size();++i){
            std::cout<<"Value: "<<pot_vals[i]<<"\t"<<pot_counts[i]<<"\n";
        }
        for(int i=0;i<pot_vals.size();++i){
            pot_counts[i]=0;
            pot_vals[i] =0;
        }
    }
    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
    void afterForces(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep) override{   
        if(simstep%sample_frequency ==0 && simstep > global_simulation->getInitStatistics()){
        measured_steps++;
        pc->traverseCells(*cell_processor);    
    }}
    void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override;
    void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{
        WriteDataToFile(particleContainer, domain);
    }
    std::string getPluginName()  {return "PotentialProfiler";}
    static PluginBase* createInstance() {return new PotentialProfiler(); }

    public:

    std::array<double,3> GetCOM(Molecule* m);
    void ProcessDistance(double distance, double pot);

    private:
    void SetBinContainer(ParticleContainer* pc);
    void WriteDataToFile(ParticleContainer* particleContainer, Domain* domain);


};

class CustomCellProcessor: public CellProcessor{

    private: 

    PotentialProfiler* const my_profiler;

    public:
    CustomCellProcessor& operator=(const CustomCellProcessor&)=delete;
    CustomCellProcessor(const double cutoff, PotentialProfiler* r);
    ~CustomCellProcessor(){};

    void initTraversal() override {}
    void preprocessCell(ParticleCell& cell) override {}
    void processCell(ParticleCell& cell) override;
    void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll=false) override;
    double processSingleMolecule(Molecule* m1, ParticleCell& cell) override {return 0.0;}
    void postprocessCell(ParticleCell& cell) override {}
    void endTraversal() override {}

    private:

    double DistanceBetweenCOMs(std::array<double,3>& com1, std::array<double,3>& com2);

    double PotentialCallBack(double eps, double sigma, double r);

};