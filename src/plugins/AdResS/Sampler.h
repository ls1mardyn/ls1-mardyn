/*
 * Created on Wed Apr 24 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

#include "Grid.h"
#include "GridHandler.h"



class Sampler{
    
    private:
    //TODO:split into node and element sample types
    struct Samples{
        std::vector<int> particles_per_cell;
        std::vector<double> material_density;
        std::vector<int> particles_per_node;
        double target_number_density;
    };

    Samples samples;
    Grid* grid;
    double measure_radius;
    GridHandler handler;

    public:

    Sampler();
    void init(Grid* grid);
    GridHandler& GetGridHandler();
    void SetMeasureRadius(double r);
    void SampleAtNodes(ParticleContainer* pc);

    //
    void SetSubsetMaterialDensityValues();
    void SetTargetValue();
    Samples& GetSamples(){
        return samples;
    }
    template<class T>
    std::ostream& WritePlaneSamples(std::ostream& out, T data);
    //TODO: really has to be templatized or so, code repetition
    std::ostream& WriteSample(std::ostream& out, std::vector<double>& smpl);
    std::ostream& WriteSample(std::ostream& out, std::vector<int>& smpl);
    std::ostream& WriteInfo(std::ostream& out);
    private:
    IdxArray GetParticleLocalCellIndices(ParticleIterator it);
    bool ParticleInsideMeasuringSpace(std::array<double, 3>, std::array<double, 3> par_pos);
    
};


template<class T>
std::ostream& Sampler::WritePlaneSamples(std::ostream& out, T data){
    std::string prefix = "//[Sampler]: ";
    std::vector<int> plane_nodes = handler.GetNodesOnPlane(grid,10);
    out<<prefix+" target value: "<<this->samples.target_number_density<<"\n";
    for(int i=0;i<plane_nodes.size();i++){
        double rx = grid->GetNodes()[plane_nodes[i]].GetPosition()[0];
        double ry = grid->GetNodes()[plane_nodes[i]].GetPosition()[1];
        double val = data[plane_nodes[i]];
        out<<plane_nodes[i]<<"\t"<<rx<<"\t"<<ry<<"\t"<<val<<"\n";
    }
    return out;
}

template<typename T>
class TimeAveraging{
    public:
    TimeAveraging():step_count(0){
        
    }
    void SetDataSize(T data){
        averaged_data.resize(data.size());
        std::fill(averaged_data.begin(),averaged_data.end(),0.0);
    }
    void AverageData(T data){
        step_count++;
        if(data.size() != averaged_data.size()){
            std::cout<<"[TimeAveraging] data structure mismatch" << std::endl;
            averaged_data.resize(data.size());
            std::fill(averaged_data.begin(),averaged_data.end(),0.0);
        }
        for(int i=0;i<data.size();i++){
            averaged_data[i] = averaged_data[i]+data[i];
        }

    }

    T& GetAveragedData(){
        return this->averaged_data;
    }

    T GetAveragedDataCopy(){
        T copy=averaged_data;
        for(int i=0;i<copy.size();i++){
            copy[i] = copy[i]/(double)step_count;
        }
        return copy;
    }

    int GetStepCount(){
        return this->step_count;
    }

    std::ostream& WriteAverage(std::ostream& out, T data){
        std::string prefix ="//[TimeAverage]: ";

        out<<prefix+"data average after: "<<step_count<<" steps"<<"\n";
        out<<prefix+"data structure with size: "<<data.size()<<"\n";
        for(int i=0;i<data.size();i++){
            out<<i<<"\t"<<data[i]<<"\t"<<data[i]/(double)step_count<<"\n";
        }
        return out;
    }



    private:
    T averaged_data;
    int step_count;
};