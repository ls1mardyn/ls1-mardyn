/*
 * Created on Wed May 08 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#include "Sampler.h"

GridSampler2::GridSampler2(Grid* g, double rad):grid(g),measure_radius(rad){

}

void GridSampler2::SampleData(ParticleContainer* pc){
    this->SampleAtNodes(pc);
}

void GridSampler2::init(){
    int total_nodes = this->grid->GetNodeInfo().total_nodes;
    std::vector<double>& vec = this->GetSampledData();
    vec.resize(total_nodes);
    std::fill(vec.begin(),vec.end(),0.0);
}

void GridSampler2::writeSample(const std::string& filename){
    std::ofstream file(filename);
    file <<"#Total nodes sampled: "<<GetSampledData().size()<<"\n";
    std::vector<Grid3D::Node>& nodes = grid->GetNodes();

    for(int i=0;i<GetSampledData().size();i++){
        std::array<double, 3>& pos = nodes[i].GetPosition();
        file<<i<<"\t"<<pos[0]<<"\t"<<pos[1]<<"\t"<<pos[2]<<"\t"<<GetSampledData()[i]<<"\n";
    }

}

void GridSampler2::SampleAtNodes(ParticleContainer* pc){

    ParticleIterator it = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    std::fill(GetSampledData().begin(), GetSampledData().end(), 0);
    std::vector<Node>& all_nodes = grid->GetNodes();
    for(it;it.isValid();++it){
        std::array<double, 3> part_pos = it->r_arr();
        //Check for every node, if a particle is inside of its "bubble"
        for(int nidx=0;nidx<GetSampledData().size();nidx++){
            
            auto nodal_pos=all_nodes[nidx].GetPosition();
            //check euclidean distance between particle as center of sphere/bubble
            if(ParticleInsideMeasuringSpace(nodal_pos, part_pos)){
                GetSampledData()[nidx] +=1;
            }
        }
    }
    //Convert to number density values
    double sphere_volume = 4.0/3.0 * M_PI * measure_radius*measure_radius*measure_radius;
    for(int nidx=0;nidx<GetSampledData().size();nidx++){
        
        GetSampledData()[nidx] = (double)GetSampledData()[nidx]/sphere_volume;
    }
}


bool GridSampler2::ParticleInsideMeasuringSpace(std::array<double, 3> nodal_pos, std::array<double, 3> par_pos){
    bool is_inside=false;

    std::array<double, 3> distance;

    distance[0] = par_pos[0]-nodal_pos[0];
    distance[1] = par_pos[1]-nodal_pos[1];
    distance[2] = par_pos[2]-nodal_pos[2];

    double norm = std::sqrt(distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2]);

    if(norm<measure_radius){
        is_inside=true;
    }
    return is_inside;

 }

 /************
  * 
 */

AveragedGridSampler::AveragedGridSampler(Grid* grid, double rad,int freq):GridSampler2{grid,rad},write_frequency(freq){

}

void AveragedGridSampler::init(){
    GridSampler2::init();
    averager.SetDataSize(GridSampler2::GetSampledData());
}

void AveragedGridSampler::SampleData(ParticleContainer* pc){
    GridSampler2::SampleData(pc);
    GridSampler2::SampleData(pc);
    averager.AverageData(GridSampler2::GetSampledData());
}

void AveragedGridSampler::writeSample(const std::string& filename){
    GridSampler2::writeSample(filename);
    
    //If averager step_count matches desired frequency
    if(averager.GetStepCount()%this->write_frequency==0){
        std::string name{"averaged"+filename};
        std::ofstream averagedfile(name);
        averager.WriteAverage(averagedfile,averager.GetSumData());
        averagedfile.close();
    }
}